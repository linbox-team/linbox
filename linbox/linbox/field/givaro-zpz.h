/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/givaro-zpz.h
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

/* WARNING this wrapper works only with an improved version of Givaro.
 * This version of givaro won't be available for public yet.
 * But it is available on my web page.
 * You can send me a mail to get it or for others details.
 */

#ifndef __FIELD_GIVARO_ZPZ
#define __FIELD_GIVARO_ZPZ


#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
#include "linbox/vector/vector-domain.h"
//-------------------------------------
// Files of Givaro library
#include <givzpz16std.h>
#include <givzpz32std.h>
#include <givzpz16table1.h>
#include <giv_randiter.h>

//--------------------------------------



// Namespace in which all LinBox code resides
namespace LinBox 
{ 

	/*  This wrappers allows to use three sorts of givaro fields :
	 *  Elements represent by a 32 bits integer
	 *  Elements represent by a 16 bits integer
	 *  Elements represent in Zech log representation by a 16 bits integer
	 *
	 *  To use this fields with the wrapper below just replace the template
	 *  parameter by the Tag appropriated.
	 *  "Std16"  for 16 bits integer
	 *  "Std32"  for 32 bits integer
	 *  "Log16"  for Zech log representation in 16 bits
	 */




	/** This template class is define just to be in phase with the LinBox
	 *  archetype. Read the archetype to know all functions are available.
	 *  Most of all methods are inherited from ZpzDom<Std16>, ZpzDom<Std32>
	 *  and ZpzDom<log16> class of Givaro.
	 *  these class allow to construct only finite field with a prime modulus.
	 */   

	template <class TAG> class GivaroZpz : public ZpzDom<TAG>, public FieldInterface
	{

	private:

		friend class DotProductDomain<GivaroZpz<TAG> > ;
				
	public:

		//typedef integer Integer;

		/** Element type.
		 *  This type is inherited from the Givaro class ZpzDom<TAG>
		 */
		typedef typename ZpzDom<TAG>::Rep Element;

		/** RandIter type
		 *  This type is inherited from the Givaro class ZpzDom<TAG>
		 */	
		typedef GIV_randIter< ZpzDom<TAG>, integer > RandIter;

		/** Constructor from an integer
		 *  this constructor use the ZpzDom<TAG> constructor
		 */
		GivaroZpz (const integer &p) : ZpzDom<TAG> (static_cast<typename TAG::type> (p))  {}

		/** Copy constructor
		 * This copy constructor use the ZpzDom<TAG> copy constructor
		 */
		GivaroZpz (const GivaroZpz<TAG>& F) : ZpzDom<TAG> (F) {}

		/** Operator =
		 */
		GivaroZpz<TAG>& operator= (const GivaroZpz<TAG>& F) {
			return (*this)=F;
		}

		/** Characteristic.
		 * Return integer representing characteristic of the domain.     
		 * @return integer representing characteristic of the domain.
		 */
		integer &characteristic (integer &c) const
			{ return c = integer (static_cast<int> (ZpzDom<TAG>::size ())); }

		/** Cardinality. 
		 * Return integer representing cardinality of the domain.     
		 * @return integer representing cardinality of the domain
		 */
		integer &cardinality (integer &c) const
			{ return c = integer (static_cast<int> (ZpzDom<TAG>::size ())); }

		/** Conversion of field base element to an integer.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to an integer.
		 * @param x integer to contain output (reference returned).
		 * @param y constant field base element.
		 */
		integer &convert (integer &x, const Element &y) const
			{ return x = integer (static_cast<int> (y)); }

		/** Initialization of field base element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.     
		 * @return reference to field base element.
		 * @param x field base element to contain output (reference returned).
		 * @param y integer.
		 */  
		Element &init (Element &x , const integer &y = 0) const
			{ return ZpzDom<TAG>::init (x, long (y)); }

	}; // class GivaroZpz<TAG>
 
	/** Specialisation of the convert function for the zech log representation
	 *	of givaro-zpz (GivaroZpz<Log16>.
	 *  this function translates the internal representation to the real 
	 *	value of the element.
	 *	This can have no sense but can be usefull
	 *  NB : the init function for this specialisation does the same thing.
	 *  the function transaltes the values to her internal representation.
	 */ 
	integer& GivaroZpz<Log16>::convert(integer& x, const Element& y) const
	{     
		int tmp = _tab_rep2value[y];
		return x = integer (tmp);
	}



	// Specialization of DotProductDomain for GivaroZpz<Std32> field

	template <>
	class DotProductDomain<GivaroZpz<Std32> > : private virtual VectorDomainBase<GivaroZpz<Std32> >
	{
	public:
		
		typedef GivaroZpz<Std32>::Element Element;
		
		DotProductDomain (const GivaroZpz<Std32> &F)
			: VectorDomainBase<GivaroZpz<Std32> > (F) ,
			  Corr(uint64(-1) % (uint64)F._p +1),
			  Max(uint64(-1))
		{}
		
	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;
		
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
		
	private:
		uint64 Corr;
		uint64 Max;
	};

	// Specialization of DotProductDomain for GivaroZpz<Std16> field

	template <>
	class DotProductDomain<GivaroZpz<Std16> > : private virtual VectorDomainBase<GivaroZpz<Std16> >
	{
	public:

		typedef GivaroZpz<Std16>::Element Element;

		DotProductDomain (const GivaroZpz<Std16> &F)
			: VectorDomainBase<GivaroZpz<Std16> > (F) ,
			  Corr(uint32(-1) % (uint32)F._p +1),
			  Max(uint32(-1))
		{}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
	
	private:
		uint32 Corr;
		uint32 Max;
};



	/* Specialization of FieldAXPY for GivaroZpz<Std32> Field */

	template <>
	class FieldAXPY<GivaroZpz<Std32> >
	{
	    public:

		typedef GivaroZpz<Std32>::Element Element;
		typedef GivaroZpz<Std32> Field;

		FieldAXPY (const Field &F) : _F (F) , Corr(uint64(-1) % (uint64)F._p +1){ _y = 0; }
		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F), _y (0) , Corr(faxpy.Corr) {}

		FieldAXPY<GivaroZpz<Std32>> &operator = (const FieldAXPY &faxpy) 
			{ _F = faxpy._F; _y = faxpy._y; Corr = faxpy.Corr; return *this; }

		inline void accumulate (const Element &a, const Element &x)
		{
			uint64 t = (uint64) a * (uint64) x;
			_y += t;

			if (_y < t)
				_y += Corr;
		}

		inline Element &get (Element &y) {
			_y %= (uint64) _F._p;
			if ((int64) _y < 0) _y += _F._p;
			y = (uint32) _y;
			return y;
		}

		inline FieldAXPY &assign (const Element y)
			{ _y = y; return *this; }

	    private:

		Field _F;
		uint64 _y;
		uint64 Corr;
	};



	/* Specialization of FieldAXPY for GivaroZpz<Std32> Field */
	
	template <>
	class FieldAXPY<GivaroZpz<Std16> >
	{
	    public:

		typedef GivaroZpz<Std16>::Element Element;
		typedef GivaroZpz<Std16> Field;

		FieldAXPY (const Field &F) : _F (F) , Corr(uint32(-1) % (uint32)F._p +1){ _y = 0; }
		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F), _y (0) , Corr(faxpy.Corr) {}

		FieldAXPY<GivaroZpz<Std16>> &operator = (const FieldAXPY &faxpy) 
			{ _F = faxpy._F; _y = faxpy._y; Corr = faxpy.Corr; return *this; }

		inline void accumulate (const Element &a, const Element &x)
		{
			uint32 t = (uint32) a * (uint32) x;
			_y += t;

			if (_y < t)
				_y += Corr;
		}

		inline Element &get (Element &y) {
			_y %= (uint32) _F._p;
			if ((int32) _y < 0) _y += _F._p;
			y = (uint16) _y;
			return y;
		}

		inline FieldAXPY &assign (const Element y)
			{ _y = y; return *this; }

	    private:

		Field _F;
		uint32 _y;
		uint32 Corr;
	};

	
} // namespace LinBox

#include "linbox/field/givaro-zpz.inl"

#endif // __FIELD_GIVARO_ZPZ
