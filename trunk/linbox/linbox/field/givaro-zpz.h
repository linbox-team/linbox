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


#include "linbox-config.h"
#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
//-------------------------------------
// Files of Givaro library
#include <givaro/givzpz16std.h>
#include <givaro/givzpz32std.h>
#include <givaro/givzpz16table1.h>
#include <givaro/giv_randiter.h>

//--------------------------------------

//--------------------------------------
// __LINBOX_XML I/O features
#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <iostream>
#include <string>
using std::istream;
using std::ostream;
using std::cout;
using std::endl;

#endif


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
        template<class Field>
                class DotProductDomain;
        template<class Field>
                class FieldAXPY;

	/** This template class is define just to be in phase with the LinBox
	 *  archetype. Read the archetype to know all functions are available.
	 *  Most of all methods are inherited from ZpzDom<Std16>, ZpzDom<Std32>
	 *  and ZpzDom<log16> class of Givaro.
	 *  these class allow to construct only finite field with a prime modulus.
	 */   

	template <class TAG> class GivaroZpz : public ZpzDom<TAG>, public FieldInterface
	{

	private:

		/*		friend class DotProductDomain<GivaroZpz<TAG> > ;
				friend class FieldAXPY<GivaroZpz<TAG> >; */
				
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

#ifdef __LINBOX_XMLENABLED
		// XML Reader constructor
		GivaroZpz(Reader &R) : ZpzDom<TAG>()
		{
			integer m;
			long e;

			if(!R.expectTagName("field") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagName("finite") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagName("characteristic") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagNum(m)) return;
			
			GivaroZpz<TAG> oth(m);
			*this = oth;


			R.upToParent();
			R.upToParent();

			if(R.getNextChild()) { // check the extension
				if(!R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagName("extension")) return;
				R.traverseChild();
				if(!R.expectTagNum(e)) return;

				if( e > 1) {
					R.setErrorString("Recieved finite field with extension greater than 1.");
					R.setErrorCode(Reader::OTHER);
					return;
				}
				R.upToParent();
				R.upToParent();
				R.getPrevChild();
			}
			R.upToParent();

			return; // field initalization complete

		}
#endif

		/** Copy constructor
		 * This copy constructor use the ZpzDom<TAG> copy constructor
		 */
		GivaroZpz (const GivaroZpz<TAG>& F) : ZpzDom<TAG> (F) {}


		// Rich Seagraves 7-16-2003
		// As is, this operator is an infinite loop
		// By not providing an operator= in GivaroZpz, 
		// the operator= in the base class (ZpzDom<TAG>) is called
		// automatically by the rules of C++, which I'm guessing is
		// the "Right Thing" for this operator
		//
		/** Operator =
		 */
		/*
		       	GivaroZpz<TAG>& operator= (const GivaroZpz<TAG>& F) {
				return (*this)=F;
			}
		*/

		/** Characteristic.
		 * Return integer representing characteristic of the domain.     
		 * @return integer representing characteristic of the domain.
		 */
		integer &characteristic (integer &c) const
			{ return c = integer (static_cast<int> (ZpzDom<TAG>::size ())); }
		long characteristic() const
			{return static_cast<int>(ZpzDom<TAG>::size());}

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
		
		double &convert (double& x, const Element& y) const
		{ return x = (double) y; }
		

		/** Initialization of field base element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.     
		 * @return reference to field base element.
		 * @param x field base element to contain output (reference returned).
		 * @param y integer.
		 */  
		Element &init (Element &x , const integer &y = 0) const
		{ return ZpzDom<TAG>::init (x, (long) (y% integer(_p))); }


		Element &init (Element &x , const double &y ) const
		{ 
			double charact = (double) _p; 
			double tmp=y;
			int sign=0;
			if (tmp < 0.0) {
				tmp=-tmp;
				sign=1;
			}	
					
			if (tmp >= charact) 
				tmp -= (charact * (double) floor( tmp/charact));
						
			if ((sign)&&(tmp))
				tmp= _p - tmp;
		
			return  x=Element(tmp);
		}
			
#ifdef __LINBOX_XMLENABLED
		ostream &write(ostream &os) const
		{
			Writer W;
			if( toTag(W))
				W.write(os);

			return os;
		}

		// this method, which is template specialized, appears at
		// the end of givaro-zpz.inl. . .
		bool toTag(Writer &) const;

		ostream &write(ostream &os, const Element &e) const
		{
			Writer W;
			if( toTag(W, e))
				W.write(os);

			return os;

		}

		bool toTag(Writer &W, const Element &e) const
		{
			string s;
			W.setTagName("cn");
			Writer::numToString(s, e);
			W.addDataChild(s);
			return true;
		}


		istream &read(istream &is, Element &e) const
		{
			Reader R(is);
			if( !fromTag(R, e)) {
				is.setstate(istream::failbit);
				if(!R.initalized())
					is.setstate(istream::badbit);
			}

			return is;
		}

		bool fromTag(Reader &R, Element &e) const
		{
			if(!R.expectTagName("cn") || !R.expectChildTextNum(e))
				return false;
			else {
				ZpzDom<TAG>::assign(e, e);
				return true;
			}
		}

#endif


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

	double& GivaroZpz<Log16>::convert(double& x, const Element& y) const
	{
		int tmp = _tab_rep2value[y];
		return x = (double) tmp;
	}

#ifdef __LINBOX_XMLENABLED
			// XML Reader constructor
	GivaroZpz<Log16>::GivaroZpz(Reader &R) : ZpzDom<Log16>()
	{
		integer m;
		long e, j;
		size_t i;
		
		if(!R.expectTagName("field") || !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("finite") || !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("characteristic") || !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagNum(m)) return;
		
		GivaroZpz<Log16> oth(m);
		// because ZpzDom<Log16> doesn't make a true copy, but instead
		// just does pointer copy, this constructor manually 
		// performs the copy
		_p = oth.characteristic() ;
		_pmone = _p - 1;
		_tab_value2rep = new Power_t[_p];
		_tab_rep2value = new Residu_t[_p];
		_tab_mul = new Power_t[4 *_p];
		
		Power_t* tab_pone = new Power_t[4 * _p];
		Power_t* tab_mone = new Power_t[4 * _p];

		_tab_addone = &tab_pone[2 * _pmone];
		_tab_subone = &tab_mone[2 * _pmone];

		for(i = 0; i < (size_t) _p; ++i) {
			_tab_value2rep[i] = oth._tab_value2rep[i];
			_tab_rep2value[i] = oth._tab_rep2value[i];
		}

		for(i = 0; i < (size_t) 4 * _p; ++i) {
			_tab_mul[i] = oth._tab_mul[i];
		}

		_tab_div = &_tab_mul[_pmone];
		_tab_neg = &_tab_neg[_pmone/2];

		// Sections of code lifted right out of 
		// src/kernel/zpz/givzpz16table1.C in Givaro source
		// distribution
		for(j = 0; j < _pmone; ++j)
			_tab_addone[i] = _tab_value2rep[1 + _tab_rep2value[j] ];
		for(j = 1 - _pmone; j < 0; ++j)
			_tab_addone[j] = _tab_value2rep[1 + _tab_rep2value[j + _pmone]];
		for(j = _pmone; j <= 2 * _pmone; ++j)
			_tab_addone[j] = 0;

		for(j = -2 * _pmone; j < 1 - _pmone; ++j)
			_tab_addone[j] = j;

		for(j=_pmone; j<=(int32)2*_pmone; j++)
			_tab_subone[j] = 0;
		for(j=-2*_pmone; j<(int32)(1-3*_pmone/2); j++)
			_tab_subone[j] = j+_pmone/2;
		for(j=-3*_pmone/2; j<(1-_pmone); j++)
			_tab_subone[j] = j-_pmone/2;
		for(j=1-_pmone; j<(1-_pmone/2); j++)
			_tab_subone[j] = _tab_addone[j + _pmone/2 + _pmone];
		for(j=_pmone/2; j<_pmone; j++)
			_tab_subone[j] = _tab_addone[j - _pmone/2];
		for(j=-_pmone/2; j<_pmone/2; j++)
			_tab_subone[j] = _tab_addone[j+_pmone/2];


		// initalize the field
		//
		

		R.upToParent();
		R.upToParent();

		if(R.getNextChild()) { // check the extension
			if(!R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagName("extension")) return;
			R.traverseChild();
			if(!R.expectTagNum(e)) return;
			
			if( e > 1) {
				R.setErrorString("Recieved finite field with extension greater than 1.");
				R.setErrorCode(Reader::OTHER);
				return;
			}
			R.upToParent();
			R.upToParent();
			R.getPrevChild();
		}
		R.upToParent();
		
		return; // field initalization complete
		
	}

#endif




	/* Specialization of FieldAXPY for GivaroZpz<Std32> Field */

	
	template <>
	class FieldAXPY<GivaroZpz<Std32> >
	{
	    public:

		typedef GivaroZpz<Std32>::Element Element;
		typedef GivaroZpz<Std32> Field;

		FieldAXPY (const Field &F) : _F (F) , Corr(uint64(-1) % (uint64)F.characteristic() +1){ _y = 0; }
		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F), _y (0) , Corr(faxpy.Corr) {}

		FieldAXPY<GivaroZpz<Std32> > &operator = (const FieldAXPY &faxpy) 
			{ _F = faxpy._F; _y = faxpy._y; Corr = faxpy.Corr; return *this; }

		inline void accumulate (const Element &a, const Element &x)
		{
			uint64 t = (uint64) a * (uint64) x;
			_y += t;

			if (_y < t)
				_y += Corr;
		}

		inline Element &get (Element &y) {
			_y %= (uint64) _F.characteristic();
			if ((int64) _y < 0) _y += _F.characteristic();
			y = (uint32) _y;
			return y;
		}

		inline FieldAXPY &assign (const Element y)
			{ _y = y; return *this; }

		inline void reset() {
			_y = 0;
		}

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

		FieldAXPY (const Field &F) : _F (F) , Corr(uint32(-1) % (uint32)F.characteristic() +1){ _y = 0; }
		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F), _y (0) , Corr(faxpy.Corr) {}

		FieldAXPY<GivaroZpz<Std16> > &operator = (const FieldAXPY &faxpy) 
			{ _F = faxpy._F; _y = faxpy._y; Corr = faxpy.Corr; return *this; }

		inline void accumulate (const Element &a, const Element &x)
		{
			uint32 t = (uint32) a * (uint32) x;
			_y += t;

			if (_y < t)
				_y += Corr;
		}

		inline Element &get (Element &y) {
			_y %= (uint32) _F.characteristic();
			if ((int32) _y < 0) _y += _F.characteristic();
			y = (uint16) _y;
			return y;
		}

		inline FieldAXPY &assign (const Element y)
			{ _y = y; return *this; }

		inline void reset() {
			_y = 0;
		}

	    private:

		Field _F;
		uint32 _y;
		uint32 Corr;
	};



	// Specialization of DotProductDomain for GivaroZpz<Std32> field
	
	template <>
	class DotProductDomain<GivaroZpz<Std32> > 
	{
	protected:
		GivaroZpz<Std32> _F;
	public:
		
		typedef GivaroZpz<Std32>::Element Element;
		
		DotProductDomain (const GivaroZpz<Std32> &F)
			: _F(F) ,
			  Corr(uint64(-1) % (uint64)F.characteristic() +1),
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
	class DotProductDomain<GivaroZpz<Std16> > 
	{
	protected:
		GivaroZpz<Std16> _F;
	public:

		typedef GivaroZpz<Std16>::Element Element;

		DotProductDomain (const GivaroZpz<Std16> &F)
			: _F(F),
			  Corr(uint32(-1) % (uint32)F.characteristic() +1),
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
		
} // namespace LinBox

#include "linbox/field/givaro-zpz.inl"

#endif // __FIELD_GIVARO_ZPZ
