/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* linbox/field/unparametric.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
 *               2005 Clement Pernet
 *
 * Written by W. J. Turner <wjturner@acm.org>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_field_unparametric_H
#define __LINBOX_field_unparametric_H
#include <typeinfo>

#include <string>
#include <algorithm>

#include "linbox/integer.h"
#include <linbox/field/field-interface.h>
#include "linbox/randiter/unparametric.h"
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>
#include <fflas-ffpack/field/unparametric.h>
//#if __LINBOX_HAVE_NTL
//#include <linbox/field/ntl-RR.h>
//#endif // __LINBOX_HAVE_NTL


namespace LinBox
{

	template <typename Target, typename Source>
	Target& Caster (Target& t, const Source& s)
	{
		return t = static_cast<Target>(s);
	}

#if 0
#if __LINBOX_HAVE_NTL
	typedef NTL::RR Targ;
	template <>
	Targ& Caster<Targ, int> (Targ& t, const int& s)
	{
		return t = s;
	}
#endif // __LINBOX_HAVE_NTL
#endif

	template <class Ring>
	struct ClassifyRing;

	template <class K>
	class UnparametricField;

	template <class K>
	struct ClassifyRing<UnparametricField<K> > {
		typedef RingCategories::GenericTag categoryTag;
	};


	/** \brief Unparameterized field adapter.
	 * \ingroup field
	 * \defgroup UnparametricField UnparametricField
	 *
	 * A field having an interface similar to that of floats is adapted to
	 * LinBox.
	 *
	 *  Used to generate efficient field classes for unparameterized fields
	 *  (or hidden parameter fields).
	 *
	 *  Some fields are implemented by definition of the C++ arithmetic
	 *  operators, such as <code>z = x*y</code>, for \c x, \c y, \c z
	 *  instances of a type \c K.  The LinBox field LinBox::Unparametric<K>
	 *  is the adaptation to LinBox.
	 *
	 *  For a typical unparametric field, some of the methods must be
	 *  defined in a specialization.
	 */
	template <class K>
	class UnparametricField : public FieldInterface,
		public  FFPACK::UnparametricField<K> {
	protected:
		integer _p; integer _card;
	public:

		/** @name Common Object Interface for a LinBox Field.
		 * These methods and member types are required of all LinBox fields.
		 * See \ref FieldArchetype  for detailed specifications.
		 */
		//@{

		/** The field's element type.
		 * Type K must provide a default constructor,
		 * a copy constructor, a destructor, and an assignment operator.
		 */

		typedef K Element;

		/// Type of random field element generators.
		typedef UnparametricRandIter<K> RandIter;

		/** @name Field Object Basics.
		*/
		//@{

		/** Builds this field to have characteristic q and cardinality q<sup>e</sup>.
		 *  This constructor must be defined in a specialization.
		 */
		UnparametricField(integer q = 0, size_t e = 1) :
			FFPACK::UnparametricField<K>((unsigned long)q,(unsigned long)e),
			_p(q), _card(q == 0 ? integer(-1) : pow(q, e) )
			{}  // assuming q is a prime or zero.

		/// construct this field as copy of F.
		UnparametricField (const UnparametricField &F) :
			FFPACK::UnparametricField<K>(F),_p(F._p), _card(F._card)
		{}


		// field/ntl-ZZ_p.h me les demande... //

		Element&inv (Element&x, const Element&y)const{return FFPACK::UnparametricField<K>::inv(x,y);}
		Element&invin (Element&x)const{return FFPACK::UnparametricField<K>::invin(x);}
		std::ostream&write (std::ostream&os)const{return FFPACK::UnparametricField<K>::write(os);}
		std::ostream&write (std::ostream&os, const Element&p)const{return FFPACK::UnparametricField<K>::write(os,p);}
		bool isZero(const Element&x)const{return FFPACK::UnparametricField<K>::isZero(x);}
		bool isOne(const Element&x)const{return FFPACK::UnparametricField<K>::isOne(x);}
		long unsigned int characteristic(long unsigned int&p)const{return FFPACK::UnparametricField<K>::characteristic(p);}
		long unsigned int characteristic()const{return FFPACK::UnparametricField<K>::characteristic();};
		template<typename Src>Element&init(Element&x, const Src&s)const{return Caster (x, s);}
		std::istream&read(std::istream&is, Element&x)const{return FFPACK::UnparametricField<K>::read(is,x);}
		std::istream&read(std::istream&is)const{return FFPACK::UnparametricField<K>::read(is);}
		template<typename T>T&convert(T&x,const Element&y)const{return Caster(x,y);}

		// fin des trucs zarbs //

		/// c := cardinality of this field (-1 if infinite).
		integer &cardinality (integer &c) const
		{
			return c = _card;
		}

		/// c := characteristic of this field (zero or prime).
		integer &characteristic (integer &c) const
		{
			return c = _p;
		}

		//@} Data Object Management




		//@} Common Object Interface

		/** @name Implementation-Specific Methods.
		 * These methods are not required of all LinBox fields
		 * and are included only for the implementation of this field
		 * template.
		 */
		//@{

		//- Default constructor
		//UnparametricField (void) {}

		/** Constructor from field object.
		 * @param  A unparameterized field object
		 */
		UnparametricField (const K &A) {}

		/** Constant access operator.
		 * @return constant reference to field object
		 */
		const K &operator () (void) const
		{
			return Element ();
		}

		/** Access operator.
		 * @return reference to field object
		 */
		K &operator () (void)
		{
			return Element ();
		}

		//@} Implementation-Specific Methods

	}; // template <class K> class UnparametricField

	template<class Field>
	class FieldAXPY;

	/*! @ingroup integers
	 * @brief NO DOc
	 */
	template<>
	class FieldAXPY<UnparametricField<integer> >  {
	public:
		typedef UnparametricField<integer> Field;
		typedef integer Element;

		/** Constructor.
		 * A faxpy object if constructed from a Field and a field element.
		 * Copies of this objects are stored in the faxpy object.
		 * @param F field F in which arithmetic is done
		 */
		FieldAXPY (const Field &F) :
			_F (F)
		{ _y = 0; }

		/** Copy constructor.
		 * @param faxpy
		 */
		FieldAXPY (const FieldAXPY<Field> &faxpy) :
			_F (faxpy._F), _y (faxpy._y)
		{}

		/** Assignment operator
		 * @param faxpy
		 */
		FieldAXPY<Field> &operator = (const FieldAXPY &faxpy)
		{ _y = faxpy._y; return *this; }

		/** Add a*x to y
		 * y += a*x.
		 * @param a constant reference to element a
		 * @param x constant reference to element x
		 * allow optimal multiplication, such as integer * int
		 */
		template<class Element1>
		inline Element&  mulacc (const Element &a, const Element1 &x)
		{
			return _y += (a * x);
		}

		template<class Element1>
		inline Element& accumulate (const Element1 &t)
		{
			return _y += t;
		}

		/** Retrieve y
		 *
		 * Performs the delayed modding out if necessary
		 */
		inline Element &get (Element &y) { y = _y; return y; }

		/** Assign method.
		 * Stores new field element for arithmetic.
		 * @return reference to self
		 * @param y constant reference to element a
		 */
		inline FieldAXPY &assign (const Element& y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			_y = 0;
		}

	private:

		/// Field in which arithmetic is done
		Field _F;

		/// Field element for arithmetic
		Element _y;

	};

} // namespace LinBox

#include "linbox/randiter/unparametric.h"

#endif // __LINBOX_field_unparametric_H

