/* linbox/field/modular.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * LargeModular is now replace by a class Modular parameterized on the element
 * type. So, the old LargeModular is equivalent to Modular<integer>. All other
 * interface details are exactly the same.
 *
 * Renamed from large-modular.h to modular.h
 * ------------------------------------
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file field/modular.h
 * @ingroup field
 * @brief A Modular field is a representations of <code>Z/mZ</code>.
 * This file groups many implementations/specialisations of modular fields.
 *   - Modular arithmetic is provided in the <code>ModularXXX<T></code> classes.
 *   - Specialisations for \ref FieldAXPY, \ref MVProductDomain, \ref DotProductDomain.
 *   - Random Iterators
 *   .
 *
 * @bug move Element& init(const Element&) to FFPACK. use using more..
 */

#ifndef __LINBOX_field_modular_H
#define __LINBOX_field_modular_H

#include <iostream>
#include <climits>
#include <cmath>

#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/util/field-axpy.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/linbox-config.h"
#include "linbox/field/field-traits.h"


// Namespace in which all LinBox code resides
namespace LinBox
{ /* Modular Base */

	template <class Element>
	class Modular;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<Modular<Element> >
	{
		typedef RingCategories::ModularTag categoryTag;
	};

	template <class Element>
	struct ClassifyRing<Modular<Element> const>
	{
		typedef RingCategories::ModularTag categoryTag;
	};


	/** * <!-- @name ModularBase-->
	 * \brief Base for prime fields where the elements are represented by various primitive types (and their operations).
	 * \ingroup field
	 * \defgroup Fields Fields <!--for old \ref Fields...-->
	 *
	 *
	 * Normally use it's children.  This class is of interest for the
	 * developer of a new field representation.
	 *
	 *
	 * This parameterized field can be used to construct any prime field.
	 * Typical use would be Modular<integer> for integers modulo a large
	 * prime, Modular<long> Modular<long long> for integers modulo a wordsize
	 * prime, etc. for integers modulo a half-wordsize prime.
	 */
	template <class _Element>
	class ModularBase {
	public:

		/*- Element type
		*/
		typedef _Element Element;

		/*- Random iterator generator type.
		 * It must meet the common object interface of random element generators
		 * as given in the the archetype RandIterArchetype.
		 */
		class RandIter;

		/*- @name Object Management
		*/
		//@{

		/*- Default constructor.
		*/
		ModularBase (void) {}

		/*- Constructor from an element type.
		 * Sets the modulus of the field throug the static member of the
		 * element type.
		 * @param modulus constant reference to integer prime modulus
		 */
		ModularBase (unsigned long modulus, int e = 1) :
			_modulus ((Element)modulus)
		{}

		/*- Constructor from an integer.
		 * Sets the modulus of the field throug the static member of the
		 * element type.
		 * @param modulus constant reference to integer prime modulus
		 */
		ModularBase (const integer &modulus) :
			_modulus ((Element) modulus)
		{}

		/*- Copy constructor.
		 * Constructs Modular object by copying the field.
		 * This is required to allow field objects to be passed by value
		 * into functions.
		 * @param  F Modular object.
		 */
		ModularBase (const ModularBase<Element> &F) :
			_modulus (F._modulus)
		{}

		/*- Conversion of field base element to a template class T.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to template class T.
		 * @param x template class T to contain output (reference returned).
		 * @param y constant field base element.
		 */
		integer &convert (integer &x, const Element &y) const
		{
			return x = y;
		}

		double &convert (double& x, const Element &y) const
		{
			return  x= (double) y;
		}

		float &convert (float& x, const Element &y) const
		{
			return  x= (float) y;
		}

		/*- Assignment of one field base element to another.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &assign (Element &x, const Element &y) const
		{
			return x = y;
		}

		/*- Cardinality.
		 * Return integer representing cardinality of the domain.
		 * Returns a non-negative integer for all domains with finite
		 * cardinality, and returns -1 to signify a domain of infinite
		 * cardinality.
		 * @return integer representing cardinality of the domain
		 */
		integer &cardinality (integer &c) const
		{
			return c = _modulus;
		}

		integer cardinality () const
		{
			return  _modulus;
		}

		/*- Characteristic.
		 * Return integer representing characteristic of the domain.
		 * Returns a positive integer to all domains with finite characteristic,
		 * and returns 0 to signify a domain of infinite characteristic.
		 * @return integer representing characteristic of the domain.
		 */
		integer &characteristic (integer &c) const
		{
			return c = _modulus;
		}

		unsigned long &characteristic (unsigned long &c) const
		{
			return c = _modulus;
		}


		integer characteristic () const
		{
			return  _modulus;
		}

		//@} Object Management

		/*- @name Arithmetic Operations
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field base elements will
		 * give undefined results.
		 */
		//@{

		/*- Equality of two elements.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return boolean true if equal, false if not.
		 * @param  x field base element
		 * @param  y field base element
		 */
		bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		/*- Zero equality.
		 * Test if field base element is equal to zero.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals zero, false if not.
		 * @param  x field base element.
		 */
		bool isZero (const Element &x) const
		{
			return x == 0;
		}

		/*- One equality.
		 * Test if field base element is equal to one.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals one, false if not.
		 * @param  x field base element.
		 */
		bool isOne (const Element &x) const
		{
			return x == 1;
		}


		//@} Arithmetic Operations

		/*- @name Input/Output Operations */
		//@{

		/*- Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream &write (std::ostream &os) const
		{
			return os << "Modular field, mod " << _modulus;
		}

		/*- Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream &read (std::istream &is) {
			return is >> _modulus;
		}


		/*- Print field base element.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return output stream to which field base element is written.
		 * @param  os  output stream to which field base element is written.
		 * @param  x   field base element.
		 */
		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << (int) x;
		}


		/*- Read field base element.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return input stream from which field base element is read.
		 * @param  is  input stream from which field base element is read.
		 * @param  x   field base element.
		 */
		std::istream &read (std::istream &is, Element &x) const
		{
			integer tmp;

			is >> tmp;

			x = abs (tmp) % integer (_modulus);
			if (tmp < 0) x = _modulus - x;

			return is;
		}

		//@}

	protected:

		/// Private (non-static) element for modulus
		Element _modulus;

	}; // class ModularBase
}

namespace LinBox
{ /* Modular */

	/* .. such comments as here should be on specialization...
	 * @param element Element type, e.g. long or integer
	 * @param Intermediate Type to use for intermediate computations. This
	 *                     should be a data type that can support integers
	 *                     twice the length of the maximal modulus used.
	 *
	 * The primality of the modulus will not be checked, so it is the
	 * programmer's responsibility to supply a prime modulus.  This class
	 * implements a field of unparameterized integers modulo a prime integer.
	 * Field has (non-static) member to contain modulus of field.
	 */

	/** @brief Prime fields of positive characteristic implemented directly in LinBox.
	 *
	 * This parameterized field can be used to construct prime fields.
	 * Typical use would be Modular<integer> for integers modulo a large
	 * prime, Modular<uint32_t>, Modular<int32_t>, or Modular<double> for
	 * integers modulo a wordsize prime.  Each of those has specialized
	 * performance features suitable to certain applications.
	 */
	template <class _Element>
	class Modular : public ModularBase<_Element> {
	public:
		typedef _Element Element;
		typedef Modular<_Element>     Self_t;
		typedef ModularBase<_Element> Father_t;
		typedef typename ModularBase<_Element>::RandIter RandIter;
		const Element zero,one, mOne;

		/*- @name Object Management
		 * @brief see \ref FieldArchetype  for member specs.
		 */
		//@{

		//private:
		/*- Default constructor.
		*/
		Modular () :
			zero(0),one(1),mOne(0)
		{}

		/*- Constructor from an element type
		 * Sets the modulus of the field throug the static member of the
		 * element type.
		 * @param modulus constant reference to integer prime modulus
		 */
		Modular (unsigned long modulus, unsigned long = 1) :
			ModularBase<_Element> (modulus),zero(0),one(1),mOne(modulus-1)
		{}


		/*- Constructor from an integer
		 * Sets the modulus of the field throug the static member of the
		 * element type.
		 * @param modulus constant reference to integer prime modulus
		 */
		Modular (const integer &modulus) :
			ModularBase<_Element> (modulus),zero(0),one(1),mOne(modulus-1)
		{
			linbox_check( !IsNegative(modulus) );
		}

		/* Assignment operator
		 * Required by the archetype
		 *
		 * @param F constant reference to Modular object
		 * @return reference to Modular object for self
		 */
		const Modular &operator=(const Modular &F)
		{
			ModularBase<Element>::_modulus = F._modulus;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);

			return *this;
		}
	public:

		//! @warning danger ! (Element = Integer or NTL stuff !):
		static inline Element getMaxModulus()
		{
			return Element((1ULL<<(sizeof(Element)*8-1))-1);
		}


		/*- Initialization of field base element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * This is not a specialization of the template function because
		 * such a specialization is not allowed inside the class declaration.
		 * @return reference to field base element.
		 * @param x field base element to contain output (reference returned).
		 * @param y integer.
		 */
		Element &init (Element &x, const integer &y ) const
		{
			x = y % ModularBase<Element>::_modulus;
			if ( IsNegative(x) )
				x += ModularBase<Element>::_modulus;
			return x;
		}

		Element &init (Element &x, const size_t &y ) const
		{
			x = (Element) y % ModularBase<Element>::_modulus;
			if ( IsNegative(x) )
				x += ModularBase<Element>::_modulus;
			return x;
		}

		Element &init (Element &x, const int y ) const
		{
			x = y % ModularBase<Element>::_modulus;
			if ( IsNegative(x) )
				x += ModularBase<Element>::_modulus;
			return x;
		}

		Element &init (Element &x, const long int y) const
		{
			x = y % ModularBase<Element>::_modulus;
			if ( IsNegative(x) )
			       	x += ModularBase<Element>::_modulus;
			return x;
		}

		/*- Initialization of field base element from a double.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * This is not a specialization of the template function because
		 * such a specialization is not allowed inside the class declaration.
		 * @return reference to field base element.
		 * @param x field base element to contain output (reference returned).
		 * @param y integer.
		 */
		Element &init (Element &x, const double &y) const
		{
			double z = fmod(y, (double)ModularBase<Element>::_modulus);
			if ( z < 0 )
				z += (double) ModularBase<Element>::_modulus;
			return x = (Element) (z+.5);
		}

		Element &init (Element &x, const float &y) const
		{
			float z = fmod(y, (float)ModularBase<Element>::_modulus);
			if ( z < 0 )
			       	z += (float) ModularBase<Element>::_modulus;
			return x = (Element) (z+.5);
		}

		Element &init(Element &x) const
		{
			return x = 0 ;
		}

		//@}
		/*- @name Arithmetic Operations
		 * @brief see \ref FieldArchetype  for member specs.
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field base elements will
		 * give undefined results.
		 */
		//@{

		/*- Addition.
		 * x = y + z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			if (x >= ModularBase<Element>::_modulus) x -= ModularBase<Element>::_modulus;
			return x;
		}

		/* Subtraction.
		 * x = y - z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if ( IsNegative(x) )
			       	x += ModularBase<Element>::_modulus;
			return x;
		}

		/* Multiplication.
		 * x = y * z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			return x = (y * z) % ModularBase<Element>::_modulus;
		}

		/* Division.
		 * x = y / z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		/* Additive Inverse (Negation).
		 * x = - y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &neg (Element &x, const Element &y) const
		{
			if (y == 0)
				return x = y;
			else
				return x = ModularBase<Element>::_modulus - y;
		}

		/* Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &inv (Element &x, const Element &y) const
		{
			// The extended Euclidean algoritm
			Element x_int, y_int, q, tx, ty, temp;
			x_int = ModularBase<Element>::_modulus;
			y_int = y;
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int;  y_int  = x_int  - q*y_int;  x_int  = temp;
				temp = ty; ty = tx - q*ty; tx = temp;
			}

			// now x_int = gcd (modulus,residue)
			x = tx;
			if ( IsNegative(x) )
			       	x += ModularBase<Element>::_modulus;

			return x;
		}

		/* Natural AXPY.
		 * r  = a * x + y
		 * This function assumes all field elements have already been
		 * constructed and initialized.
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 * @param  y field element.
		 */
		Element &axpy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			r = (a * x + y) % ModularBase<Element>::_modulus;
			if ( IsNegative(r) )
			       	r += ModularBase<Element>::_modulus;
			return r;
		}

		//@} Arithmetic Operations

		/*- @name Inplace Arithmetic Operations
		 * @brief see \ref FieldArchetype  for member specs.
		 * x <- x op y; x <- op x
		 */
		//@{

		/*- Inplace Addition.
		 * x += y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if (x >= ModularBase<Element>::_modulus) x -= ModularBase<Element>::_modulus;
			return x;
		}

		/* Inplace Subtraction.
		 * x -= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if ( IsNegative(x) )
				x += ModularBase<Element>::_modulus;
			return x;
		}

		/* Inplace Multiplication.
		 * x *= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &mulin (Element &x, const Element &y) const
		{
			x *= y;
			x %= ModularBase<Element>::_modulus;
			return x;
		}

		/* Inplace Division.
		 * x /= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &divin (Element &x, const Element &y) const
		{
			Element temp;
			inv (temp, y);
			return mulin (x, temp);
		}

		/* Inplace Additive Inverse (Inplace Negation).
		 * x = - x
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element &negin (Element &x) const
		{
			if (x == 0)
				return x;
			else
				return x = ModularBase<Element>::_modulus - x;
		}

		/* Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes the field base elementhas already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		/* Inplace AXPY.
		 * r  += a * x
		 * This function assumes all field elements have already been
		 * constructed and initialized.
		 * Purely virtual
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 */
		Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r = (r + a * x) % ModularBase<Element>::_modulus;
			if ( IsNegative(r) )
			       	r += ModularBase<Element>::_modulus;
			return r;
		}

		//@} Inplace Arithmetic Operations

	private:

		friend class FieldAXPY<Modular<Element> >;

	}; // class Modular


	/*! Specialization of FieldAXPY for parameterized modular field */

	template <class _Element>
	class FieldAXPY<Modular<_Element> > {
	public:

		typedef _Element Element;
		typedef Modular<_Element> Field;

		FieldAXPY (const Field &F) :
			_field (F)
		{ _y = 0; }
		FieldAXPY (const FieldAXPY<Modular<Element> > &faxpy) :
			_field (faxpy._field), _y (faxpy._y)
		{}

		FieldAXPY<Modular <Element> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			return *this;
		}

		inline Element& mulacc (const Element &a, const Element &x)
		{
			return accumulate(a * x);
		}

		inline Element& accumulate (const Element &t)
		{
			return _y+=t;
		}

		inline Element &get (Element &y) { _y %= _field._modulus; y = _y; return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			_field.init(_y, 0);
		}

	private:

		Field _field;
		Element _y;
	};


	template <>
	inline std::ostream& ModularBase<Integer>::write (std::ostream &os) const
	{
		return os << "GMP integers mod " << _modulus;
	}

	template <>
	inline integer& Modular<integer>::init (integer& x, const double& y) const
	{
		integer tmp = (integer)y % _modulus;
		if (tmp<0) tmp += _modulus;
		return x = tmp;
	}

	template<>
	integer Modular<integer>::getMaxModulus()
	{
		return -1 ;
	}



} // namespace LinBox

#include "linbox/field/Modular/modular-unsigned.h"
#include "linbox/randiter/modular.h"
#include "linbox/field/Modular/modular-int32.h"
#ifdef __LINBOX_HAVE_INT64
#include "linbox/field/Modular/modular-int64.h"
#endif
#include "linbox/field/Modular/modular-short.h"
#include "linbox/field/Modular/modular-byte.h"
#include "linbox/field/Modular/modular-double.h"
#include "linbox/field/Modular/modular-float.h"

#endif // __LINBOX_field_modular_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

