/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
 * See COPYING for license information.
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
{
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
		ModularBase (unsigned long modulus) :
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
		typedef typename ModularBase<_Element>::RandIter RandIter;
		const Element zero,one;
		Element mone;

		/*- @name Object Management
		 * @brief see \ref FieldArchetype  for member specs.
		 */
		//@{

		//private:
		/*- Default constructor.
		*/
		Modular () :
			zero(0),one(1)
		{}

		/*- Constructor from an element type
		 * Sets the modulus of the field throug the static member of the
		 * element type.
		 * @param modulus constant reference to integer prime modulus
		 */
		Modular (unsigned long modulus) :
			ModularBase<_Element> (modulus),zero(0),one(1)
		{}

		/*- Constructor from an integer
		 * Sets the modulus of the field throug the static member of the
		 * element type.
		 * @param modulus constant reference to integer prime modulus
		 */
		Modular (const integer &modulus) :
			ModularBase<_Element> (modulus),zero(0),one(1)
		{}

		/* Assignment operator
		 * Required by the archetype
		 *
		 * @param F constant reference to Modular object
		 * @return reference to Modular object for self
		 */
		const Modular &operator=(const Modular &F)
		{
			ModularBase<Element>::_modulus = F._modulus;
			return *this;
		}
	public:

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
		Element &init (Element &x, const integer &y = 0) const
		{
			x = y % ModularBase<Element>::_modulus;
			if (x < 0) x += ModularBase<Element>::_modulus;
			return x;
		}

		Element &init (Element &x, const size_t &y = 0) const
		{
			x = (Element) y % ModularBase<Element>::_modulus;
			if (x < 0) x += ModularBase<Element>::_modulus;
			return x;
		}

		Element &init (Element &x, const int y =0) const
		{
			x = y % ModularBase<Element>::_modulus;
			if (x < 0) x += ModularBase<Element>::_modulus;
			return x;
		}

		Element &init (Element &x, const long int y =0) const
		{
			x = y % ModularBase<Element>::_modulus;
			if (x < 0) x += ModularBase<Element>::_modulus;
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
		Element &init (Element &x, const double &y=0) const
		{
			double z = fmod(y, (double)ModularBase<Element>::_modulus);
			if (z < 0) z += (double) ModularBase<Element>::_modulus;
			return x = (Element) (z+.5);
		}

		Element &init (Element &x, const float &y=0) const
		{
			float z = fmod(y, (float)ModularBase<Element>::_modulus);
			if (z < 0) z += (float) ModularBase<Element>::_modulus;
			return x = (Element) (z+.5);
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
			if (x < 0) x += ModularBase<Element>::_modulus;
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
			if (x < 0) x += ModularBase<Element>::_modulus;

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
			if (r < 0) r += ModularBase<Element>::_modulus;
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
			if (x < 0) x += ModularBase<Element>::_modulus;
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
			if (r < 0) r += ModularBase<Element>::_modulus;
			return r;
		}

		//@} Inplace Arithmetic Operations

	private:

		friend class FieldAXPY<Modular<Element> >;

	}; // class Modular

	/** @brief Allows compact storage when the modulus is less than 2^8.
	 *
	 *  Requires <code>1 < modulus < 2^8</code>, normally prime.  See \ref
	 *  FieldArchetype for member specifications.
	 */
	template <>
	class Modular<uint8_t> : public FieldInterface, public ModularBase<uint8_t> {
	public:
		typedef uint8_t Element;
		const Element zero,one;
		Element mone;

		Modular () :
			zero(0),one(1),_k (0)
		{}
		Modular (uint32_t modulus) :
			ModularBase<Element> (modulus),
			zero(0),one(1),mone((Element)(modulus-1)),
			_k (((uint64_t) -1LL) / ((modulus - 1) * (modulus - 1))),
			_pinv (1.0 / (double) ((Element) modulus))
		{
			linbox_check(modulus < UINT8_MAX);
		}
		Modular (const integer &modulus) :
			ModularBase<Element> ((long) modulus),
			zero(0),one(1),mone(modulus-1),
			_k (((uint64_t) -1LL) / (((Element)modulus - 1) * ((Element)modulus - 1))),
			_pinv (1.0 / (double) ((Element) modulus))
		{
			linbox_check(modulus < UINT8_MAX);
		}

		const Modular &operator=(const Modular &F)
		{
			ModularBase<Element>::_modulus = F._modulus;
			_k = F._k;
			_pinv = F._pinv;
			mone = F.mone;
			return *this;
		}

		Element &init (Element &x, const integer &y = 0) const
		{
			x = (Element) (abs (y) % integer (ModularBase<Element>::_modulus));
			if (y < 0)
				x = Element(ModularBase<Element>::_modulus - x);
			return x;
		}

		Element &init (Element &x, const double &y) const
		{
			double z = fmod(y, (double)_modulus);
			if (z < 0) z += (double) _modulus;
			return x = (Element) (z);
		}

		/*! add elements
		 * @todo is it faster to use uint32 and multiple casts ?
		 */
		Element &add (Element &x, const Element &y, const Element &z) const
		{
			uint32_t t = (uint32_t) y + (uint32_t) z;
			if (t >= (uint32_t) ModularBase<Element>::_modulus)
				t -= ModularBase<Element>::_modulus;
			return x = (Element)t;
		}

		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			int32_t t = (int32_t) y - (int32_t) z;
			if (t < 0)
				t += ModularBase<Element>::_modulus;
			return x =  (Element)t;
		}

		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			return x = Element( ((uint32_t) y * (uint32_t) z) % (uint32_t) ModularBase<Element>::_modulus );
		}

		Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		Element &neg (Element &x, const Element &y) const
		{
			if (y == 0)
				return x = y;
			else
				return x = (Element) (ModularBase<Element>::_modulus - y);
		}

		Element &inv (Element &x, const Element &y) const
		{
			// The extended Euclidean algoritm
			int32_t x_int, y_int, q, tx, ty, temp;
			x_int = ModularBase<Element>::_modulus;
			y_int = y;
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}

			if (tx < 0) tx += ModularBase<Element>::_modulus;

			// now x_int = gcd (modulus,residue)
			return x = (Element) tx;
		}

		Element &axpy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			r = Element(((uint32_t) a * (uint32_t) x + (uint32_t) y) % (uint32_t) ModularBase<Element>::_modulus) ;
			return r;
		}

		Element &addin (Element &x, const Element &y) const
		{
			uint32_t t = uint32_t((long) x + (long) y);
			if (t >= (uint32_t) ModularBase<Element>::_modulus) t -= ModularBase<Element>::_modulus;
			return x = (Element) t;
		}

		/*! subin.
		 * @todo why \c long here ?
		 */
		Element &subin (Element &x, const Element &y) const
		{
			long t = x - y;
			if (t < 0) t += ModularBase<Element>::_modulus;
			return x = (Element) t;
		}

		Element &mulin (Element &x, const Element &y) const
		{
			x = (Element)( ((uint32_t) x * (uint32_t) y) % (uint32_t) ModularBase<Element>::_modulus );
			return x;
		}

		Element &divin (Element &x, const Element &y) const
		{
			Element temp;
			inv (temp, y);
			return mulin (x, temp);
		}

		Element &negin (Element &x) const
		{
			if (x == 0)
				return x;
			else
				return x = Element(ModularBase<Element>::_modulus - x);
		}

		Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r = (Element)( ((uint32_t) r + (uint32_t) a * (uint32_t) x) % (uint32_t) ModularBase<Element>::_modulus);
			return r;
		}

	private:

		friend class FieldAXPY<Modular<Element> >;
		friend class DotProductDomain<Modular<Element> >;
		friend class MVProductDomain<Modular<Element> >;

		// Number of times one can perform an axpy into a long long
		// before modding out is mandatory.
		uint64_t _k;

		// Inverse of modulus in floating point
		double _pinv;

	}; // class Modular<uint8_t>

	/** @brief Specialization of class Modular for uint16_t element type */
	template <>
	class Modular<uint16_t> : public FieldInterface, public ModularBase<uint16_t> {
	public:

		typedef uint16_t Element;

		const Element zero,one;
		Element mone;

		Modular () :
			zero(0),one(1),_k (0)
		{}
		Modular (uint32_t modulus) :
			ModularBase<Element> (modulus),
			zero(0),one(1),mone((Element)(modulus-1)),
			_k (((uint64_t) -1LL) / ((ModularBase<Element>::_modulus - 1) * (ModularBase<Element>::_modulus - 1))),
			_pinv (1.0 / (double) ((Element) ModularBase<Element>::_modulus))
		{
			linbox_check(modulus<UINT16_MAX);
		}
		Modular (const integer &modulus) :
			ModularBase<Element> ((long) modulus),
			zero(0),one(1),mone(Element(modulus-1)),
			_k (((uint64_t) -1LL) / ((ModularBase<Element>::_modulus - 1) * (ModularBase<Element>::_modulus - 1))),
			_pinv (1.0 / (double) ((Element) ModularBase<Element>::_modulus))
		{
			linbox_check(modulus<UINT16_MAX);
		}

		const Modular &operator=(const Modular &F)
		{
			ModularBase<Element>::_modulus = F._modulus;
			_k = F._k;
			_pinv = F._pinv;
			return *this;
		}

		Element &init (Element &x, const integer &y = 0) const
		{
			x = abs (y) % integer (ModularBase<Element>::_modulus);
			if (y < 0)
				x = Element(ModularBase<Element>::_modulus - x);
			return x;
		}

		Element &init (Element &x, const double &y) const
		{
			double z = fmod(y, (double)_modulus);
			if (z < 0)
				z += (double) _modulus;
			return x = (Element) (z);
		}

		Element &add (Element &x, const Element &y, const Element &z) const
		{
			uint32_t t = (uint32_t) y + (uint32_t) z;
			if (t >= (uint32_t) ModularBase<Element>::_modulus)
				t -= ModularBase<Element>::_modulus;
			return x = (Element) t;
		}

		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			int32_t t = (int32_t) y - (int32_t) z;
			if (t < 0) t += ModularBase<Element>::_modulus;
			return x =  (Element) t;
		}

		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			return x = (Element) ( ((uint32_t) y * (uint32_t) z) % (uint32_t) ModularBase<Element>::_modulus);
		}

		Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		Element &neg (Element &x, const Element &y) const
		{
			if (y == 0)
				return x = y;
			else
				return x = (Element)(  ModularBase<Element>::_modulus - y);
		}

		Element &inv (Element &x, const Element &y) const
		{
			// The extended Euclidean algoritm
			int32_t x_int, y_int, q, tx, ty, temp;
			x_int = ModularBase<Element>::_modulus;
			y_int = y;
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}

			if (tx < 0) tx += ModularBase<Element>::_modulus;

			// now x_int = gcd (modulus,residue)
			return x = (Element)  tx;
		}

		Element &axpy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			r =  (Element)( ((uint32_t) a * (uint32_t) x + (uint32_t) y) % (uint32_t) ModularBase<Element>::_modulus );
			return r;
		}

		Element &addin (Element &x, const Element &y) const
		{
			uint32_t t = uint32_t( (long) x + (long) y );
			if (t >= (uint32_t) ModularBase<Element>::_modulus) t -= ModularBase<Element>::_modulus;
			return x = (Element)  t;
		}

		Element &subin (Element &x, const Element &y) const
		{
			long t = x - y;
			if (t < 0) t += ModularBase<Element>::_modulus;
			return x =  (Element) t;
		}

		Element &mulin (Element &x, const Element &y) const
		{
			x =  (Element)( ((uint32_t) x * (uint32_t) y) % (uint32_t) ModularBase<Element>::_modulus);
			return x;
		}

		Element &divin (Element &x, const Element &y) const
		{
			Element temp;
			inv (temp, y);
			return mulin (x, temp);
		}

		Element &negin (Element &x) const
		{
			if (x == 0)
				return x;
			else
				return x = (Element) ( ModularBase<Element>::_modulus - x);
		}

		Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r = (Element) ( ((uint32_t) r + (uint32_t) a * (uint32_t) x) % (uint32_t) ModularBase<Element>::_modulus);
			return r;
		}

	private:

		friend class FieldAXPY<Modular<Element> >;
		friend class DotProductDomain<Modular<Element> >;
		friend class MVProductDomain<Modular<Element> >;

		// Number of times one can perform an axpy into a long long
		// before modding out is mandatory.
		uint64_t _k;

		// Inverse of modulus in floating point
		double _pinv;

	}; // class Modular<uint16_t>

	/** @brief Specialization of class Modular for uint32_t element type */
	template <>
	class Modular<uint32_t> : public FieldInterface, public ModularBase<uint32_t> {
	public:

		typedef uint32_t Element;

		const Element zero,one;
		Element mone;

		Modular () :
			zero(0),one(1)
		{}
		Modular (uint32_t modulus)  :
			ModularBase<uint32_t> (modulus),zero(0),one(1),mone(modulus-1)
		{ init_two_64 (); }
		Modular (const integer &modulus) :
			ModularBase<uint32_t> (modulus),zero(0),one(1),mone(modulus-1)
		{ init_two_64 (); }

		const Modular &operator=(const Modular &F)
		{
			ModularBase<Element>::_modulus = F._modulus;
			_two_64 = F._two_64;
			mone = F.mone;
			return *this;
		}

		Element &init (Element &x, const integer &y ) const
		{
			x = abs (y) % integer (ModularBase<Element>::_modulus);
			if (y < 0) x = ModularBase<Element>::_modulus - x;
			return x;
		}

		Element &init (Element &x, const long int &y ) const
		{
			x = abs (y) % integer (ModularBase<Element>::_modulus);
			if (y < 0) x = ModularBase<Element>::_modulus - x;
			return x;
		}

		Element &init (Element &x, const int &y ) const
		{
			x = abs (y) % integer (ModularBase<Element>::_modulus);
			if (y < 0) x = ModularBase<Element>::_modulus - x;
			return x;
		}

		Element &init (Element &x, const long unsigned int &y ) const
		{
			x = Element(y %  (ModularBase<Element>::_modulus));
			return x;
		}

		Element &init (Element &x, const double &y) const
		{
			double z = fmod(y, (double)_modulus);
			if (z < 0) z += (double) _modulus;
			return x = (Element) (z);
		}

		template< class XXX>
		Element& init(Element & x, const XXX & y) const
		{
			return init(x,double(y));
		}

		Element &init (Element &x) const
		{
			return x = zero ;
		}


		Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = y + z;
			if ((uint32_t) x >= (uint32_t) ModularBase<Element>::_modulus) x -= ModularBase<Element>::_modulus;
			return x;
		}

		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = y - z;
			if ((int32_t) x < 0) x += ModularBase<Element>::_modulus;
			return x;
		}

		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			return x = Element( ((uint64_t) y * (uint64_t) z) % (uint64_t) ModularBase<Element>::_modulus);
		}

		Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		Element &neg (Element &x, const Element &y) const
		{
			if (y == 0)
				return x = y;
			else
				return x = ModularBase<Element>::_modulus - y;
		}

		Element &inv (Element &x, const Element &y) const
		{
			// The extended Euclidean algoritm
			int64_t x_int, y_int, q, tx, ty, temp;
			x_int = ModularBase<Element>::_modulus;
			y_int = y;
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int;  y_int  = x_int  - q * y_int;
				x_int  = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}

			if (tx < 0) tx += ModularBase<Element>::_modulus;

			// now x_int = gcd (modulus,residue)
			return x = Element(tx);
		}

		Element &axpy (Element &r,
			       const Element &a,
			       const Element &x,
			       const Element &y) const
		{
			r = Element( ((uint64_t) a * (uint64_t) x + (uint64_t) y) % (uint64_t) ModularBase<Element>::_modulus );
			if ((int32_t) r < 0) r += ModularBase<Element>::_modulus;
			return r;
		}

		Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if ((uint32_t) x >= (uint32_t) ModularBase<Element>::_modulus) x -= ModularBase<Element>::_modulus;
			return x;
		}

		Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if ((int32_t) x < 0) x += ModularBase<Element>::_modulus;
			return x;
		}

		Element &mulin (Element &x, const Element &y) const
		{
			x = Element( ((uint64_t) x * (uint64_t) y) % (uint64_t) ModularBase<Element>::_modulus );
			return x;
		}

		Element &divin (Element &x, const Element &y) const
		{
			Element temp;
			inv (temp, y);
			return mulin (x, temp);
		}

		Element &negin (Element &x) const
		{
			if (x == 0)
				return x;
			else
				return x = ModularBase<Element>::_modulus - x;
		}

		Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r = Element( ((uint64_t) r + (uint64_t) a * (uint64_t) x) % (uint64_t) ModularBase<Element>::_modulus );
			if ((int32_t) r < 0) r += ModularBase<Element>::_modulus;
			return r;
		}

	private:

		void init_two_64 ()
		{
			uint64_t two_64 = 2;

			for (int i = 0; i < 6; ++i)
				two_64 = (two_64 * two_64) % ModularBase<Element>::_modulus;

			_two_64 = (Element) two_64;
		}

		friend class FieldAXPY<Modular<uint32_t> >;
		friend class DotProductDomain<Modular<uint32_t> >;
		friend class MVProductDomain<Modular<uint32_t> >;

		Element _two_64;

	}; // class Modular<uint32_t>

	/*! Specialization of FieldAXPY for parameterized modular field */

	template <class _Element>
	class FieldAXPY<Modular<_Element> > {
	public:

		typedef _Element Element;
		typedef Modular<_Element> Field;

		FieldAXPY (const Field &F) :
			_F (F)
		{ _y = 0; }
		FieldAXPY (const FieldAXPY<Modular<Element> > &faxpy) :
			_F (faxpy._F), _y (faxpy._y)
		{}

		FieldAXPY<Modular <Element> > &operator = (const FieldAXPY &faxpy)
		{
			_F = faxpy._F;
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

		inline Element &get (Element &y) { _y %= _F._modulus; y = _y; return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			_F.init(_y, 0);
		}

	private:

		Field _F;
		Element _y;
	};

	/*! Specialization of FieldAXPY for uint8_t modular field */

	template <>
	class FieldAXPY<Modular<uint8_t> > {
	public:

		typedef uint8_t Element;
		typedef Modular<uint8_t> Field;

		FieldAXPY (const Field &F) :
			_F (F),
			i ( (int)F._k)
		{
			_y = 0;
		}

		FieldAXPY (const FieldAXPY &faxpy) :
			_F (faxpy._F),
			_y (0),
			i ((int)faxpy._F._k)
		{}

		FieldAXPY<Modular<uint8_t> > &operator = (const FieldAXPY &faxpy)
		{
			_F = faxpy._F;
			_y = faxpy._y;
			return *this;
		}

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint32_t t = (uint32_t) a * (uint32_t) x;

			if (!i--) {
				i = int(_F._k);
				return _y = _y % (uint32_t) _F._modulus + t;
			}
			else
				return _y += t;
		}

		inline uint64_t& accumulate (const Element &t)
		{

			if (!i--) {
				i = int( _F._k );
				return _y = _y % (uint32_t) _F._modulus + t;
			}
			else
				return _y += t;
		}

		inline Element &get (Element &y)
		{
			_y %= (uint32_t) _F._modulus;
			if ((int32_t) _y < 0) _y += _F._modulus;
			y = (uint8_t) _y;
			i = int(_F._k);
			return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			i = int(_F._k);
			return *this;
		}

		inline void reset()
		{
			_y = 0;
		}

	private:

		Field _F;
		uint64_t _y;
		int i;
	};

	/*! Specialization of FieldAXPY for uint16_t modular field */
	template <>
	class FieldAXPY<Modular<uint16_t> > {
	public:

		typedef uint16_t Element;
		typedef Modular<uint16_t> Field;

		FieldAXPY (const Field &F) :
			_F (F),
			i ((int)F._k)
		{ _y = 0; }
		FieldAXPY (const FieldAXPY &faxpy) :
			_F (faxpy._F), _y (0), i ((int) faxpy._F._k)
		{}

		FieldAXPY<Modular<uint16_t> > &operator = (const FieldAXPY &faxpy)
		{
			_F = faxpy._F;
			_y = faxpy._y;
			return *this;
		}

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (long long) a * (long long) x;

			if (!i--) {
				i = (int) _F._k;
				return _y = _y % (uint64_t) _F._modulus + t;
			}
			else
				return _y += t;
		}

		inline uint64_t& accumulate (const Element &t)
		{
			if (!i--) {
				i = (int) _F._k;
				return _y = _y % (uint64_t) _F._modulus + t;
			}
			else
				return _y += t;
		}

		inline Element &get (Element &y)
		{
			_y %= (uint64_t) _F._modulus;
			if ((int64_t) _y < 0) _y += _F._modulus;
			y = (uint16_t) _y;
			i = (int) _F._k;
			return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			i = (int) _F._k;
			return *this;
		}

		inline void reset()
		{
			_y = 0;
		}

	private:

		Field _F;
		uint64_t _y;
		int i;
	};

	/*! Specialization of FieldAXPY for unsigned short modular field */

	template <>
	class FieldAXPY<Modular<uint32_t> > {
	public:

		typedef uint32_t Element;
		typedef Modular<uint32_t> Field;

		FieldAXPY (const Field &F) :
			_F (F), _y(0)
		{ }

		FieldAXPY (const FieldAXPY &faxpy) :
			_F (faxpy._F), _y (0)
		{}

		FieldAXPY<Modular<uint32_t> > &operator = (const FieldAXPY &faxpy)
		{
			_F = faxpy._F;
			_y = faxpy._y;
			return *this;
		}

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (uint64_t) a * (uint64_t) x;
			_y += t;

			if (_y < t)
				return _y += _F._two_64;
			else
				return _y;
		}

		inline uint64_t& accumulate (const Element &t)
		{
			_y += t;

			if (_y < t)
				return _y += _F._two_64;
			else
				return _y;
		}

		inline uint64_t& accumulate_special (const Element &t)
		{
			return _y += t;
		}

		inline Element &get (Element &y) {
			_y %= (uint64_t) _F._modulus;
			//if ((int64_t) _y < 0) _y += _F._modulus;
			return y = (uint32_t) _y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset() {
			_y = 0;
		}

	private:

		Field _F;
		uint64_t _y;
	};

	//! Specialization of DotProductDomain for unsigned short modular field

	template <>
	class DotProductDomain<Modular<uint8_t> > : private virtual VectorDomainBase<Modular<uint8_t> > {
	public:

		typedef uint8_t Element;

		DotProductDomain (const Modular<uint8_t> &F) :
			VectorDomainBase<Modular<uint8_t> > (F)
		{}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
	};

	//! Specialization of DotProductDomain for unsigned short modular field

	template <>
	class DotProductDomain<Modular<uint16_t> > : private virtual VectorDomainBase<Modular<uint16_t> > {
	public:

		typedef uint16_t Element;

		DotProductDomain (const Modular<uint16_t> &F) :
			VectorDomainBase<Modular<uint16_t> > (F)
		{}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
	};

	//! Specialization of DotProductDomain for uint32_t modular field

	template <>
	class DotProductDomain<Modular<uint32_t> > : private virtual VectorDomainBase<Modular<uint32_t> > {
	public:

		typedef uint32_t Element;

		DotProductDomain (const Modular<uint32_t> &F) :
			VectorDomainBase<Modular<uint32_t> > (F)
		{}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const;

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const;
	};

	//! Specialization of MVProductDomain for uint8_t modular field

	template <>
	class MVProductDomain<Modular<uint8_t> > {
	public:

		typedef uint8_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized (VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint32_t> _tmp;
	};

	//! Specialization of MVProductDomain for uint16_t modular field

	template <>
	class MVProductDomain<Modular<uint16_t> > {
	public:

		typedef uint16_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized (VD, w, A, v, VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint16_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
	};

	//! Specialization of MVProductDomain for uint32_t modular field

	template <>
	class MVProductDomain<Modular<uint32_t> > {
	public:

		typedef uint32_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized (VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<uint32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
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


} // namespace LinBox

#include "linbox/field/Modular/modular.inl"
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

