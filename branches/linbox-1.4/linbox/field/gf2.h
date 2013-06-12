/* linbox/field/gf2.h
 * Copyright (C) 2003-2007 The LinBox group
 *
 * Authors : B. Hovinen, JG Dumas, C. Pernet
 *
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

#ifndef __LINBOX_field_gf2_H
#define __LINBOX_field_gf2_H

#include <iostream>
#include <climits>
#include <cmath>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
#include "linbox/vector/bit-vector.h"
#include "linbox/field/field-traits.h"
// #include "linbox/vector/vector-domain.h"

#ifndef __PATHCC__
#define stdBitReference std::_Bit_reference
#else
#define stdBitReference std::vector<bool>::reference
#endif

// Namespace in which all LinBox code resides
namespace LinBox
{

	class GF2RandIter;

	/**
	 * \brief Integers modulo 2
	 *
	 * This is a tuned implementation of the field of integers modulo
	 * 2. In particular, when one constructs a VectorDomain object over
	 * this field, highly optimized bit operations will be used to make
	 * vector arithmetic very fast.
	 \ingroup field
	 */

	template <class Ring>
	struct ClassifyRing;

	class GF2;

	template<>
	struct ClassifyRing<GF2> {
		typedef RingCategories::ModularTag categoryTag;
	};

	class GF2 : public FieldInterface {
	public:
		const bool zero,one,mOne;


		/** Element type
		*/
		typedef bool Element;

		/** Random iterator generator type.
		 * It must meet the common object interface of random element generators
		 * as given in the the archetype RandIterArchetype.
		 */
		typedef GF2RandIter RandIter;

		/** @name Object Management
		*/
		//@{

		/** Default constructor.
		*/
		GF2 () :
			zero(false),one(true),mOne(true)
		{}
		GF2 (int p, int exp = 1) :
			zero(false),one(true),mOne(true)
		{
			if(p != 2) throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be 2");
			if(exp != 1) throw PreconditionFailed(__func__,__FILE__,__LINE__,"exponent must be 1");
		}

		/** Copy constructor.
		 * Constructs Modular object by copying the field.
		 * This is required to allow field objects to be passed by value
		 * into functions.
		 * @param  F Modular object.
		 */
		GF2 (const GF2 & F ) :
			zero(false),one(true),mOne(true) {}

		/** Assignment operator.
		 * Required by the archetype
		 *
		 * @param F constant reference to Modular object
		 * @return reference to Modular object for self
		 */
		const GF2 &operator = (const GF2 &F)
		{
			return *this;
		}

		/** Initialization of field base element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * This is not a specialization of the template function because
		 * such a specialization is not allowed inside the class declaration.
		 * @return reference to field base element.
		 * @param x field base element to contain output (reference returned).
		 * @param y integer.
		 */
		Element &init (Element &x, const int &y ) const
		{
			return x = y & 1;
		}

		Element &init (Element &x, const unsigned int &y ) const
		{
			return x = y & 1;
		}

		Element &init (Element &x, const long &y ) const
		{
			return x = y & 1;
		}

		Element &init (Element &x, const unsigned long &y ) const
		{
			return x = y & 1;
		}

		Element &init (Element &x, const float &y) const
		{
			return x = static_cast<unsigned char>(y) & 1;
		}

		Element &init (Element &x, const double &y) const
		{
			return x = static_cast<unsigned char>(y) & 1;
		}

		Element &init (Element &x, const integer &y) const
		{
			return x = static_cast<long>(y) & 1;
		}


		Element &init(Element&x) const
		{
			return x = false;
		}

		BitVector::reference init (BitVector::reference x, const integer &y = 0) const
		{
			return x = long (y) & 1;
		}

		stdBitReference init (stdBitReference x, const integer &y = 0) const
		{
			return x = long (y) & 1;
		}

		/** Conversion of field base element to a template class T.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to template class T.
		 * @param x template class T to contain output (reference returned).
		 * @param y constant field base element.
		 */
		integer &convert (integer &x, Element y) const
		{
			return x = y;
		}

		stdBitReference convert (stdBitReference x, Element y) const
		{
			return x = y;
		}

		template<class XXX>
		XXX& convert (XXX& x, Element y) const
		{
			return x = static_cast<XXX>(y);
		}

#if 0
		unsigned int &convert (unsigned int &x, Element y) const
		{
			return x = static_cast<unsigned int>(y);
		}

		int &convert (int &x, Element y) const
		{
			return x = static_cast<int>(y);
		}

		unsigned long &convert (unsigned long &x, Element y) const
		{
			return x = static_cast<unsigned long>(y);
		}

		long &convert (long &x, Element y) const
		{
			return x = static_cast<int>(y);
		}

		float &convert (float &x, Element y) const
		{
			return x = static_cast<float>(y);
		}

		double &convert (double &x, Element y) const
		{
			return x = static_cast<double>(y);
		}
#endif

		/** Assignment of one field base element to another.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &assign (Element &x, Element y) const
		{
			return x = y;
		}

		BitVector::reference assign (BitVector::reference x, Element y) const
		{
			return x = y;
		}

		stdBitReference assign (stdBitReference x, Element y) const
		{
			return x = y;
		}

		/** Cardinality.
		 * Return integer representing cardinality of the domain.
		 * Returns a non-negative integer for all domains with finite
		 * cardinality, and returns -1 to signify a domain of infinite
		 * cardinality.
		 * @return integer representing cardinality of the domain
		 */
		integer &cardinality (integer &c) const
		{
			return c = 2;
		}

		/** Characteristic.
		 * Return integer representing characteristic of the domain.
		 * Returns a positive integer to all domains with finite characteristic,
		 * and returns 0 to signify a domain of infinite characteristic.
		 * @return integer representing characteristic of the domain.
		 */
		integer &characteristic (integer &c) const
		{
			return c = 2;
		}

		//@} Object Management

		/** @name Arithmetic Operations
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field base elements will
		 * give undefined results.
		 */
		//@{

		/** Equality of two elements.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return boolean true if equal, false if not.
		 * @param  x field base element
		 * @param  y field base element
		 */
		bool areEqual (Element x, Element y) const
		{
			return x == y;
		}

		/** Zero equality.
		 * Test if field base element is equal to zero.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals zero, false if not.
		 * @param  x field base element.
		 */
		bool isZero (Element x) const
		{
			return !x;
		}

		/** One equality.
		 * Test if field base element is equal to one.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals one, false if not.
		 * @param  x field base element.
		 */
		bool isOne (Element x) const
		{
			return x;
		}

		//@} Arithmetic Operations

		/** @name Input/Output Operations */
		//@{

		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream &write (std::ostream &os) const
		{
			return os << "integers mod 2";
		}

		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream &read (std::istream &is)
		{
			return is;
		}

		/** Print field base element.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return output stream to which field base element is written.
		 * @param  os  output stream to which field base element is written.
		 * @param  x   field base element.
		 */
		std::ostream &write (std::ostream &os, Element x) const
		{
			return os << x;
		}

		/** Read field base element.
		 * @pre This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return input stream from which field base element is read.
		 * @param  is  input stream from which field base element is read.
		 * @param  x   field base element.
		 */
		std::istream &read (std::istream &is, Element &x) const
		{ is >> x; return is;
		}

		/** Read field base element.
		 * @param is input stream
		 * @param x
		 * @return  \c is
		 */
		std::istream &read (std::istream &is, BitVector::reference x) const
		{ is >> x; return is;
		}

		/** Read field base element.
		 * @param is input stream
		 * @param x
		 * @return  \c is
		 */
		std::istream &read (std::istream &is, stdBitReference x) const
		{ bool a; is >> a; x=a; return is;
		}

		//@}

		/** @name Arithmetic Operations
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field base elements will
		 * give undefined results.
		 */
		//@{

		/** Addition.
		 * x = y + z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &add (Element &x, Element y, Element z) const
		{
			return x = y ^ z;
		}

		/** Addition.
		 * @param x
		 * @param y
		 * @param z
		 */
		BitVector::reference add (BitVector::reference x, Element y, Element z) const
		{
			return x = y ^ z;
		}

		/** Addition.
		 * @param x
		 * @param y
		 * @param z
		 */
		stdBitReference add (stdBitReference x, Element y, Element z) const
		{
			return x = y ^ z;
		}

		/** Subtraction.
		 * x = y - z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &sub (Element &x, Element y, Element z) const
		{
			return x = y ^ z;
		}

		/** Subtraction.
		 * @param x
		 * @param y
		 * @param z
		 */
		BitVector::reference sub (BitVector::reference x, Element y, Element z) const
		{
			return x = y ^ z;
		}

		/** Subtraction.
		 * @param x
		 * @param y
		 * @param z
		 */
		stdBitReference sub (stdBitReference x, Element y, Element z) const
		{
			return x = y ^ z;
		}

		/** Multiplication.
		 * x = y * z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element &mul (Element &x, Element y, Element z) const
		{
			return x = y & z;
		}

		/** Multiplication.
		 * @param x
		 * @param y
		 * @param z
		 */
		BitVector::reference mul (BitVector::reference x, Element y, Element z) const
		{
			return x = y & z;
		}

		/** Multiplication.
		 * @param x
		 * @param y
		 * @param z
		 */
		stdBitReference mul (stdBitReference x, Element y, Element z) const
		{
			return x = y & z;
		}

		/** Division.
		 * x = y / z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 * @bug z is unused
		 */
		Element &div (Element &x, Element y, Element z ) const
		{
			return x = y;
		}

		/** Division.
		 * @param x
		 * @param y
		 * @param z
		 */
		BitVector::reference div (BitVector::reference x, Element y, Element z ) const
		{
			return x = y;
		}

		/** Division.
		 * @param x
		 * @param y
		 * @param z
		 */
		stdBitReference div (stdBitReference x, Element y, Element z ) const
		{
			return x = y;
		}

		/** Additive Inverse (Negation).
		 * x = - y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &neg (Element &x, Element y) const
		{
			return x = y;
		}

		/** Additive Inverse (Negation).
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		BitVector::reference neg (BitVector::reference x, Element y) const
		{
			return x = y;
		}

		/** Additive Inverse (Negation).
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		stdBitReference neg (stdBitReference x, Element y) const
		{
			return x = y;
		}

		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &inv (Element &x, Element y) const
		{
			return x = y;
		}

		/** Multiplicative Inverse.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		BitVector::reference inv (BitVector::reference x, Element y) const
		{
			return x = y;
		}

		/** Multiplicative Inverse.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		stdBitReference inv (stdBitReference x, Element y) const
		{
			return x = y;
		}

		/** Natural AXPY.
		 * r  = a * x + y
		 * This function assumes all field elements have already been
		 * constructed and initialized.
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 * @param  y field element.
		 */
		BitVector::reference axpy (BitVector::reference r,
					   Element a,
					   Element x,
					   Element y) const
		{
			return r = (a & x) ^ y;
		}

		/** Natural AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 * @param  y
		 */
		stdBitReference axpy (stdBitReference r,
					  Element a,
					  Element x,
					  Element y) const
		{
			return r = (a & x) ^ y;
		}

		/** Natural AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 * @param  y
		 */
		Element &axpy (Element &r, Element a, Element x, Element y) const
		{
			return r = (a & x) ^ y;
		}

		//@} Arithmetic Operations

		/** @name Inplace Arithmetic Operations
		 * x <- x op y; x <- op x
		 */
		//@{

		/** Inplace Addition.
		 * x += y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &addin (Element &x, Element y) const
		{
			return x ^= y;
		}

		/** Inplace Addition.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		BitVector::reference addin (BitVector::reference x, Element y) const
		{
			return x ^= y;
		}

		/** Inplace Addition.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		stdBitReference addin (stdBitReference x, Element y) const
		{
			return x = x ^ y;
		}

		/** Inplace Subtraction.
		 * x -= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &subin (Element &x, Element y) const
		{
			return x ^= y;
		}

		/** Inplace Subtraction.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		BitVector::reference subin (BitVector::reference x, Element y) const
		{
			return x ^= y;
		}

		/** Inplace Subtraction.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		stdBitReference subin (stdBitReference x, Element y) const
		{
			return x = x ^ y;
		}

		/** Inplace Multiplication.
		 * x *= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element &mulin (Element &x, Element y) const
		{
			return x &= y;
		}

		/** Inplace Multiplication.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		BitVector::reference mulin (BitVector::reference x, Element y) const
		{
			return x &= y;
		}

		/** Inplace Multiplication.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 */
		stdBitReference mulin (stdBitReference x, Element y) const
		{
			return x = (bool)x & y;
		}

		/** Inplace Division.
		 * x /= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @bug y is unused
		 */
		Element &divin (Element &x, Element y ) const
		{
			return x;
		}

		/** Inplace Division.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 * @bug y is unused
		 */
		BitVector::reference divin (BitVector::reference x, Element y ) const
		{
			return x;
		}

		/** Inplace Division.
		 * @return reference to x.
		 * @param  x
		 * @param  y
		 * @bug y is unused
		 */
		stdBitReference divin (stdBitReference x, Element y ) const
		{
			return x;
		}

		/** Inplace Additive Inverse (Inplace Negation).
		 * x = - x
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element &negin (Element &x) const
		{
			return x;
		}

		/** Inplace Additive Inplace (Inplace Negation).
		 * @return reference to x.
		 * @param  x
		 * @bug y is unused
		 */
		BitVector::reference negin (BitVector::reference x) const
		{
			return x;
		}

		/** Inplace Additive Inplace (Inplace Negation).
		 * @return reference to x.
		 * @param  x
		 * @bug y is unused
		 */
		stdBitReference negin (stdBitReference x) const
		{
			return x;
		}

		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes the field base elementhas already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element &invin (Element &x) const
		{
			return x;
		}

		/** Inplace Multiplicative Inverse.
		 * @return reference to x.
		 * @param  x
		 */
		BitVector::reference invin (BitVector::reference x) const
		{
			return x;
		}

		/** Inplace Multiplicative Inverse.
		 * @return reference to x.
		 * @param  x
		 */
		stdBitReference invin (stdBitReference x) const
		{
			return x;
		}

		/** Inplace AXPY.
		 * r  += a * x
		 * This function assumes all field elements have already been
		 * constructed and initialized.
		 * Purely virtual
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 */
		Element &axpyin (Element &r, Element a, Element x) const
		{
			return r ^= a & x;
		}

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 */
		BitVector::reference axpyin (BitVector::reference r, Element a, Element x) const
		{
			return r ^= a & x;
		}

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 */
		stdBitReference axpyin (stdBitReference r, Element a, Element x) const
		{
			return r = r ^ (a & x);
		}

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 */
		Element &axpyin (Element &r, const stdBitReference a, Element x) const
		{
			return r ^= a & x;
		}

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 */
		stdBitReference axpyin (stdBitReference r, const stdBitReference a, Element x) const
		{
			return r = r ^ (a & x);
		}

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 */
		Element &axpyin (Element &r, Element a, const stdBitReference x) const
		{
			return r ^= a & static_cast<bool>(x);
		}

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 */
		stdBitReference axpyin (stdBitReference r, Element a, const stdBitReference x) const
		{
			return r = r ^ (a & static_cast<bool>(x));
		}

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 */
		Element &axpyin (Element &r, const stdBitReference a, const stdBitReference x) const
		{
			return r ^= a & static_cast<bool>(x);
		}

		/** Inplace AXPY.
		 * @return reference to r.
		 * @param  r
		 * @param  a
		 * @param  x
		 */
		stdBitReference axpyin (stdBitReference r, const stdBitReference a, const stdBitReference x) const
		{
			return r = r ^ (a & static_cast<bool>(x));
		}

		//@} Inplace Arithmetic Operations

		static inline int getMaxModulus()
		{
			return 2;
		}

	}; // class GF2

} // namespace LinBox

// #define LINBOX_field_gf2_H
// #include "linbox/vector/vector-domain.h"


// Specialization of homomorphism for basefield
#include "linbox/randiter/gf2.h"

// #include <bits/stl_bvector.h>
namespace std
{
	//! @todo JGD 05.11.2009 : it should be in bits/stl_bvector.h  ...
	inline void swap(stdBitReference __x, stdBitReference __y)
	{
		bool __tmp = __x;
		__x = __y;
		__y = __tmp;
	}
}


#include "linbox/field/gf2.inl"

#endif // __LINBOX_field_gf2_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

