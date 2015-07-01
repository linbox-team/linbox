/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/gf2.h
 * Copyright (C) 2003 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Evolved from modular.h by Bradford Hovinen
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __FIELD_GF2_H
#define __FIELD_GF2_H

#include <iostream>
#include <climits>
#include <cmath>

#include "linbox/integer.h"
#include "linbox/field/field-interface.h"
#include "linbox/vector/bit-vector.h"
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <string>
using std::istream;
using std::ostream;
using std::string;

#endif

// Namespace in which all LinBox code resides
namespace LinBox 
{

class GF2RandIter;

/** @name GF2
 * @memo Integers modulo 2
 *
 * @doc
 * This is a tuned implementation of the field of integers modulo
 * 2. In particular, when one constructs a VectorDomain object over
 * this field, highly optimized bit operations will be used to make
 * vector arithmetic very fast.
 */

class GF2 : public FieldInterface
{
    public:

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
	GF2 () {}

	/** Copy constructor.
	 * Constructs Modular object by copying the field.
	 * This is required to allow field objects to be passed by value
	 * into functions.
	 * @param  F Modular object.
	 */
	GF2 (const GF2 &F) {}

#ifdef __LINBOX_XMLENABLED
	// XML Reader Constructor
	GF2(Reader &R) 
	{
		int m;

		if(!R.expectTagName("field") || !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("finite") || !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("characteristic") || !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagNum(m)) return;
		if(m != 2) {
			R.setErrorString("Tried to initalize field GF2, but got wrong characteristic");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		R.upToParent();
		
		R.upToParent();
		if(R.getNextChild()) {
			R.traverseChild();
			if(!R.expectTagName("extension") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagNum(m)) return;
			if(m > 1) {
				R.setErrorString("Could not create GF2, got wrong extension degree.");
				R.setErrorCode(Reader::OTHER);
				return;
			}
			R.upToParent();
			R.upToParent();
			R.getPrevChild();
		}
		R.upToParent();

		return;
	}
 

#endif
	/** Assignment operator
	 * Required by the archetype
	 *
	 * @param F constant reference to Modular object
	 * @return reference to Modular object for self
	 */
	const GF2 &operator = (const GF2 &F) 
		{ return *this; }

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
	Element &init (Element &x, const integer &y = 0) const
		{ return x = long (y) & 1; }

	BitVector::reference init (BitVector::reference x, const integer &y = 0) const
		{ return x = long (y) & 1; }

	/** Conversion of field base element to a template class T.
	 * This function assumes the output field base element x has already been
	 * constructed, but that it is not already initialized.
	 * @return reference to template class T.
	 * @param x template class T to contain output (reference returned).
	 * @param y constant field base element.
	 */
	integer &convert (integer &x, Element y) const
		{ return x = y; }
 
	/** Assignment of one field base element to another.
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &assign (Element &x, Element y) const
		{ return x = y; }

	BitVector::reference assign (BitVector::reference x, Element y) const
		{ return x = y; }

	/** Cardinality.
	 * Return integer representing cardinality of the domain.
	 * Returns a non-negative integer for all domains with finite
	 * cardinality, and returns -1 to signify a domain of infinite
	 * cardinality.
	 * @return integer representing cardinality of the domain
	 */
	integer &cardinality (integer &c) const
		{ return c = 2; }

	/** Characteristic.
	 * Return integer representing characteristic of the domain.
	 * Returns a positive integer to all domains with finite characteristic,
	 * and returns 0 to signify a domain of infinite characteristic.
	 * @return integer representing characteristic of the domain.
	 */
	integer &characteristic (integer &c) const
		{ return c = 2; }

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
		{ return x == y; }

	/** Zero equality.
	 * Test if field base element is equal to zero.
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return boolean true if equals zero, false if not.
	 * @param  x field base element.
	 */
	bool isZero (Element x) const
		{ return !x; }
 
	/** One equality.
	 * Test if field base element is equal to one.
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return boolean true if equals one, false if not.
	 * @param  x field base element.
	 */
	bool isOne (Element x) const
		{ return x; }

	//@} Arithmetic Operations

#ifndef __LINBOX_XMLENABLED
	/** @name Input/Output Operations */
	//@{

	/** Print field.
	 * @return output stream to which field is written.
	 * @param  os  output stream to which field is written.
	 */
	std::ostream &write (std::ostream &os) const 
		{ return os << "integers mod 2"; }

	/** Read field.
	 * @return input stream from which field is read.
	 * @param  is  input stream from which field is read.
	 */
	std::istream &read (std::istream &is)
		{ return is; }

	/** Print field base element.
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return output stream to which field base element is written.
	 * @param  os  output stream to which field base element is written.
	 * @param  x   field base element.
	 */
	std::ostream &write (std::ostream &os, Element x) const
		{ return os << x; }
 
	/** Read field base element.
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return input stream from which field base element is read.
	 * @param  is  input stream from which field base element is read.
	 * @param  x   field base element.
	 */
	std::istream &read (std::istream &is, Element &x) const
		{ is >> x; return is; }

	std::istream &read (std::istream &is, BitVector::reference x) const
		{ is >> x; return is; }

	//@}

#else
	ostream &write(ostream &os) const
	{
		Writer W;
		if( toTag(W))
			W.write(os);

		return os;
	}

	bool toTag(Writer &W) const
	{

		W.setTagName("field");
		W.setAttribute("implDetail", "gf2");
		W.setAttribute("cardinality", "2");

		W.addTagChild();
		W.setTagName("finite");

		W.addTagChild();
		W.setTagName("characteristic");
		W.addNum(2);
		W.upToParent();

		W.upToParent();

		return true;
	}

	ostream &write(ostream &os, const Element &e) const 
	{
		Writer W;
		if( toTag(W, e))
			W.write(os);

		return os;
	}

	bool toTag(Writer &W, const Element &e) const
	{
		// the manual method
		W.setTagName("cn");
		if(e) 
			W.addDataChild("1");
		else
			W.addDataChild("0");

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
		int m;
		if(!R.expectTagNum(m)) return false;
		if(m == 0) 
			e = false;
		else if(m == 1) 
			e = true;
		else {
			R.setErrorString("Tried to initalize element of GF2, but element was out of range");
			R.setErrorCode(Reader::OTHER);
			return false;
		}

		return true;
	}

#endif


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
		{ return x = y ^ z; }

	BitVector::reference add (BitVector::reference x, Element y, Element z) const
		{ return x = y ^ z; }
 
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
		{ return x = y ^ z; }

	BitVector::reference sub (BitVector::reference x, Element y, Element z) const
		{ return x = y ^ z; }
 
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
		{ return x = y & z; }

	BitVector::reference mul (BitVector::reference x, Element y, Element z) const
		{ return x = y & z; }
 
	/** Division.
	 * x = y / z
	 * This function assumes all the field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 * @param  z field base element.
	 */
	Element &div (Element &x, Element y, Element z) const
		{ return x = y; }

	BitVector::reference div (BitVector::reference x, Element y, Element z) const
		{ return x = y; }
 
	/** Additive Inverse (Negation).
	 * x = - y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &neg (Element &x, Element y) const
		{ return x = y; }

	BitVector::reference neg (BitVector::reference x, Element y) const
		{ return x = y; }
 
	/** Multiplicative Inverse.
	 * x = 1 / y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &inv (Element &x, Element y) const
		{ return x = y; }

	BitVector::reference inv (BitVector::reference x, Element y) const
		{ return x = y; }

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
		{ return r = (a & x) ^ y; }

	Element &axpy (Element &r, Element a, Element x, Element y) const
		{ return r = (a & x) ^ y; }

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
		{ return x ^= y; }

	BitVector::reference addin (BitVector::reference x, Element y) const
		{ return x ^= y; }
 
	/** Inplace Subtraction.
	 * x -= y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &subin (Element &x, Element y) const
		{ return x ^= y; }

	BitVector::reference subin (BitVector::reference x, Element y) const
		{ return x ^= y; }
 
	/** Inplace Multiplication.
	 * x *= y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &mulin (Element &x, Element y) const
		{ return x &= y; }

	BitVector::reference mulin (BitVector::reference x, Element y) const
		{ return x &= y; }
 
	/** Inplace Division.
	 * x /= y
	 * This function assumes both field base elements have already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 * @param  y field base element.
	 */
	Element &divin (Element &x, Element y) const
		{ return x; }

	BitVector::reference divin (BitVector::reference x, Element y) const
		{ return x; }
 
	/** Inplace Additive Inverse (Inplace Negation).
	 * x = - x
	 * This function assumes the field base element has already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 */
	Element &negin (Element &x) const
		{ return x; }

	BitVector::reference negin (BitVector::reference x) const
		{ return x; }
 
	/** Inplace Multiplicative Inverse.
	 * x = 1 / x
	 * This function assumes the field base elementhas already been
	 * constructed and initialized.
	 * @return reference to x.
	 * @param  x field base element (reference returned).
	 */
	Element &invin (Element &x) const
		{ return x; }

	BitVector::reference invin (BitVector::reference x) const
		{ return x; }

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
		{ return r ^= a & x; }

	BitVector::reference axpyin (BitVector::reference r, Element a, Element x) const
		{ return r ^= a & x; }

	//@} Inplace Arithmetic Operations

}; // class GF2

} // namespace LinBox

#include "linbox/randiter/gf2.h"
#include "linbox/field/gf2.inl"

#endif // __FIELD_GF2_H
