/* -*- mode: c; style: linux -*- */

/* linbox/src/field/unparametric.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
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
 
#ifndef __FIELD_UNPARAMETRIC_H
#define __FIELD_UNPARAMETRIC_H

#include <string>
#include <algorithm>

#include "linbox/integer.h"
#include "linbox/randiter/unparametric.h"

namespace LinBox
{
	/** Unparameterized field template.
	 * Implements LinBox field common object interface for unparameterized 
	 * fields.
	 * Used to generate efficient field classes for unparameterized fields.
	 * Constructs LinBox unparameterized fields from field types K.
	 * In particular, constructs LinBox fields for
	 * unparameterized fields from field types that
	 * adhere to the operations for double, for
	 * example unparam_field< float >.
	 * Can be used as a pattern to write a particular
	 * field interface, such as, unparam_field< SaclibQ > as
	 * a template specialization.
	 * @param  K unparameterized field class
	 */
	template <class K>
	class UnparametricField
	{
	    public:
    
		/** @name Common Object Interface for a LinBox Field.
		 * These methods are required of all LinBox fields.
		 */
		//@{
    
		/** Field element type.
		 * The field element must contain a default constructor, 
		 * a copy constructor, a destructor, and an assignment operator.
		 */
		typedef K Element;    

		/// Random field element generator type.
		typedef UnparametricRandIter<K> RandIter;

		/** @name Object Management.
		 * x <- convert (y)
		 */
		//@{
    
		/** Copy constructor.
		 * Constructs UnparametricField object by copying the field.
		 * This is required to allow field objects to be passed by value
		 * into functions.
		 * @param  F UnparametricField object.
		 */
		UnparametricField (const UnparametricField &F) {}
    
		/** Destructor.
		 * This destructs the field object, but it does not destroy the field 
		 * element objects.  The destructor for each field element must also 
		 * be called.
		 * _elem_ptr points.
		 */
		~UnparametricField () {}
    
		/** Assignment operator.
		 * Assigns UnparametricField object F to field.
		 * @param  F UnparametricField object.
		 */
		UnparametricField &operator=(const UnparametricField &F) { return *this; }
    
		/** Initialization of field element from an integer.
		 * Behaves like C++ allocator construct.
		 * This function assumes the output field element x has already been 
		 * constructed, but that it is not already initialized.
		 * For now, this is done by converting the integer type to a C++ 
		 * long and then to the element type through the use of static casts.
		 * This, of course, assumes such static casts are possible.
		 * This function should be changed in the future to avoid using long.
		 * @return reference to field element.
		 * @param x field element to contain output (reference returned).
		 * @param y integer.
		 */
		Element &init (Element &x, const Integer &y=0) const 
			{ return x = static_cast<const element&> (static_cast<const long&> (y)); }
    
		/** Conversion of field element to an integer.
		 * This function assumes the output field element x has already been 
		 * constructed, but that it is not already initialized.
		 * For now, this is done by converting the element type to a C++ 
		 * long and then to the integer type through the use of static casts.
		 * This, of course, assumes such static casts are possible.
		 * This function should be changed in the future to avoid using long.
		 * @return reference to integer.
		 * @param x reference to integer to contain output (reference returned).
		 * @param y constant reference to field element.
		 */
		Integer &convert (Integer &x, const Element &y) const 
		{ 
			Element temp (y);
			return x = static_cast<long> (temp); 
		}
    
		/** Assignment of one field element to another.
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 * @return reference to x
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &assign (Element &x, const Element &y) const { return x = y; }
    
		/** Cardinality.
		 * Return integer representing cardinality of the field.
		 * Returns a non-negative integer for all fields with finite
		 * cardinality, and returns -1 to signify a field of infinite 
		 * cardinality.
		 * The default behavior of this function is to return -1 to signify
		 * a field of infinite cardinality.  This should be changed via
		 * partial template specialization for fields of other cardinalities.
		 * @return integer representing cardinality of the field
		 */
		Integer &cardinality (Integer &c) const { return c = -1; }
    
		/** Characteristic.
		 * Return integer representing characteristic of the field.
		 * Returns a positive integer to all fields with finite characteristic,
		 * and returns 0 to signify a field of infinite characteristic.
		 * The default behavior of this function is to return 0 to signify
		 * a field of infinite characteristic.  This should be changed via
		 * partial template specialization for fields of other characteristics.
		 * @return integer representing characteristic of the field.
		 */
		Integer &characteristic (Integer &c) const { return c = 0; }
    
		//@} Object Management
    
		/** @name Arithmetic Operations 
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field elements will
		 * give undefined results.
		 */
		//@{
    
		/** Equality of two elements.
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 * @return boolean true if equal, false if not.
		 * @param  x field element
		 * @param  y field element
		 */
		bool areEqual (const Element &x, const Element &y) const { return x == y; }
    
		/** Addition.
		 * x = y + z
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 * @param  z field element.
		 */
		Element &add (Element &x, const Element &y, const Element &z) const
			{ return x = y + z; }
    
		/** Subtraction.
		 * x = y - z
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 * @param  z field element.
		 */
		Element &sub (Element &x, const Element &y, const Element &z) const
			{ return x = y - z; }
    
		/** Multiplication.
		 * x = y * z
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 * @param  z field element.
		 */
		Element &mul (Element &x, const Element &y, const Element &z) const
			{ return x = y * z; }
    
		/** Division.
		 * x = y / z
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 * @param  z field element.
		 */
		Element &div (Element &x, const Element &y, const Element &z) const
			{ return x = y / z; }
    
		/** Additive Inverse (Negation).
		 * x = - y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &neg (Element &x, const Element &y) const { return x = - y; }
    
		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &inv (Element &x, const Element &y) const 
			{ return x = Element (1) / y; }
    
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
		Element &axpy (Element &r, 
			       const Element &a, 
			       const Element &x, 
			       const Element &y) const
			{ return r = a * x + y; }
 
		//@} Arithmetic Operations
    
		/** @name Inplace Arithmetic Operations 
		 * x <- x op y; x <- op x
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field elements will
		 * give undefined results.
		 */
		//@{
    
		/** Zero equality.
		 * Test if field element is equal to zero.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * @return boolean true if equals zero, false if not.
		 * @param  x field element.
		 */
		bool isZero (const Element &x) const { return x == element (0); }
    
		/** One equality.
		 * Test if field element is equal to one.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * @return boolean true if equals one, false if not.
		 * @param  x field element.
		 */
		bool isOne (const Element &x) const { return x == element (1); }
    
		/** Inplace Addition.
		 * x += y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &addin (Element &x, const Element &y) const { return x += y; }
    
		/** Inplace Subtraction.
		 * x -= y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &subin (Element &x, const Element &y) const { return x -= y; }
    
		/** Inplace Multiplication.
		 * x *= y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &mulin (Element &x, const Element &y) const { return x *= y; }
    
		/** Inplace Division.
		 * x /= y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &divin (Element &x, const Element &y) const { return x /= y; }
    
		/** Inplace Additive Inverse (Inplace Negation).
		 * x = - x
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 */
		Element &negin (Element &x) const { return x = - x; }
    
		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes the field elementhas already been 
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 */
		Element &invin (Element &x) const { return x = element (1) / x; }
    
		/** Inplace AXPY.
		 * r  += a * x
		 * This function assumes all field elements have already been 
		 * constructed and initialized.
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 */
		Element &axpyin (Element &r, const Element &a, const Element &x) const
			{ return r += a * x; }
 
		//@} Inplace Arithmetic Operations
    
		/** @name Input/Output Operations */
		//@{
    
		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		ostream &write (ostream &os) const { return os << "unparamterized field"; }
    
		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		istream &read (istream &is) const { return is; }
    
		/** Print field element.
		 * @return output stream to which field element is written.
		 * @param  os  output stream to which field element is written.
		 * @param  x   field element.
		 */
		ostream &write (ostream &os, const Element &x) const { return os << x; }
    
		/** Read field element.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		istream &read (istream &is, Element &x) const { return is >> x; }
    
		//@}
    
		//@} Common Object Interface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all LinBox fields
		 * and are included only for the implementation of this field
		 * template.
		 */
		//@{
    
		/// Default constructor
		UnparametricField (void) {}
    
		/** Constructor from field object.
		 * @param  A unparameterized field object
		 */
		UnparametricField (const K &A) {} 
    
		/** Constant access operator.
		 * @return constant reference to field object
		 */
		const K &operator () (void) const { return Element (); }
    
		/** Access operator.
		 * @return reference to field object
		 */
		K &operator () (void) { return Element (); }
    
		//@} Implementation-Specific Methods

	}; // template <class K> class UnparametricField
  
} // namespace LinBox

#include "linbox/randiter/unparametric.h"

#endif // __FIELD_UNPARAMETRIC_H_
