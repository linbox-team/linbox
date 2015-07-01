/* -*- mode: c; style: linux -*- */

/* linbox/src/field/abstract.h
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

#ifndef __FIELD_ENVELOPE_H
#define __FIELD_ENVELOPE_H

#include <iostream>

#include "linbox/integer.h"
#include "linbox/element/envelope.h"
#include "linbox/field/abstract.h"
#include "linbox/element/abstract.h"
#include "linbox/randiter/abstract.h"
#include "linbox/randiter/envelope.h"

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
	// Forward declarations
	template <class Field> class RandIter_envelope;

	/** Field Envelope Template.
	 * Derived class used to implement the field archetype to minimize
	 * code bloat.  This class implements all purely virtual member functions
	 * of the abstract base class.  This class is used to wrap a
	 * \Ref{LinBox}
	 * field so that it might be used with the Field archetype.
	 */
	template <class Field>
	class Field_envelope : public Field_abstract
	{
	    public:

		/** Element type.
		 * It is derived from the class Element_abstract, and it must contain
		 * a wrapped field element.
		 */
		typedef Element_envelope<Field> element;

		/** Random iterator generator type.
		 * It is derived from the class RandIter_abstract, and it must contain
		 * a wrapped field random iterator generator.
		 */
		typedef RandIter_envelope<Field> randIter;

		/** @name Object Management
		 */
		//@{
 
		/** Default constructor.
		 * In this implementation, this means copying the field E._field.
		 */
		Field_envelope (void) {}

		/** Constructor from field to be wrapped.
		 * @param F Field object to be wrapped.
		 */
		Field_envelope (const Field& F) : _field (F) {}
 
		/** Copy constructor.
		 * Constructs Field_envelope object by copying the field.
		 * This is required to allow field objects to be passed by value
		 * into functions.
		 * In this implementation, this means copying the field E._field.
		 * @param  E Field_envelope object.
		 */
		Field_envelope (const Field_envelope& E) : _field (E._field) {}
 
		/** Virtual copy constructor.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * This function is not part of the common object interface.
		 * @return pointer to new object in dynamic memory.
		 */
		Field_abstract* clone () const
			{ return new Field_envelope (*this); }

		/** Assignment operator.
		 * Required by abstract base class.
		 * @return reference to Field_abstract object for self
		 * @param F constant reference to Field_abstract object
		 */
		Field_abstract& operator= (const Field_abstract& F)
		{
			if (this != &F) // guard against self-assignment
				_field = static_cast<const Field_envelope&> (F)._field;

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
		Element_abstract& init (Element_abstract& x, const Integer& y = 0) const
		{
			_field.init (static_cast<Element_envelope<Field>&> (x)._elem, y);
			return x;
		}
 
		/** Conversion of field base element to a template class T.
		 * This function assumes the output field base element x has already been
		 * constructed, but that it is not already initialized.
		 * @return reference to template class T.
		 * @param x template class T to contain output (reference returned).
		 * @param y constant field base element.
		 */
		Integer& convert (Integer& x, const Element_abstract& y) const
		{
			_field.convert (x, static_cast<const Element_envelope<Field>&> (y)._elem);
			return x;
		}
 
		/** Assignment of one field base element to another.
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element_abstract& assign (Element_abstract& x, const Element_abstract& y) const
		{
			_field.assign (static_cast<Element_envelope<Field>&> (x)._elem,
				       static_cast<const Element_envelope<Field>&> (y)._elem);
			return x;
		}

		/** Cardinality.
		 * Return integer representing cardinality of the domain.
		 * Returns a non-negative integer for all domains with finite
		 * cardinality, and returns -1 to signify a domain of infinite
		 * cardinality.
		 * @return integer representing cardinality of the domain
		 */
		Integer& cardinality (Integer& c) const
			{ return _field.cardinality (c); }
 
		/** Characteristic.
		 * Return integer representing characteristic of the domain.
		 * Returns a positive integer to all domains with finite characteristic,
		 * and returns 0 to signify a domain of infinite characteristic.
		 * @return integer representing characteristic of the domain.
		 */
		Integer& characteristic (Integer& c) const
			{ return _field.characteristic (c); }

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
		bool areEqual (const Element_abstract& x, const Element_abstract& y) const
		{
			return _field.areEqual (static_cast<const Element_envelope<Field>&> (x)._elem,
						static_cast<const Element_envelope<Field>&> (y)._elem);
		}

		/** Addition.
		 * x = y + z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element_abstract& add (Element_abstract& x,
				       const Element_abstract& y,
				       const Element_abstract& z) const
		{
			_field.add (static_cast<Element_envelope<Field>&> (x)._elem,
				    static_cast<const Element_envelope<Field>&> (y)._elem,
				    static_cast<const Element_envelope<Field>&> (z)._elem);
			return x;
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
		Element_abstract& sub (Element_abstract& x,
				       const Element_abstract& y,
				       const Element_abstract& z) const
		{
			_field.sub (static_cast<Element_envelope<Field>&> (x)._elem,
				    static_cast<const Element_envelope<Field>&> (y)._elem,
				    static_cast<const Element_envelope<Field>&> (z)._elem);
			return x;
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
		Element_abstract& mul (Element_abstract& x,
				       const Element_abstract& y,
				       const Element_abstract& z) const
		{
			_field.mul (static_cast<Element_envelope<Field>&> (x)._elem,
				    static_cast<const Element_envelope<Field>&> (y)._elem,
				    static_cast<const Element_envelope<Field>&> (z)._elem);
			return x;
		}
 
		/** Division.
		 * x = y / z
		 * This function assumes all the field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 * @param  z field base element.
		 */
		Element_abstract& div (Element_abstract& x,
				       const Element_abstract& y,
				       const Element_abstract& z) const
		{
			_field.div (static_cast<Element_envelope<Field>&> (x)._elem,
				    static_cast<const Element_envelope<Field>&> (y)._elem,
				    static_cast<const Element_envelope<Field>&> (z)._elem);
			return x;
		}
 
		/** Additive Inverse (Negation).
		 * x = - y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element_abstract& neg (Element_abstract& x, const Element_abstract& y) const
		{
			_field.neg (static_cast<Element_envelope<Field>&> (x)._elem,
				    static_cast<const Element_envelope<Field>&> (y)._elem);
			return x;
		}
 
		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element_abstract& inv (Element_abstract& x, const Element_abstract& y) const
		{
			_field.inv (static_cast<Element_envelope<Field>&> (x)._elem,
				    static_cast<const Element_envelope<Field>&> (y)._elem);
			return x;
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
		Element_abstract& axpy (Element_abstract& r, 
					const Element_abstract& a, 
					const Element_abstract& x, 
					const Element_abstract& y) const
		{
			_field.axpy (static_cast<Element_envelope<Field>&> (r)._elem,
				     static_cast<const Element_envelope<Field>&> (a)._elem,
				     static_cast<const Element_envelope<Field>&> (x)._elem,
				     static_cast<const Element_envelope<Field>&> (y)._elem);
			return r;
		}
 
		//@} Arithmetic Operations
 
		/** @name Inplace Arithmetic Operations
		 * x <- x op y; x <- op x
		 */
		//@{

		/** Zero equality.
		 * Test if field base element is equal to zero.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals zero, false if not.
		 * @param  x field base element.
		 */
		bool isZero (const Element_abstract& x) const
			{ return _field.isZero (static_cast<const Element_envelope<Field>&> (x)._elem); }
 
		/** One equality.
		 * Test if field base element is equal to one.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return boolean true if equals one, false if not.
		 * @param  x field base element.
		 */
		bool isOne (const Element_abstract& x) const
			{ return _field.isOne (static_cast<const Element_envelope<Field>&> (x)._elem); }

		/** Inplace Addition.
		 * x += y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element_abstract& addin (Element_abstract& x, const Element_abstract& y) const
		{
			_field.addin (static_cast<Element_envelope<Field>&> (x)._elem,
				      static_cast<const Element_envelope<Field>&> (y)._elem);
			return x;
		}
 
		/** Inplace Subtraction.
		 * x -= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element_abstract& subin (Element_abstract& x, const Element_abstract& y) const
		{
			_field.subin (static_cast<Element_envelope<Field>&> (x)._elem,
				      static_cast<const Element_envelope<Field>&> (y)._elem);
			return x;
		}
 
		/** Inplace Multiplication.
		 * x *= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element_abstract& mulin (Element_abstract& x, const Element_abstract& y) const
		{
			_field.mulin (static_cast<Element_envelope<Field>&> (x)._elem,
				      static_cast<const Element_envelope<Field>&> (y)._elem);
			return x;
		}

		/** Inplace Division.
		 * x /= y
		 * This function assumes both field base elements have already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 * @param  y field base element.
		 */
		Element_abstract& divin (Element_abstract& x, 
					 const Element_abstract& y) const
		{
			_field.divin (static_cast<Element_envelope<Field>&> (x)._elem,
				      static_cast<const Element_envelope<Field>&> (y)._elem);
			return x;
		}
 
		/** Inplace Additive Inverse (Inplace Negation).
		 * x = - x
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element_abstract& negin (Element_abstract& x) const
		{
			_field.negin (static_cast<Element_envelope<Field>&> (x)._elem);
			return x;
		}
 
		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes the field base elementhas already been
		 * constructed and initialized.
		 * @return reference to x.
		 * @param  x field base element (reference returned).
		 */
		Element_abstract& invin (Element_abstract& x) const
		{
			_field.invin (static_cast<Element_envelope<Field>&> (x)._elem);
			return x;
		}

		/** Inplace AXPY.
		 * r  += a * x
		 * This function assumes all field elements have already been 
		 * constructed and initialized.
		 * @return reference to r.
		 * @param  r field element (reference returned).
		 * @param  a field element.
		 * @param  x field element.
		 */
		Element_abstract& axpyin (Element_abstract& r, 
					  const Element_abstract& a, 
					  const Element_abstract& x) const
		{
			_field.axpyin (static_cast<Element_envelope<Field>&> (r)._elem,
				       static_cast<const Element_envelope<Field>&> (a)._elem,
				       static_cast<const Element_envelope<Field>&> (x)._elem);
			return r;
		}
 
		//@} Inplace Arithmetic Operations

		/** @name Input/Output Operations */
		//@{

		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		ostream& write (ostream& os) const { return _field.write (os); }
 
		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		istream& read (istream& is) { return _field.read (is); }

		/** Print field base element.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return output stream to which field base element is written.
		 * @param  os  output stream to which field base element is written.
		 * @param  x   field base element.
		 */
		ostream& write (ostream& os, const Element_abstract& x) const
			{ return _field.write (os, static_cast<const Element_envelope<Field>&> (x)._elem); }
 
		/** Read field base element.
		 * This function assumes the field base element has already been
		 * constructed and initialized.
		 * @return input stream from which field base element is read.
		 * @param  is  input stream from which field base element is read.
		 * @param  x   field base element.
		 */
		istream& read (istream& is, Element_abstract& x) const
			{ return _field.read (is, static_cast<Element_envelope<Field>&> (x)._elem); }

		//@}

	    private:

		friend class RandIter_envelope<Field>;

		/// Wrapped field.
		Field _field;

	}; // class Field_envelope

} // namespace LinBox

#include "linbox/randiter/envelope.h"

#endif // __FIELD_ENVELOPE_H
