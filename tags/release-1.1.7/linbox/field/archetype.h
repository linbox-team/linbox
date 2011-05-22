/* linbox/field/archetype.h
 * Copyright (C) 1999-2005 William J Turner,
 *               2001 Bradford Hovinen
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
 *
 * ------------------------------------
 * 2002-05-14 William J. Turner <wjturner@acm.org>
 * 
 * changed randIter to RandIter.
 * ------------------------------------
 * 2005-06-24 William J. Turner <wjturner@acm.org>
 * 
 * Removed using declarations.
 * ------------------------------------
 */


#ifndef __LINBOX_field_archetype_H
#define __LINBOX_field_archetype_H

#include <iostream>
#include "linbox/field/field-interface.h"
#include "linbox/field/abstract.h"
#include "linbox/field/envelope.h"
#include "linbox/element/archetype.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/randiter/abstract.h"
#include "linbox/randiter/envelope.h"
#include "linbox/integer.h"
#include "linbox/linbox-config.h"
#include "linbox/util/error.h"

namespace LinBox
{
	// Forward declarations
	class RandIterArchetype;
	/** \brief field specification and archetypical instance.
	\ingroup field
	 *
	 * The %FieldArchetype and its encapsulated
	 * element class contain pointers to the \ref FieldAbstract 
	 * and its encapsulated field element, respectively.
	 * %FieldAbstract  then uses virtual member functions to
	 * define operations on its encapsulated field element.  This field 
	 * element has no knowledge of the field properties being used on it 
	 * which means the field object must supply these operations.
	 *
	 * It does not contain elements zero and one because they can be created 
	 * whenever necessary, although it might be beneficial from an efficiency
	 * stand point to include them.  However, because of archetype use three,
	 * the elements themselves cannot be contained, but rather pointers to them.
	 */
	class FieldArchetype : public FieldInterface
	{
	    public:

		/** @name Common Object Interface for a LinBox Field.
		 * These methods are required of all \ref{LinBox} fields.
		 */
		//@{
    
		/// the type in which field elements are represented.
		typedef ElementArchetype Element;

		/// An object of this type is a generator of random field elements.
		typedef RandIterArchetype RandIter;
    
		/// @name Object Management
		//@{
    
		/** \brief Copy constructor.
		 *
		 * Each field class is expected to provide a copy constructor.
		 * This is required to allow field objects to be passed by value into functions.
		 *
		 * In this archetype implementation, this means copying the
		 * field to which <tt> F._field_ptr</tt> points, the
		 * element to which <tt> F._elem_ptr</tt> points, and the
		 * random element generator to which
		 * <tt> F._randIter_ptr</tt> points.
		 */
		FieldArchetype (const FieldArchetype &F) 
		{ 
			if (F._field_ptr != 0) _field_ptr = F._field_ptr->clone (); 
			if (F._elem_ptr != 0) _elem_ptr = F._elem_ptr->clone ();
			if (F._randIter_ptr != 0) _randIter_ptr = F._randIter_ptr->clone ();
		}

		/** \brief Destructor.
		 *
		 * This destroys the field object, but it does not
		 * destroy any field element objects.  
		 *
		 * In this archetype implementation, destruction is deletion of
		 * the field object to which <tt> _field_ptr</tt>
		 * points, the field element to which <tt>
		 * _elem_ptr</tt> points, and the random element
		 * generator to which <tt> _randIter_ptr</tt> points.
		 */
		~FieldArchetype (void) 
		{
			if (_field_ptr != 0) delete _field_ptr;
			if (_elem_ptr != 0) delete _elem_ptr; 
			if (_randIter_ptr != 0) delete _randIter_ptr;
		}
    
		/** \brief Assignment operator.
		 *
		 * In this archetype implementation, this means copying the field
		 * to which <tt> F._field_ptr</tt> points, the element to which
		 * <tt> F._elem_ptr</tt> points, and the random element
		 * generator to which <tt> F._randIter_ptr</tt> points.
		 *
		 * @param F <tt> FieldArchetype</tt> object.
		 */
		FieldArchetype &operator=(const FieldArchetype &F)
		{
			if (this != &F) { // guard against self-assignment
				if (_field_ptr != 0) delete _field_ptr;
				if (_elem_ptr != 0) delete _elem_ptr;
				if (_randIter_ptr != 0) delete _randIter_ptr;
				if (F._field_ptr != 0) _field_ptr = F._field_ptr->clone (); 
				if (F._elem_ptr != 0) _elem_ptr = F._elem_ptr->clone ();
				if (F._randIter_ptr != 0) _randIter_ptr = F._randIter_ptr->clone ();
			}

			return *this;
		}
    
		/** \brief Initialization of field element from an integer.
		 *
		 * x becomes the image of n under the natural map from the integers
		 * to the prime subfield.  It is the result obtained from adding n 1's 
		 * in the field.

		 * This function assumes the output field element x
		 * has already been constructed, but that it is not
		 * necessarily already initialized. In this
		 * archetype implementation, this means the <tt> _elem_ptr</tt> of
		 * x exists, but that it may be the null pointer.
		 *
		 * @return reference to x.
		 * @param x output field element.
		 * @param n input integer.
		 */
		Element &init (Element &x, const integer &n = 0 ) const
		{
			if (x._elem_ptr != 0) delete x._elem_ptr;
			x._elem_ptr = _elem_ptr->clone ();
			_field_ptr->init (*x._elem_ptr, n);
			return x;
		}
  
		/** \brief Conversion of field element to an integer.
		 *
		The meaning of conversion is specific to each field class.
		However, if x is in the prime subfield, the integer n returned is such 
		that an init from n will reproduce x.  Most often, 0 &leq; n &lt; characteristic.

		 *
		 * @return reference to n.
		 * @param n output integer.
		 * @param x input field element.
		 */
		integer &convert (integer &n, const Element &y = 0) const
		{
			_field_ptr->convert (n, *y._elem_ptr);
			return n;
		}
    
		/** \brief  Assignment of one field element to another.
		 *
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this archetype implementation, this means for both x and
		 * y, <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x
		 * @param  x destination field element.
		 * @param  y source field element.
		 */
		Element &assign (Element &x, const Element &y) const
		{
			if (x._elem_ptr == 0) 
				x._elem_ptr = _elem_ptr->clone ();

			_field_ptr->assign (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
    
		/** \brief Cardinality.
		 *
		 * Return c, integer representing cardinality of the field.
		 * c becomes a non-negative integer for all fields with finite
		 * cardinality, and -1 to signify a field of infinite cardinality.
		 */
		integer &cardinality (integer &c) const 
			{ return _field_ptr->cardinality (c); }
    
		/** \brief Characteristic.
		 *
		 * Return c, integer representing characteristic of the field
		 * (the least positive n such that the sum of n copies of x is 0 for all field elements x).
		 * c becomes a positive integer for all fields with finite characteristic,
		 * and 0 to signify a field of infinite characteristic.
		 */
		integer &characteristic (integer &c) const
			{ return _field_ptr->characteristic (c); }
    
		//@} Object Management
    
		/** @name Arithmetic Operations 
		 * x <- y op z; x <- op y
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field elements will
		 * give undefined results.
		 */
		//@{
    
		/** \brief Equality of two elements.
		 *
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and
		 * y, <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return boolean true if equal, false if not.
		 * @param  x field element
		 * @param  y field element
		 */
		bool areEqual (const Element &x, const Element &y) const
			{ return _field_ptr->areEqual (*x._elem_ptr, *y._elem_ptr); }
    
		/** \brief Addition, x <-- y + z.
		 *
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for x, y, and z, 
		 * <tt> _elem_ptr</tt> exists and does not point to null.
		 *
		 * @return reference to x.
		 */
		Element &add (Element &x, const Element &y, const Element &z) const
		{
			_field_ptr->add (*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
			return x;
		}
    
		/** \brief Subtraction, x <-- y - z.
		 *
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for x, y, and z,
		 * <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 */
		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			_field_ptr->sub (*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
			return x;
		}
    
		/** \brief Multiplication, x <-- y * z.
		 *
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for x, y, and z,
		 * <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 */
		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			_field_ptr->mul (*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
			return x;
		}
    
		/** Division, x <-- y / z.
		 *
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for x, y, and z,
		 * <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 */
		Element &div (Element &x, const Element &y, const Element &z) const
		{
			_field_ptr->div (*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
			return x;
		}
    
		/** \brief Additive Inverse (Negation), x <-- - y.
		 *
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 */
		Element &neg (Element &x, const Element &y) const
		{
			_field_ptr->neg (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
    
		/** \brief Multiplicative Inverse, x <-- 1 / y.
		 *
		 * Requires that y is a unit (i.e. nonzero in a field).
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 */
		Element &inv (Element &x, const Element &y) const
		{
			_field_ptr->inv (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
    
    
		/** \brief Field element AXPY, r  <-- a * x + y.
		 *
		 * This function assumes all field elements have already been 
		 * constructed and initialized.
		 * @return reference to r.
		 */
		Element &axpy (Element       &r, 
			       const Element &a,
			       const Element &x, 
			       const Element &y) const
		{
			_field_ptr->axpy (*r._elem_ptr, *a._elem_ptr, *x._elem_ptr,  *y._elem_ptr);
			return r;
		}

		//@} Arithmetic Operations
    
		/** @name Predicates
		*/
		//@{
		/** Zero equality.
		 * Test if field element is equal to zero.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means the <tt>_elem_ptr</tt> 
		 * of x exists and does not point to null.
		 *
		 * @return boolean true if equals zero, false if not.
		 * @param  x field element.
		 */
		bool isZero (const Element &x) const 
			{ return _field_ptr->isZero (*x._elem_ptr); }
    
		/** One equality.
		 * Test if field element is equal to one.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means the <tt> _elem_ptr</tt> 
		 *of x exists and does not point to null.
		 *
		 * @return boolean true if equals one, false if not.
		 * @param  x field element.
		 */
		bool isOne (const Element &x) const 
			{ return _field_ptr->isOne (*x._elem_ptr); }
		//@}

		/** @name Inplace Arithmetic Operations 
		 * x <- x op y; x <- op x
		 * These operations require all elements, including x, to be initialized
		 * before the operation is called.  Uninitialized field elements will
		 * give undefined results.
		 */
		//@{
    
		/** Inplace Addition.
		 * x += y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * <tt> _elem_ptr</tt> exists and does not point to null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &addin (Element &x, const Element &y) const
		{
			_field_ptr->addin (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
    
		/** Inplace Subtraction.
		 * x -= y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &subin (Element &x, const Element &y) const
		{
			_field_ptr->subin (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
 
		/** Inplace Multiplication.
		 * x *= y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &mulin (Element &x, const Element &y) const
		{
			_field_ptr->mulin (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
   
		/** Inplace Division.
		 * x /= y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * <tt> _elem_ptr</tt> exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &divin (Element &x, const Element &y) const
		{
			_field_ptr->divin (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
    
		/** Inplace Additive Inverse (Inplace Negation).
		 * x = - x
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means the <tt>
		 * _elem_ptr</tt> of x exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 */
		Element &negin (Element &x) const
		{
			_field_ptr->negin (*x._elem_ptr);
			return x;
		}
    
		/** Inplace Multiplicative Inverse.
		 * x = 1 / x
		 * This function assumes the field elementhas already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means the <tt>
		 * _elem_ptr</tt> of x exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 */
		Element &invin (Element &x) const
		{
			_field_ptr->invin (*x._elem_ptr);
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
		Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			_field_ptr->axpyin (*r._elem_ptr, *a._elem_ptr, *x._elem_ptr);
			return r;
		}

		//@} Inplace Arithmetic Operations
    
		/** @name Input/Output Operations */
		//@{
    
		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream &write (std::ostream &os) const { return _field_ptr->write (os); }
    
		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream &read (std::istream &is) { return _field_ptr->read (is); }
    
		/** Print field element.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for the <tt>
		 * _elem_ptr</tt> for x exists and does not point to
		 * null.
		 *
		 * @return output stream to which field element is written.
		 * @param  os  output stream to which field element is written.
		 * @param  x   field element.
		 */
		std::ostream &write (std::ostream &os, const Element &x) const 
			{ return _field_ptr->write (os, *x._elem_ptr); }
    
		/** Read field element.
		 * This function assumes the field element has already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for the <tt>
		 * _elem_ptr</tt> for x exists and does not point to
		 * null.
		 *
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		std::istream &read (std::istream &is, Element &x) const
			{ return _field_ptr->read (is, *x._elem_ptr); }
    
		//@} Input/Output Operations
		//@} Common Object Interface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all \ref{LinBox Fields}
		 * and are included only for this implementation of the archetype.
		 */
		//@{

		/** Constructor.
		 * Constructs field from pointer to \ref{FieldAbstract} and its
		 * encapsulated element and random element generator.
		 * Not part of the interface.
		 * Creates new copies of field, element, and random iterator generator
		 * objects in dynamic memory.
		 * @param  field_ptr pointer to \ref{FieldAbstract}.
		 * @param  elem_ptr  pointer to \ref{ElementAbstract}, which is the
		 *                    encapsulated element of \ref{FieldAbstract}.
		 * @param  randIter_ptr  pointer to \ref{RandIterAbstract}, which is the
		 *                        encapsulated random iterator generator
		 *                        of \ref{FieldAbstract}.
		 */
		FieldArchetype (FieldAbstract    *field_ptr,
				 ElementAbstract  *elem_ptr,
				 RandIterAbstract *randIter_ptr = 0)
			: _field_ptr (field_ptr->clone ()), 
			  _elem_ptr (elem_ptr->clone ())
		{
			if (randIter_ptr != 0) _randIter_ptr = randIter_ptr->clone ();
		}

    
		/** Constructor.
		 * Constructs field from ANYTHING matching the interface
		 * using the enveloppe as a \ref{FieldAbstract} and its
		 * encapsulated element and random element generator if needed.
		 * @param  field_ptr pointer to field matching the interface
		 * @param  elem_ptr  pointer to element matching the interface
		 * @param  randIter_ptr  pointer to random matching the interface
		 */
		template<class Field_qcq>
			FieldArchetype (Field_qcq *f) { constructor (f, f); }
	
		//@} Implementation-Specific Methods
    
	    protected:
    
		friend class ElementArchetype;
		friend class RandIterArchetype;
    
		/** Pointer to FieldAbstract object.
		 * Not part of the interface.
		 * Included to allow for archetype use three.
		 */
		mutable FieldAbstract *_field_ptr;
    
		/** Pointer to ElementAbstract object.
		 * Not part of the interface.
		 * Included to allow for archetype use three.
		 */
		mutable ElementAbstract *_elem_ptr;
    
		/** Pointer to RandIterAbstract object.
		 * Not part of the interface.
		 * Included to allow for archetype use three.
		 */
		mutable RandIterAbstract *_randIter_ptr;

		/** Template method for constructing archetype from a derived class of 
		 * FieldAbstract.
		 * This class is needed to help the constructor differentiate between 
		 * classes derived from FieldAbstract and classes that aren't.
		 * Should be called with the same argument to both parameters?
		 * @param	trait	pointer to FieldAbstract or class derived from it
		 * @param	field_ptr	pointer to class derived from FieldAbstract
		 */
		template<class Field_qcq>
		void constructor (FieldAbstract *trait, 
				  Field_qcq      *field_ptr)
		{
			_field_ptr    = field_ptr->clone ();
			_elem_ptr     = static_cast<ElementAbstract*>  (new typename Field_qcq::Element ());
			_randIter_ptr = static_cast<RandIterAbstract*> (new typename Field_qcq::RandIter (*field_ptr));
		}
	 
		/** Template method for constructing archetype from a class not derived 
		 * from FieldAbstract.
		 * This class is needed to help the constructor differentiate between 
		 * classes derived from FieldAbstract and classes that aren't.
		 * Should be called with the same argument to both parameters?
		 * @param	trait	pointer to class not derived from FieldAbstract
		 * @param	field_ptr	pointer to class not derived from FieldAbstract
		 */
		template<class Field_qcq>
		void constructor (void      *trait, 
				  Field_qcq *field_ptr)
		{
			FieldEnvelope< Field_qcq > EnvF (*field_ptr);
			constructor (static_cast<FieldAbstract*> (&EnvF), &EnvF) ;
		}


		/** Only authorize inhertied classes to use the empty constructor
		 **/
		FieldArchetype() {}
		

	}; // class FieldArchetype
  
} // namespace LinBox

#include "linbox/randiter/archetype.h"

#endif // __LINBOX_field_archetype_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
