/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/archetype.h
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
 *
 * ------------------------------------
 * 2002-05-14 William J. Turner <wjturner@acm.org>
 * 
 * changed randIter to RandIter.
 * ------------------------------------
 */


#ifndef __FIELD_ARCHETYPE_H
#define __FIELD_ARCHETYPE_H

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
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <iostream>
#include <string>
using std::string;
using std::ostream;

#endif

#include "linbox/util/error.h"

namespace LinBox
{
	// Forward declarations
	class RandIterArchetype;
/*	
#ifdef __LINBOX_XMLENABLED
	class FieldArchetypeFTor;
#endif
*/

	/** Field Archetype.
	 * Archetype for the field common object interface to \Ref{LinBox}.
	 *
	 * The \Ref{FieldArchetype} and its encapsulated
	 * element class contain pointers to the \Ref{FieldAbstract}
	 * and its encapsulated field element, respectively.
	 * \Ref{FieldAbstract} then uses virtual member functions to
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
		 * These methods are required of all \Ref{LinBox} fields.
		 */
		//@{
    
		/// element type.
		typedef ElementArchetype Element;

		/// Random iterator generator type.
		typedef RandIterArchetype RandIter;
    
		/// @name Object Management
		//@{
    
		/** Copy constructor.
		 *
		 * Constructs FieldArchetype object by copying the
		 * field.  This is required to allow field objects to
		 * be passed by value into functions.
		 *
		 * In this implementation, this means copying the
		 * field to which {\tt F.\_field\_ptr} points, the
		 * element to which {\tt F.\_elem\_ptr} points, and the
		 * random element generator to which
		 * {\tt F.\_randIter\_ptr} points.
		 *
		 * @param F {\tt FieldArchetype} object.
		 */
		FieldArchetype (const FieldArchetype &F) 
		{ 
			if (F._field_ptr != 0) _field_ptr = F._field_ptr->clone (); 
			if (F._elem_ptr != 0) _elem_ptr = F._elem_ptr->clone ();
			if (F._randIter_ptr != 0) _randIter_ptr = F._randIter_ptr->clone ();
		}
/*		
#ifdef __LINBOX_XMLENABLED

		friend class FieldArchetypeFtor;
		FieldArchetype(Reader &R);
#endif
*/		


		/** Destructor.
		 *
		 * This destroys the field object, but it does not
		 * destroy any field element objects.  The destructor
		 * for each field element must also be called.
		 *
		 * In this implementation, destruction is deletion of
		 * the field object to which {\tt \_field\_ptr}
		 * points, the field element to which {\tt
		 * \_elem\_ptr} points, and the random element
		 * generator to which {\tt \_randIter\_ptr} points.
		 */
		~FieldArchetype (void) 
		{
			if (_field_ptr != 0) delete _field_ptr;
			if (_elem_ptr != 0) delete _elem_ptr; 
			if (_randIter_ptr != 0) delete _randIter_ptr;
		}
    
		/** Assignment operator.
		 *
		 * Assigns {\tt FieldArchetype} object F to field.  In
		 * this implementation, this means copying the field
		 * to which {\tt F.\_field\_ptr} points, the element to which
		 * {\tt F.\_elem\_ptr} points, and the random element
		 * generator to which {\tt F.\_randIter\_ptr} points.
		 *
		 * @param F {\tt FieldArchetype} object.
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
    
		/** Initialization of field element from an integer.
		 * Behaves like C++ allocator construct.
		 *
		 * This function assumes the output field element x
		 * has already been constructed, but that it is not
		 * necessarily already initialized. In this
		 * implementation, this means the {\tt \_elem\_ptr} of
		 * x exists, but that it may be the null pointer.
		 *
		 * @return reference to field element.
		 * @param x field element to contain output (reference returned).
		 * @param y constant reference to integer.
		 */
		Element &init (Element &x, const integer &y = 0 ) const
		{
			if (x._elem_ptr != 0) delete x._elem_ptr;
			x._elem_ptr = _elem_ptr->clone ();
			_field_ptr->init (*x._elem_ptr, y);
			return x;
		}
  
		/** Conversion of field element to an integer.
		 *
		 * This function assumes the output field element x
		 * has already been constructed, but that it is not
		 * already initialized. In this implementation, this
		 * means the {\tt \_elem\_ptr} of y exists, and that
		 * it is not the null pointer.
		 *
		 * @return reference to integer.
		 * @param x reference to integer to contain output (reference returned).
		 * @param y constant reference to field element.
		 */
		integer &convert (integer &x, const Element &y = 0) const
		{
			_field_ptr->convert (x, *y._elem_ptr);
			return x;
		}
    
		/** Assignment of one field element to another.
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and
		 * y, {\tt \_elem\_ptr} exists and does not point to
		 * null.
		 *
		 * @return reference to x
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &assign (Element &x, const Element &y) const
		{
			if (x._elem_ptr == 0) 
				x._elem_ptr = _elem_ptr->clone ();

			_field_ptr->assign (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
    
		/** Cardinality.
		 * Return integer representing cardinality of the field.
		 * Returns a non-negative integer for all fields with finite
		 * cardinality, and returns -1 to signify a field of infinite 
		 * cardinality.
		 * @return constant reference to integer representing cardinality 
		 *	       of the field
		 */
		integer &cardinality (integer &c) const 
			{ return _field_ptr->cardinality (c); }
    
		/** Characteristic.
		 * Return integer representing characteristic of the field.
		 * Returns a positive integer to all fields with finite characteristic,
		 * and returns 0 to signify a field of infinite characteristic.
		 * @return constant reference to integer representing characteristic 
		 * 	       of the field.
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
    
		/** Equality of two elements.
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and
		 * y, {\tt \_elem\_ptr} exists and does not point to
		 * null.
		 *
		 * @return boolean true if equal, false if not.
		 * @param  x field element
		 * @param  y field element
		 */
		bool areEqual (const Element &x, const Element &y) const
			{ return _field_ptr->areEqual (*x._elem_ptr, *y._elem_ptr); }
    
		/** Addition.
		 * x = y + z
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for x, y, and z, 
		 * {\tt \_elem\_ptr} exists and does not point to null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 * @param  z field element.
		 */
		Element &add (Element &x, const Element &y, const Element &z) const
		{
			_field_ptr->add (*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
			return x;
		}
    
		/** Subtraction.
		 * x = y - z
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for x, y, and z,
		 * {\tt \_elem\_ptr} exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 * @param  z field element.
		 */
		Element &sub (Element &x, const Element &y, const Element &z) const
		{
			_field_ptr->sub (*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
			return x;
		}
    
		/** Multiplication.
		 * x = y * z
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for x, y, and z,
		 * {\tt \_elem\_ptr} exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 * @param  z field element.
		 */
		Element &mul (Element &x, const Element &y, const Element &z) const
		{
			_field_ptr->mul (*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
			return x;
		}
    
		/** Division.
		 * x = y / z
		 * This function assumes all the field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for x, y, and z,
		 * {\tt \_elem\_ptr} exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 * @param  z field element.
		 */
		Element &div (Element &x, const Element &y, const Element &z) const
		{
			_field_ptr->div (*x._elem_ptr, *y._elem_ptr, *z._elem_ptr);
			return x;
		}
    
		/** Additive Inverse (Negation).
		 * x = - y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * {\tt \_elem\_ptr} exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &neg (Element &x, const Element &y) const
		{
			_field_ptr->neg (*x._elem_ptr, *y._elem_ptr);
			return x;
		}
    
		/** Multiplicative Inverse.
		 * x = 1 / y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * {\tt \_elem\_ptr} exists and does not point to
		 * null.
		 *
		 * @return reference to x.
		 * @param  x field element (reference returned).
		 * @param  y field element.
		 */
		Element &inv (Element &x, const Element &y) const
		{
			_field_ptr->inv (*x._elem_ptr, *y._elem_ptr);
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
		Element &axpy (Element       &r, 
			       const Element &a,
			       const Element &x, 
			       const Element &y) const
		{
			_field_ptr->axpy (*r._elem_ptr, *a._elem_ptr, *x._elem_ptr,  *y._elem_ptr);
			return r;
		}

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
		 *
		 * In this implementation, this means the {\tt
		 * \_elem\_ptr} of x exists and does not point to
		 * null.
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
		 * In this implementation, this means the {\tt
		 * \_elem\_ptr} of x exists and does not point to
		 * null.
		 *
		 * @return boolean true if equals one, false if not.
		 * @param  x field element.
		 */
		bool isOne (const Element &x) const 
			{ return _field_ptr->isOne (*x._elem_ptr); }
    
		/** Inplace Addition.
		 * x += y
		 * This function assumes both field elements have already been 
		 * constructed and initialized.
		 *
		 * In this implementation, this means for both x and y
		 * {\tt \_elem\_ptr} exists and does not point to
		 * null.
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
		 * {\tt \_elem\_ptr} exists and does not point to
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
		 * {\tt \_elem\_ptr} exists and does not point to
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
		 * {\tt \_elem\_ptr} exists and does not point to
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
		 * In this implementation, this means the {\tt
		 * \_elem\_ptr} of x exists and does not point to
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
		 * In this implementation, this means the {\tt
		 * \_elem\_ptr} of x exists and does not point to
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

#ifndef __LINBOX_XMLENABLED 
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
		 * In this implementation, this means for the {\tt
		 * \_elem\_ptr} for x exists and does not point to
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
		 * In this implementation, this means for the {\tt
		 * \_elem\_ptr} for x exists and does not point to
		 * null.
		 *
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		std::istream &read (std::istream &is, Element &x) const
			{ return _field_ptr->read (is, *x._elem_ptr); }
    
		//@} Input/Output Operations
#else
		ostream &write(ostream &os) const
		{
			if(_field_ptr)
				return _field_ptr->write(os);
			else
				return os;
		}

	        bool toTag(Writer &W) const
		{
			if(_field_ptr)
				return _field_ptr->toTag(W);
			else 
				return false;
			
		}

		ostream &write(ostream &os, const Element &e) const
		{
			if(_field_ptr)
				return _field_ptr->write(os, *(e._elem_ptr));
			else
				return os;
		}

		bool toTag(Writer &W, const Element &e) const
		{
			if(_field_ptr)
				return _field_ptr->toTag(W, *(e._elem_ptr));
			else 
				return false;
			
		}

		istream &read(istream &is, Element &e) const
		{
			if(_field_ptr)
				return _field_ptr->read(is, *(e._elem_ptr));
			else
				return is;
		}

		bool fromTag(Reader &R, Element &e) const
		{
			if(_field_ptr)
				return _field_ptr->fromTag(R, *(e._elem_ptr));
			else
				return false;
		}
#endif    
		//@} Common Object Interface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all \Ref{LinBox Fields}
		 * and are included only for this implementation of the archetype.
		 */
		//@{

		/** Constructor.
		 * Constructs field from pointer to \Ref{FieldAbstract} and its
		 * encapsulated element and random element generator.
		 * Not part of the interface.
		 * Creates new copies of field, element, and random iterator generator
		 * objects in dynamic memory.
		 * @param  field\_ptr pointer to \Ref{FieldAbstract}.
		 * @param  elem\_ptr  pointer to \Ref{ElementAbstract}, which is the
		 *                    encapsulated element of \Ref{FieldAbstract}.
		 * @param  randIter\_ptr  pointer to \Ref{RandIterAbstract}, which is the
		 *                        encapsulated random iterator generator
		 *                        of \Ref{FieldAbstract}.
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
		 * using the enveloppe as a \Ref{FieldAbstract} and its
		 * encapsulated element and random element generator if needed.
		 * @param  field\_ptr pointer to field matching the interface
		 * @param  elem\_ptr  pointer to element matching the interface
		 * @param  randIter\_ptr  pointer to random matching the interface
		 */
		template<class Field_qcq>
			FieldArchetype (Field_qcq *f) { constructor (f, f); }
	
		//@} Implementation-Specific Methods
    
	    private:
    
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
		 * @param	field\_ptr	pointer to class derived from FieldAbstract
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
		 * @param	field\_ptr	pointer to class not derived from FieldAbstract
		 */
		template<class Field_qcq>
		void constructor (void      *trait, 
				  Field_qcq *field_ptr)
		{
			FieldEnvelope< Field_qcq > EnvF (*field_ptr);
			constructor (static_cast<FieldAbstract*> (&EnvF), &EnvF) ;
		}

	}; // class FieldArchetype
  
} // namespace LinBox
/*
#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/field-reader-analyzer.h"

namespace LinBox {


	class FieldArchetypeFtor {
		public:
			FieldArchetypeFtor(FieldArchetype &A) : _A(A) {}
			
			template<class Field>
			void* operator()(Field* F) 
			{
				_A.constructor(F, F);
				
				//				FieldEnvelope<Field> Envelope(*F);
				delete F;

				//				_A._field_ptr = new FieldEnvelope<Field>(Envelope);
				return NULL;
			}
			
			FieldArchetype &_A;
	};


		// Note - This constructor is potentially dangerous as it 
		// makes a call to FieldReaderAnalyzer's makeField method, 
		// which can throw errors.  This class does not catch the
		// errors, as they are not the responsibility of FieldArchetype
		// but of the user code
		//
	FieldArchetype::FieldArchetype(Reader &R)
	{
		FieldReaderAnalyzer Analyzer(R);
		FieldArchetypeFtor Functor(*this);
		
		Analyzer.makeField(Functor);
		
		return;
	}
}
#endif
*/


#include "linbox/randiter/archetype.h"

#endif // __FIELD_ARCHETYPE_H
