/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/unparametric.h
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
 */
 
#ifndef __FIELD_UNPARAMETRIC_H
#define __FIELD_UNPARAMETRIC_H

#include <string>
#include <algorithm>

#include "linbox/integer.h"
#include <linbox/field/field-interface.h>
#include "linbox/randiter/unparametric.h"
#include "linbox-config.h"
#include <linbox/field/field-traits.h>

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#include <iostream>

#endif

namespace LinBox 
{
	/** \brief Unparameterized field adapter.
	\ingroup field

	A field having an interface similar to that of floats is adapted to LinBox.

	Used to generate efficient field classes for unparameterized fields (or hidden parameter fields).
	 
	Some fields are implemented by definition of the C++ arithmetic operators, such as z = x*y, 
	for z, y, z instances of a type K.   The LinBox field 
	Unparametric<K> is the adaptation to LinBox.

	For a typical unparametric field, some of the methods must be defined in a specialization. 

	*/

	template <class Ring>
	struct ClassifyRing;

	template <class K>
	class UnparametricField;

	template <class K>
	struct ClassifyRing<UnparametricField<K> > { 
		typedef RingCategories::GenericTag categoryTag;
	};

	template <class K>
	class UnparametricField : public FieldInterface
	{
	    protected: integer _p; integer _card;
	    public:
    
		/** @name Common Object Interface for a LinBox Field.
		 * These methods and member types are required of all LinBox fields.
		 * See \ref{FieldArchetype} for detailed specifications.
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
		UnparametricField(integer q = 0, size_t e = 1)
		: _p(q), _card(q == 0 ? integer(-1) : pow(q, e) ) {}  // assuming q is a prime or zero.

		/// construct this field as copy of F.
		UnparametricField (const UnparametricField &F) : _p(F._p), _card(F._card){}
    
		/** XML constructor
		 * Takes in an XML reader and attempts to convert the
		 * XML representation over to a valid field.  As this class
		 * is mostly a wrapper for a particular field type, the
		 * XML does little more than encode the cardinality of this
		 * field, and perhaps the characteristic
		 *
		 */
#ifdef __LINBOX_XMLENABLED

		// The default XML constructor for the unparametric field
		// Doesn't really do anything to set the characteristic and
		// cardinality of the field.  This constructor will not carry
		// over to template specializations of this field, so
		// this constructor must be overloaded
		//
		UnparametricField(LinBox::Reader &R)
		{
			if(!R.expectTagName("field")) return;
			if(!R.expectAttributeNum("cardinality", _card)) return;
			if(!R.expectChildTag()) return;
			R.traverseChild();

			// infinite fields
			if(R.checkTagName("integer") || R.checkTagName("rational") || R.checkTagName("real")) {
				_p =   Integer(0);
				R.upToParent();
			}
			// finite fields
			else if(R.checkTagName("finite")) {
				if(!R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagName("characteristic") || !R.expectChildTextNum(_p)) return;
				R.upToParent();
				R.upToParent();
			}
			else if(!R.expectTagName("unknown"))
				return;
			else {
				if(!R.expectChildTag()) return;
				R.traverseChild();
				if(!R.expectTagName("characteristic") || !R.expectChildTextNum(_p)) return;
				return;
			}
		}
#endif

		/// 
		~UnparametricField () {}
    
		/* Assignment operator.
		 * Assigns UnparametricField object F to field.
		 * @param  F UnparametricField object.
		 */
		 // I believe this should be virtual -bds
		///
		const UnparametricField &operator=(const UnparametricField &F) const { return *this; }
		//@} Field Object Basics.
    
		/** @name Data Object Management.
		 * first argument is set and the value is also returned.
		 */
	        //@{

		/// x := y.  Caution: it is via cast to long.  Good candidate for specialization.
		Element &init (Element &x, const integer &y=0) const 
			{ return x = static_cast<const Element&> (static_cast<const long&> (y)); }

    
		/// x :=  y.  Caution: it is via cast to long.  Good candidate for specialization.
		integer &convert (integer &x, const Element &y) const 
		{ 
			Element temp (y);
			//Dan Roche changed this from long to integer.
			return x = static_cast<integer> (temp); 
		}
    
		/// x :=  y.  Caution: it is via cast to long.  Good candidate for specialization. --dpritcha
		double &convert (double &x, const Element &y) const 
		{ 
			Element temp (y);
			return x = static_cast<double> (temp); 
		}
    
		///
		Element &assign (Element &x, const Element &y) const { return x = y; }
    
		/// c := cardinality of this field (-1 if infinite).
		integer &cardinality (integer &c) const { return c = _card; }
    
		/// c := characteristic of this field (zero or prime).
		integer &characteristic (integer &c) const { return c = _p; }
    
		//@} Data Object Management
    
		/// @name Comparison Predicates
		//@{
		///  x == y
		bool areEqual (const Element &x, const Element &y) const { return x == y; }
    
		///  x == 0
		bool isZero (const Element &x) const { return x == Element (0); }
    
		///  x == 1
		bool isOne (const Element &x) const { return x == Element (1); }
		//@} Comparison Predicates
    
		 
		/** @name Arithmetic Operations 
		 * The first argument is set and is also the return value.
		 */
		//@{
    
		/// x := y + z
		Element &add (Element &x, const Element &y, const Element &z) const
			{ return x = y + z; }
    
		/// x := y - z
		Element &sub (Element &x, const Element &y, const Element &z) const
			{ return x = y - z; }
    
		/// x := y*z
		Element &mul (Element &x, const Element &y, const Element &z) const
			{ return x = y * z; }
    
		/// x := y/z
		Element &div (Element &x, const Element &y, const Element &z) const
			{ return x = y / z; }
    
		/// x := -y
		Element &neg (Element &x, const Element &y) const { return x = - y; }
    
		/// x := 1/y
		Element &inv (Element &x, const Element &y) const 
			{ return x = Element (1) / y; }
    
		/// z := a*x + y 
		// more optimal implementation, if available, can be defined in a template specialization.
		Element &axpy (Element &z, 
			       const Element &a, 
			       const Element &x, 
			       const Element &y) const
			{ return z = a * x + y; }
 
		//@} Arithmetic Operations
    
		/** @name Inplace Arithmetic Operations 
		 * The first argument is modified and the result is the return value.
		 */
		//@{
    
		/// x := x + y
		Element &addin (Element &x, const Element &y) const { return x += y; }
    
		/// x := x - y
		Element &subin (Element &x, const Element &y) const { return x -= y; }
    
		/// x := x*y
		Element &mulin (Element &x, const Element &y) const { return x *= y; }
    
		/// x := x/y
		Element &divin (Element &x, const Element &y) const { return x /= y; }
    
		/// x := -x
		Element &negin (Element &x) const { return x = - x; }
    
		/// x := 1/x
		Element &invin (Element &x) const { return x = Element (1) / x; }
    
		/// y := a*x + y
		Element &axpyin (Element &y, const Element &a, const Element &x) const
			{ return y += a * x; }
 
		//@} Inplace Arithmetic Operations

#ifndef __LINBOX_XMLENABLED
    
		/** @name Input/Output Operations */
		//@{
    
		/** Print field.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 */
		std::ostream &write (std::ostream &os) const { return os << "unparamterized field"; }
    
		/** Read field.
		 * @return input stream from which field is read.
		 * @param  is  input stream from which field is read.
		 */
		std::istream &read (std::istream &is) const { return is; }
    
		/** Print field element.
		 * @return output stream to which field element is written.
		 * @param  os  output stream to which field element is written.
		 * @param  x   field element.
		 */
		std::ostream &write (std::ostream &os, const Element &x) const { return os << x; }
    
		/** Read field element.
		 * @return input stream from which field element is read.
		 * @param  is  input stream from which field element is read.
		 * @param  x   field element.
		 */
		std::istream &read (std::istream &is, Element &x) const { return is >> x; }
    
		//@}
#else // These are the XML writing methods.  They presume nothing about what 
      // the field is, but merely write stuff out.  Note, for the elements,
      // the writer relies on the operator<< function of the ostringstream
      // class to produce a proper std::string for use with the XML

		
		std::ostream &write(std::ostream &out) const 
		{
			LinBox::Writer W;
			if( toTag(W))
				W.write(out);

			return out;
		}

		std::ostream &write(std::ostream &out, const Element &e) const 
		{
			LinBox::Writer W;
			if( toTag(W, e))
				W.write(out);

			return out;
              }
			
		// for this function, the type of the field is unknown, so
		// just print an unknown field
		bool toTag(LinBox::Writer &W) const
		{
			std::string s;
			W.setTagName("field");
			W.setAttribute("cardinality", LinBox::Writer::numToString(s, _card));
			W.setAttribute("implDetail", "unknown");

			W.addTagChild();
			W.setTagName("unknown");

			W.addTagChild();
			W.setTagName("characteristic");
			W.addDataChild(LinBox::Writer::numToString(s, _p));
			W.upToParent();
			
			W.upToParent();

			return true;
		}

		bool toTag(LinBox::Writer &W, const Element &e) const
		{
			std::string s;

			W.setTagName("cn");

			// use the ostrinstream operator<<, let the C++
			// constructor sort it out.  If Element is a built-in 
			// type,this works.  If it isn't, you should have
			// made a template specialization.  Sorry :-)

			W.addDataChild(LinBox::Writer::numToString(s, e));

			return true;
		}

		std::istream &read(std::istream &in, Element &e) const
		{
			LinBox::Reader R(in);
			if( !fromTag(R, e)) {
				in.setstate(std::istream::failbit);
				if(!R.initalized()) {
					in.setstate(std::istream::badbit);
				}
			}

			return in;
		}

		// this is the "shot-in-the-dark" method of element
		// initalization.  Essentially, read in the element tag and
		// attempt to write 
		bool fromTag(LinBox::Reader &R, Element &e) const
		{
			std::string s;
			if(!R.expectTagName("cn") || !R.expectChildTextNum(s)) return false;

			return true;
		}

			
#endif

    
		//@} Common Object Interface
    
		/** @name Implementation-Specific Methods.
		 * These methods are not required of all LinBox fields
		 * and are included only for the implementation of this field
		 * template.
		 */
		//@{
    
		/// Default constructor
		//UnparametricField (void) {}
    
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

	template<class Field>
	class FieldAXPY;

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
                FieldAXPY (const Field &F) : _F (F) { _y = 0; }
 
                /** Copy constructor.
                 * @param faxpy
                 */
                FieldAXPY (const FieldAXPY<Field> &faxpy) : _F (faxpy._F), _y (faxpy._y) {}
 
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
                 * @param y_init constant reference to element a
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

#endif // __FIELD_UNPARAMETRIC_H_
