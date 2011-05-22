/* linbox/ring/archetype.h
 * Copyright(C) LinBox
 * Written by J-G Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Clement Pernet <Clement.Pernet@imag.fr>
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
 */


#ifndef __LINBOX_ring_archetype_H
#define __LINBOX_ring_archetype_H

#include <iostream>
#include "linbox/field/archetype.h"
#include "linbox/ring/ring-interface.h"
#include "linbox/ring/abstract.h"
#include "linbox/ring/envelope.h"
#include "linbox/element/archetype.h"
#include "linbox/element/abstract.h"
#include "linbox/element/envelope.h"
#include "linbox/randiter/abstract.h"
#include "linbox/randiter/envelope.h"
#include "linbox/randiter/archetype.h"
#include "linbox/integer.h"
#include "linbox/linbox-config.h"

#include "linbox/util/error.h"

namespace LinBox
{ 
	// Forward declarations
	class RandIterArchetype;

	/** 
	 * \brief specification and archetypic instance for the ring interface 
	\ingroup ring
	 *
	 * The \Ref{RingArchetype} and its encapsulated
	 * element class contain pointers to the \Ref{RingAbstract}
	 * and its encapsulated ring element, respectively.
	 * \Ref{RingAbstract} then uses virtual member functions to
	 * define operations on its encapsulated ring element.  This ring 
	 * element has no knowledge of the ring properties being used on it 
	 * which means the ring object must supply these operations.
	 *
	 * It does not contain elements zero and one because they can be created 
	 * whenever necessary, although it might be beneficial from an efficiency
	 * stand point to include them.  However, because of archetype use three,
	 * the elements themselves cannot be contained, but rather pointers to them.
	 */
	class RingArchetype : public virtual FieldArchetype
	{
	    public:

		/** @name Common Object Interface for a LinBox Ring.
		 * These methods are required of all \Ref{LinBox} rings.
		 */
		//@{
    
		/// element type.
/* 		typedef ElementArchetype Element; */
		typedef FieldArchetype::Element Element;
		/// Random iterator generator type.
/* 		typedef RandIterArchetype RandIter; */
		typedef FieldArchetype::RandIter RandIter;
		/// @name Object Management
		//@{
    
		/** Copy constructor.
		 *
		 * Constructs RingArchetype object by copying the
		 * ring.  This is required to allow ring objects to
		 * be passed by value into functions.
		 *
		 * In this implementation, this means copying the
		 * ring to which {\tt F.\_ring\_ptr} points, the
		 * element to which {\tt F.\_elem\_ptr} points, and the
		 * random element generator to which
		 * {\tt F.\_randIter\_ptr} points.
		 *
		 * @param F {\tt RingArchetype} object.
		 */
	  RingArchetype (const RingArchetype &F) : FieldArchetype ( F )	{ }
    
	  /** \brief Invertibility test.
		 * Test if ring element is invertible.
		 * This function assumes the ring element has already been
		 * constructed and initialized.
		 * In this implementation, this means the {\tt
		 * \_elem\_ptr} of x exists and does not point to
		 * null.
		 *
		 * @return boolean true if equals zero, false if not.
		 * @param  x ring element.
		 */
	  bool isUnit (const Element &x) const 
	    { return _ring_ptr->isUnit (*x._elem_ptr); }
    
		/** Divisibility of zero test.
		 * Test if ring element is a zero divisor.
		 * This function assumes the ring element has already been
		 * constructed and initialized.
		 *
		 * In this implementation, this means the {\tt
		 * \_elem\_ptr} of x exists and does not point to
		 * null.
		 *
		 * @return boolean true if divides zero, false if not.
		 * @param  x ring element.
		 */
		bool isZeroDivisor (const Element &x) const 
			{ return _ring_ptr->isZeroDivisor (*x._elem_ptr); }
    

		/** Constructor.
		 * Constructs ring from pointer to \Ref{RingAbstract} and its
		 * encapsulated element and random element generator.
		 * Not part of the interface.
		 * Creates new copies of ring, element, and random iterator generator
		 * objects in dynamic memory.
		 * @param  ring\_ptr pointer to \Ref{RingAbstract}.
		 * @param  elem\_ptr  pointer to \Ref{ElementAbstract}, which is the
		 *                    encapsulated element of \Ref{RingAbstract}.
		 * @param  randIter\_ptr  pointer to \Ref{RandIterAbstract}, which is the
		 *                        encapsulated random iterator generator
		 *                        of \Ref{RingAbstract}.
		 */
		RingArchetype (RingAbstract    *ring_ptr,
				 ElementAbstract  *elem_ptr,
				 RandIterAbstract *randIter_ptr = 0)
		  : FieldArchetype( static_cast<FieldAbstract*>(ring_ptr->clone()),
				    elem_ptr, randIter_ptr ),
		  _ring_ptr (dynamic_cast<RingAbstract*>(ring_ptr->clone ()))
		{
		}

    
		/** Constructor.
		 * Constructs ring from ANYTHING matching the interface
		 * using the enveloppe as a \Ref{RingAbstract} and its
		 * encapsulated element and random element generator if needed.
		 * @param  ring\_ptr pointer to ring matching the interface
		 * @param  elem\_ptr  pointer to element matching the interface
		 * @param  randIter\_ptr  pointer to random matching the interface
		 */
		template<class Ring_qcq>
		RingArchetype (Ring_qcq *f)
		{ Ring_constructor (f, f); }
	
		//@} Implementation-Specific Methods
    
	    private:
    
		friend class ElementArchetype;
		friend class RandIterArchetype;
    
		/** Pointer to RingAbstract object.
		 * Not part of the interface.
		 * Included to allow for archetype use three.
		 */
		mutable RingAbstract *_ring_ptr;
    

		/** Template method for constructing archetype from a derived class of 
		 * RingAbstract.
		 * This class is needed to help the constructor differentiate between 
		 * classes derived from RingAbstract and classes that aren't.
		 * Should be called with the same argument to both parameters?
		 * @param	trait	pointer to RingAbstract or class derived from it
		 * @param	ring\_ptr	pointer to class derived from RingAbstract
		 */
		template<class Ring_qcq>
		void Ring_constructor (RingAbstract *trait, 
				  Ring_qcq     *ring_ptr)
		{
			constructor( static_cast<FieldAbstract*>(trait), ring_ptr);
			_ring_ptr    = dynamic_cast<RingAbstract*>(ring_ptr->clone ());

		}
	 
		/** Template method for constructing archetype from a class not derived 
		 * from RingAbstract.
		 * This class is needed to help the constructor differentiate between 
		 * classes derived from RingAbstract and classes that aren't.
		 * Should be called with the same argument to both parameters?
		 * @param	trait	pointer to class not derived from RingAbstract
		 * @param	ring\_ptr	pointer to class not derived from RingAbstract
		 */
		template<class Ring_qcq>
		void Ring_constructor (void      *trait, 
				  Ring_qcq *ring_ptr)
		{
			RingEnvelope< Ring_qcq > EnvF (*ring_ptr);
			Ring_constructor (static_cast<RingAbstract*> (&EnvF), &EnvF) ;
		}

	}; // class RingArchetype
  
}  // namespace LinBox


#endif // __LINBOX_ring_archetype_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
