/* linbox/ring/abstract.h
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
 */

#ifndef __LINBOX_ring_abstract_H
#define __LINBOX_ring_abstract_H

#include <iostream>

#include "linbox/element/abstract.h"
#include "linbox/randiter/abstract.h"
#include "linbox/integer.h"
#include "linbox/linbox-config.h"
#include "linbox/field/abstract.h"
#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#endif

namespace LinBox 
{ 

	/** \brief Abstract ring base class.
	\ingroup ring 

	 * Found in the file \URL{linbox/ring/abstract.h}.
	 * Abstract base class used to implement the ring archetype to minimize
	 * code bloat.  All public member functions of this class are purely
	 * virtual and must be implemented by all derived classes.
	 *
	 * If a template is instantiated on the ring archetype, we can change the
	 * ring it is using by changing the derived class of this class.  This allows
	 * us to change the ring used in a template without having to reinstantiate
	 * it.  This minimizes code bloat, but it also introduces indirection through
	 * the use of pointers and virtual functions which is inefficient.
	 */
  class RingAbstract : public FieldAbstract
	{
	    public:

/* 		/// element type. */
/* 		typedef ElementAbstract Element; */

/* 		/// Random iterator generator type. */
/* 		typedef RandIterAbstract RandIter; */
		typedef FieldAbstract::Element Element;
		typedef FieldAbstract::RandIter RandIter;
		/** @name Object Management
		 * There are no public constructors for this class.
		 * It should only be used in tandem with \Ref{FieldArchetype}.
		 */
		//@{

		/** Destructor.
		 * Required because of virtual member functions.
		 * Virtual.
		 */
	  virtual ~RingAbstract (void) {}

	  /** Virtual copy constructor.
	   * Required because constructors cannot be virtual.
	   * Passes construction on to derived classes.
	   * Purely virtual.
	   * This function is not part of the common object interface.
	   * @return pointer to new object in dynamic memory.
	   */
	  //virtual RingAbstract* clone () const = 0;

		/** Invertibility test.
		 * Test if ring element is invertible.
		 * This function assumes the ring element has already been
		 * constructed and initialized.
		 * Purely virtual.
		 * @return boolean true if invertible, false if not.
		 * @param  x ring element.
		 */
		virtual bool isUnit (const Element &x) const = 0;
 
		/** Divisibility of zero test.
		 * Test if ring element is a zero divisor.
		 * This function assumes the ring element has already been
		 * constructed and initialized.
		 * Purely virtual.
		 * @return boolean true if divides zero, false if not.
		 * @param  x ring element.
		 */
		virtual bool isZeroDivisor (const Element &x) const = 0;


	    private:

		/// FieldArchetype is friend.
		friend class RingArchetype;

	}; // class FieldAbstract

}  // namespace LinBox

#endif // __LINBOX_ring_abstract_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
