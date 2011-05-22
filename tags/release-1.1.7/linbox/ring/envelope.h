/* linbox/ring/envelope.h
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
 * ------------------------------------
 */

#ifndef __LINBOX_ring_envelope_H
#define __LINBOX_ring_envelope_H

#include <iostream>

#include "linbox/integer.h"
#include "linbox/element/envelope.h"
#include "linbox/ring/abstract.h"
#include "linbox/element/abstract.h"
#include "linbox/randiter/abstract.h"
#include "linbox/randiter/envelope.h"

#include "linbox/linbox-config.h"
#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

#endif

// Namespace in which all LinBox code resides
namespace LinBox 
{ 
	// Forward declarations
	template <class Ring> class RandIterEnvelope;

	/** 
	 * \brief implement the ring archetype to minimize code bloat.  
	\ingroup ring

	This class implements all purely virtual member functions
	 * of the abstract base class.  This class is used to wrap a
	 * \Ref{LinBox}
	 * ring so that it might be used with the Ring archetype.
	 */
	template <class Ring>
	  class RingEnvelope : public virtual RingAbstract, public virtual FieldEnvelope<Ring>
	{
	    public:

		/** element type.
		 * It is derived from the class ElementAbstract, and it must contain
		 * a wrapped ring element.
		 */
/* 		typedef ElementEnvelope<Ring> Element; */
		typedef typename FieldEnvelope<Ring>::Element Element;
		/** Random iterator generator type.
		 * It is derived from the class RandIterAbstract, and it must contain
		 * a wrapped ring random iterator generator.
		 */
/* 		typedef RandIterEnvelope<Ring> RandIter; */
		typedef typename FieldEnvelope<Ring>::RandIter RandIter;

		/** @name Object Management
		 */
		//@{
 
		/** Default constructor.
		 * In this implementation, this means copying the ring {\tt E.\_field}.
		 */
		RingEnvelope (void) {}

		/** Constructor from ring to be wrapped.
		 * @param F Ring object to be wrapped.
		 */
		RingEnvelope (const Ring& F) : FieldEnvelope<Ring> (F) {}
 
		/** Copy constructor.
		 * Constructs RingEnvelope object by copying the ring.
		 * This is required to allow ring objects to be passed by value
		 * into functions.
		 * In this implementation, this means copying the ring {\tt E.\_field}.
		 * @param  E RingEnvelope object.
		 */
		RingEnvelope (const RingEnvelope& E) : FieldEnvelope<Ring> (E._field) {}

#ifdef __LINBOX_XMLENABLED
		RingEnvelope(Reader &R) : FieldEnvelope<Ring>(R) {}
#endif

 
		/** Virtual copy constructor.
		 * Required because constructors cannot be virtual.
		 * Passes construction on to derived classes.
		 * This function is not part of the common object interface.
		 * @return pointer to new object in dynamic memory.
		 */
// 		RingAbstract* clone () const
// 			{ return new RingEnvelope (*this); }

		FieldAbstract* clone () const
			{ return static_cast<RingAbstract*>( new RingEnvelope (*this)); }

		/** Assignment operator.
		 * Required by abstract base class.
		 * @return reference to RingAbstract object for self
		 * @param F constant reference to RingAbstract object
		 */
		RingAbstract& operator= (const RingAbstract& F)
		{
			if (this != &F) // guard against self-assignment
				_field = static_cast<const RingEnvelope&> (F)._field;

			return *this;
		}

		/** Assignment operator.
		 * Required by abstract base class.
		 * @return reference to RingAbstract object for self
		 * @param F constant reference to RingAbstract object
		 */
		FieldAbstract& operator= (const FieldAbstract& F)
		{
			if (static_cast<RingAbstract*>(this) != &F) // guard against self-assignment
				_field = reinterpret_cast<const RingEnvelope&> (F)._field;

			return * static_cast<RingAbstract*>(this);
		}

		ElementAbstract& init (ElementAbstract& x, const integer& y = 0) const{
			return FieldEnvelope<Ring>::init( x,y);
		}

		integer& convert (integer& x, const ElementAbstract& y) const{
			return FieldEnvelope<Ring>::convert( x,y);
		}

		ElementAbstract& assign (ElementAbstract& x, const ElementAbstract& y) const{
			return FieldEnvelope<Ring>::assign(x,y);
		}

		ElementAbstract& neg (ElementAbstract& x, const ElementAbstract& y) const{
			return FieldEnvelope<Ring>::neg(x,y);
		}
		ElementAbstract& inv (ElementAbstract& x, const ElementAbstract& y) const{
			return FieldEnvelope<Ring>::inv(x,y);
		}
		ElementAbstract& negin (ElementAbstract& x) const{
			return FieldEnvelope<Ring>::negin(x);
		}
		ElementAbstract& invin (ElementAbstract& x) const{
			return FieldEnvelope<Ring>::invin(x);
		}
		integer& cardinality (integer& c) const{
			return FieldEnvelope<Ring>::cardinality(c);
			
		}
		integer& characteristic (integer& c) const{
			return FieldEnvelope<Ring>::characteristic(c);
			
		}
		bool areEqual (const ElementAbstract& x, const ElementAbstract& y) const{
			
			return FieldEnvelope<Ring>::areEqual(x,y);
		}

		ElementAbstract& add (ElementAbstract& x,
				       const ElementAbstract& y,
				      const ElementAbstract& z) const{
			return FieldEnvelope<Ring>::add(x,y,z);
		}
		ElementAbstract& sub (ElementAbstract& x,
				       const ElementAbstract& y,
				      const ElementAbstract& z) const{
			return FieldEnvelope<Ring>::sub(x,y,z);
		}
		ElementAbstract& mul (ElementAbstract& x,
				       const ElementAbstract& y,
				      const ElementAbstract& z) const{
			return FieldEnvelope<Ring>::mul(x,y,z);
		}

		ElementAbstract& div (ElementAbstract& x,
				       const ElementAbstract& y,
				      const ElementAbstract& z) const{
			return FieldEnvelope<Ring>::div(x,y,z);
		}
		
		ElementAbstract& axpy (ElementAbstract& r, 
					const ElementAbstract& a, 
					const ElementAbstract& x, 
					const ElementAbstract& y) const
		{
			return FieldEnvelope<Ring>::axpy(r,a,x,y);
		}
		ElementAbstract& addin (ElementAbstract& x,
					const ElementAbstract& z) const{
			return FieldEnvelope<Ring>::addin(x,z);
		}
		ElementAbstract& subin (ElementAbstract& x,
				      const ElementAbstract& z) const{
			return FieldEnvelope<Ring>::subin(x,z);
		}
		ElementAbstract& mulin (ElementAbstract& x,
				      const ElementAbstract& z) const{
			return FieldEnvelope<Ring>::mulin(x,z);
		}

		ElementAbstract& divin (ElementAbstract& x,
					const ElementAbstract& z) const{
			return FieldEnvelope<Ring>::divin(x,z);
		}
		
		ElementAbstract& axpyin(ElementAbstract& r, 
					const ElementAbstract& x, 
					const ElementAbstract& y) const
		{
			return FieldEnvelope<Ring>::axpyin(r,x,y);
		}
		bool isZero (const ElementAbstract& x) const{
			return FieldEnvelope<Ring>::isZero(x);
		}

		
		bool isOne (const ElementAbstract& x) const{
			return FieldEnvelope<Ring>::isOne(x);
		}

		/** Invertibility test.
		 * Test if ring element is invertible.
		 * This function assumes the ring element has already been
		 * constructed and initialized.
		 *
		 * @return boolean true if equals zero, false if not.
		 * @param  x ring element.
		 */

		bool isUnit (const ElementAbstract& x) const
			{ return _field.isUnit (static_cast<const ElementEnvelope<Ring>&> (x)._elem); }
 
		/** Divisibility of zero test.
		 * Test if ring element is a zero divisor.
		 * This function assumes the ring element has already been
		 * constructed and initialized.
		 *
		 * @return boolean true if divides zero, false if not.
		 * @param  x ring element.
		 */
		bool isZeroDivisor (const ElementAbstract& x) const
			{ return _field.isZeroDivisor (static_cast<const ElementEnvelope<Ring>&> (x)._elem); }


		std::ostream& write (std::ostream& os) const 
		{ return  FieldEnvelope<Ring>::write (os); }
		
		std::istream& read (std::istream& is) { return FieldEnvelope<Ring>::read (is); }

		std::ostream& write (std::ostream& os, const ElementAbstract& x) const
			{ return FieldEnvelope<Ring>::write (os, x); }
 
		std::istream& read (std::istream& is, ElementAbstract& x) const
			{ return FieldEnvelope<Ring>::read (is, x); }
	    private:

		friend class RandIterEnvelope<Ring>;


	}; // class RingEnvelope

} // namespace LinBox

#include "linbox/randiter/envelope.h"

#endif // __LINBOX_ring_envelope_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
