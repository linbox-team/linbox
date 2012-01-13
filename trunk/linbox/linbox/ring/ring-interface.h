/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/ring/ring-interface.h
 * Copyright(C) LinBox
 * Written by
 *  Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 *  Clément Pernet <Clement.Pernet@imag.fr>
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

 */

#ifndef __LINBOX_ring_interface_H
#define __LINBOX_ring_interface_H

namespace LinBox
{
	// LinBox Ring Interface
	/*
	 * The LinBox @link Rings ring@endlink common object @link Interfaces interface@endlink.
	 * The ring interface includes the following public members:
	 *
	 * Types: \c Element and \c RandIter.
	 *
	 * Object management member functions:
	 *   null constructor, copy constructor, destructor, assignment operator,
	 *   \c convert(), \c init(), \c assign(), \c characteristic(),\c cardinality().
	 *
	 * Predicates on ring elements:
	 *   \c areEqual(), \c isZero(), \c isOne().
	 *
	 * Basic arithmetic functions:
	 *   \c axpy(), \c add(), \c neg(), \c sub(), \c mul(), \c inv(), \c div().
	 *
	 * Inplace arithmetic functions:
	 *   \c axpyin(), \c addin(), \c negin(), \c subin(), \c mulin(), \c invin(), \c divin().
	 *
	 * I/O functions:
	 *   \c read() and \c write() for I/O of the ring itself and for I/O of its elements.
	 *
	 * The ring archetype class is is the reference instantiation of this
	 * interface and contains the generic specifications of the member functions.
	 * Documentation in other ring classes is more limited. It serves primarily to explain special properties
	 * specific to the class of the interface member functions and to explain any constructors
	 * or other functionality unique to the class.
	 *
	 *  @see Interfaces
	 */
	/**
	 * @brief This ring base class exists solely to aid documentation organization.

	 *  For the general ring member function documentation consult the
	 *  @link RingArchetype RingArchetype@endlink. For specific properties
	 *  of individual representations consult the specific ring classes.
	 */
	class RingInterface {
		/*
		   public:
		// this just demo's that some declarations could be here.
		typedef ElementArchetype Element;
		virtual Element& mul(Element& c, const Element& a, const Element& b) const = 0;
		*/
	};// empty class so doc++ makes a nice hierarchy.

} // namespace LinBox

#endif // __LINBOX_ring_interface_H

