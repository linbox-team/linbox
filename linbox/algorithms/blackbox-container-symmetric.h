/* linbox/algorithms/blackbox-container.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *
 * ------------------------------------
 * Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Modifications to incorporate into linbox-new tree
 * ------------------------------------
 *
 * 
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
 *.
 */

// ================================================================
// LinBox Project 1999
// Black Box iterator and container
// For symmetric matrix with same left and right vector
// the sequence is this->u^t this->v, this->u^t A this->v, ...,  this->u^t A^n this->v,
// Time-stamp: <25 Jan 02 16:04:24 Jean-Guillaume.Dumas@imag.fr>
// ================================================================

#ifndef __LINBOX_blackbox_container_symmetric_H
#define __LINBOX_blackbox_container_symmetric_H

#include "linbox/randiter/archetype.h"

#include "linbox/algorithms/blackbox-container.h"

namespace LinBox
{

	/// \brief See base class for doc.
	template<class Field, class _Blackbox, class RandIter = typename Field::RandIter>
	class BlackboxContainerSymmetric : public BlackboxContainerBase<Field, _Blackbox> {
	public:
		typedef _Blackbox Blackbox;

		BlackboxContainerSymmetric () {}

		template<class Vector>
		BlackboxContainerSymmetric (const Blackbox *D, const Field &F, const Vector &u0) :
			BlackboxContainerBase<Field, _Blackbox> (D, F)
		{ this->init (u0, u0); }
		BlackboxContainerSymmetric (const Blackbox *D, const Field &F, RandIter &g) :
			BlackboxContainerBase<Field, _Blackbox> (D, F)
		{ this->init (g); }

	protected:

		void _launch ()
		{
			if (this->casenumber > 0) {
				if (this->casenumber == 1) {
					this->casenumber = 2;
					this->_BB->apply (this->v, this->u);                // this->v <- B(B^i u_0) = B^(i+1) u_0
					this->_VD.dot (this->_value, this->u, this->v);     // t <- this->u^t this->v = u_0^t B^(2i+1) u_0
				}
				else {
					this->casenumber = -1;
					this->_VD.dot (this->_value, this->v, this->v);     // t <- this->v^t this->v = u_0^t B^(2i+2) u_0
				}
			}
			else {
				if (this->casenumber == 0) {
					this->casenumber = 1;
					this->_VD.dot (this->_value, this->u, this->u);     // t <- this->u^t this->u = u_0^t B^(2i+4) u_0
				}
				else {
					this->casenumber = 0;
					this->_BB->apply (this->u, this->v);                // this->u <- B(B^(i+1) u_0) = B^(i+2) u_0
					this->_VD.dot (this->_value, this->v, this->u);     // t <- this->v^t this->u = u_0^t B^(2i+3) u_0
				}
			}
		}

		void _wait () {}
	};

}

#endif // __LINBOX_blackbox_container_symmetric_H

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

