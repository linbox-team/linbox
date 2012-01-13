/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/blackbox-container-symmetrize.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *
 * ------------------------------------
 * Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Ported to linbox-new framework; replaced use of DOTPROD, etc. with
 * VectorDomain
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

/*! @internal
 * @file algorithms/blackbox-containter-symetrize.h
 * @ingroup algorithms
 * @brief Symmetrizing iterator for rank computations
 */

#ifndef __LINBOX_blackbox_container_symmetrize_H
#define __LINBOX_blackbox_container_symmetrize_H

#include "linbox/algorithms/blackbox-container-base.h"

namespace LinBox
{

	/** \brief Symmetrizing iterator (for rank computations).
	 * @ingroup algorithms
	 *
	 * Symmetrizing iterator (for rank computations)
	 * Same left and right vector
	 * A is supposed to have tranpose-vector product
	 * the sequence is
	 *
	 \code
	 this->u^t this->u
	 (A this->u)^t (A this->u) = this->u^t (A^t A) this->u
	 (A^t (A this->u))^t (A^t (A this->u)) = this->u^t (A^t A)^2 this->u
	 etc.
	 \endcode
	 */

	template<class Field, class _Blackbox, class RandIter = typename Field::RandIter>
	class BlackboxContainerSymmetrize : public BlackboxContainerBase<Field, _Blackbox> {
	public:

		typedef _Blackbox Blackbox;

		BlackboxContainerSymmetrize () {}

		template<class Vector>
		BlackboxContainerSymmetrize (const Blackbox *D, const Field &F, const Vector &u0) :
			BlackboxContainerBase<Field, Blackbox> (D, F)
		{
		       	init (u0);
	       	}

		//BlackboxContainerSymmetrize (const Blackbox *D, const Field &F, RandIter &g = typename Field::RandIter(_field) )
		BlackboxContainerSymmetrize (const Blackbox *D, const Field &F, RandIter &g = typename Field::RandIter() ) :
			BlackboxContainerBase<Field, Blackbox> (D, F)
		{
		       	init (g);
	       	}

	private:
		void _launch ()
		{
			if (this->casenumber) {
				this->casenumber = 0;
				this->_BB->apply (this->v, this->u);
				this->_VD.dot (this->_value, this->v, this->v);
			}
			else {
				this->casenumber = 1;
				this->_BB->applyTranspose (this->u, this->v);
				this->_VD.dot (this->_value, this->u, this->u);
			}
		}

		void _wait () {}
	};

}

#endif // __LINBOX_blackbox_container_symmetrize_H
