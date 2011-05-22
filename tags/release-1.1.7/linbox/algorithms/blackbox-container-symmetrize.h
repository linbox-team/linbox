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
 * See COPYING for license information.
 */

#ifndef __LINBOX_blackbox_container_symmetrize_H
#define __LINBOX_blackbox_container_symmetrize_H

#include <linbox/algorithms/blackbox-container-base.h>

namespace LinBox 
{

/** \brief Symmetrizing iterator (for rank computations).

 #
//================================================================
// LinBox Project 1999
// Symmetrizing iterator (for rank computations)
// Same left and right vector
// A is supposed to have tranpose-vector product
// the sequence is this->u^t this->u, (A this->u)^t (A this->u) = this->u^t (A^t A) this->u, 
// (A^t (A this->u))^t (A^t (A this->u)) = this->u^t (A^t A)^2 this->u , etc.
// Time-stamp: <13 Jun 02 18:16:43 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================
 #
 */

	template<class Field, class _Blackbox, class RandIter = typename Field::RandIter>
	class BlackboxContainerSymmetrize : public BlackboxContainerBase<Field, _Blackbox> {
	    public:

                typedef _Blackbox Blackbox;

		BlackboxContainerSymmetrize () {} 

		template<class Vector>
		BlackboxContainerSymmetrize (const Blackbox *D, const Field &F, const Vector &u0) 
			: BlackboxContainerBase<Field, Blackbox> (D, F) { init (u0); }
    
		//BlackboxContainerSymmetrize (const Blackbox *D, const Field &F, RandIter &g = typename Field::RandIter(_F) ) 
		BlackboxContainerSymmetrize (const Blackbox *D, const Field &F, RandIter &g = typename Field::RandIter() ) 
			: BlackboxContainerBase<Field, Blackbox> (D, F) { init (g); }

	    private:
		void _launch () {
			if (this->casenumber) {
				this->casenumber = 0;
				this->_BB->apply (this->v, this->u);
				this->_VD.dot (this->_value, this->v, this->v); 
			} else {
				this->casenumber = 1;
				this->_BB->applyTranspose (this->u, this->v); 
				this->_VD.dot (this->_value, this->u, this->u);
			}
		}

		void _wait () {}
	};
 
}

#endif // __LINBOX_blackbox_container_symmetrize_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
