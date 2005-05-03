/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#ifndef __BLACKBOX_CONTAINER_SYMMETRIZE_H
#define __BLACKBOX_CONTAINER_SYMMETRIZE_H

#include <linbox/algorithms/blackbox-container-base.h>

namespace LinBox 
{

/** \brief Symmetrizing iterator (for rank computations).
@doc
 #
//================================================================
// LinBox Project 1999
// Symmetrizing iterator (for rank computations)
// Same left and right vector
// A is supposed to have tranpose-vector product
// the sequence is u^t u, (A u)^t (A u) = u^t (A^t A) u, 
// (A^t (A u))^t (A^t (A u)) = u^t (A^t A)^2 u , etc.
// Time-stamp: <13 Jun 02 18:16:43 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================
 #
 */

	template<class Field, class _Blackbox, class RandIter = typename Field::RandIter>
	class BlackboxContainerSymmetrize : public BlackboxContainerBase<Field, _Blackbox> {
		using BlackboxContainerBase<Field, _Blackbox>::u; 
		using BlackboxContainerBase<Field, _Blackbox>::v; 
		using BlackboxContainerBase<Field, _Blackbox>::_BB; 
		using BlackboxContainerBase<Field, _Blackbox>::_VD; 
		using BlackboxContainerBase<Field, _Blackbox>::casenumber; 
		using BlackboxContainerBase<Field, _Blackbox>::_value; 
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
			if (casenumber) {
				casenumber = 0;
				_BB->apply (v, u);
				_VD.dot (_value, v, v); 
			} else {
				casenumber = 1;
				_BB->applyTranspose (u, v); 
				_VD.dot (_value, u, u);
			}
		}

		void _wait () {}
	};
 
};

#endif // __BLACKBOX_CONTAINER_SYMMETRIZE_H
