/* -*- mode: c++; style: linux -*- */

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
 * See COPYING for license information.
 */

// ================================================================
// LinBox Project 1999
// Black Box iterator and container 
// For symmetric matrix with same left and right vector
// the sequence is u^t v, u^t A v, ...,  u^t A^n v,  
// Time-stamp: <25 Jan 02 16:04:24 Jean-Guillaume.Dumas@imag.fr> 
// ================================================================

#ifndef __BLACKBOX_CONTAINER_SYMMETRIC_H
#define __BLACKBOX_CONTAINER_SYMMETRIC_H

#include "linbox/field/archetype.h"
#include "linbox/randiter/archetype.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/algorithms/blackbox-container.h"

namespace LinBox
{

template<class Field, class Vector, class RandIter = typename Field::RandIter>
class BlackboxContainerSymmetric : public BlackboxContainerBase<Field, Vector>
{
    public:
	BlackboxContainerSymmetric () {} 
	BlackboxContainerSymmetric (BlackboxArchetype<Vector> *D, const Vector &u0)
		: BlackboxContainerBase<Field, Vector> (D)
		{ init (u0, u0); }
	BlackboxContainerSymmetric (BlackboxArchetype<Vector> *D, RandIter &g)
		: BlackboxContainerBase<Field, Vector> (D)
		{ init (g); }

    protected:

	void _launch () {
		if (even > 0) {
			if (even == 1) {
				even = 2;
				_BB->Apply (v, u);          // v <- B(B^i u_0) = B^(i+1) u_0
				_VD.dot (_value, u, v);     // t <- u^t v = u_0^t B^(2i+1) u_0
			} else {
				even = -1;
				_VD.dot (_value, v, v);     // t <- v^t v = u_0^t B^(2i+2) u_0
			}
		} else {
			if (even == 0) {
				even = 1;
				_VD.dot (_value, u, u);     // t <- u^t u = u_0^t B^(2i+4) u_0
			} else {
				even = 0;
				_BB->apply (u, v);          // u <- B(B^(i+1) u_0) = B^(i+2) u_0
				_VD.dot (_value, v, u);     // t <- v^t u = u_0^t B^(2i+3) u_0
			}   
		}
	}
    
	void _wait () {}
};
 
}

#endif // __BBContainer_SYMMETRIC_H__
