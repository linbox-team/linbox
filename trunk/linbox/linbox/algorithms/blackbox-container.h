/* -*- mode: c; style: linux -*- */

/* linbox/algorithms/blackbox-container.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __BLACKBOX_CONTAINER_H
#define __BLACKBOX_CONTAINER_H

#include <linbox/randiter/archetype.h>
#include <linbox/algorithms/blackbox-container-base.h>

namespace LinBox 
{

template<class Field, class Vector, class RandIter = typename Field::RandIter>
class BlackboxContainer : public BlackboxContainerBase<Field, Vector> {
    public:
	BlackboxContainer () {} 

	BlackboxContainer(const Blackbox * D, const Field &F, const Vector &u0) 
		: BlackboxContainerBase<Field, Vector> (D, F)
		{ init (u0, u0); w = u; }
    
	BlackboxContainer(const Blackbox * D, const Field &F, const Vector &u0, const Vector& v0) 
		: BlackboxContainerBase<Field, Vector> (D, F)
		{ init (u0, v0); w = u;}
    
	BlackboxContainer(const Blackbox * D, const Field &F, RandIter &g) 
		: BlackboxContainerBase<Field, Vector> (D, F)
		{ init (g); w = u; }
    
    protected:
	Vector w;

	void _launch () {
		if (even) {
			_BB->apply (v, w);  // GV
			_VD.dot (_value, u, v);  // GV 
			even = 0;
		} else {
			_BB->apply (w, v);  // GV
			_VD.dot (_value, u, w);  // GV
			even = 1;
		}  
	}

	void _wait () {}
};
 
}

#endif // __BLACKBOX_CONTAINER_H

