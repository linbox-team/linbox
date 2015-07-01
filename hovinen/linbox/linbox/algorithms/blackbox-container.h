/* -*- mode: c; style: linux -*- */

/* linbox/src/library/objects/algorithms/blackbox/blackbox-container.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas, Bradford Hovinen
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __BLACKBOX_CONTAINER_H
#define __BLACKBOX_CONTAINER_H

#include <linbox/randiter/archetype.h>
#include <linbox/algorithms/blackbox-container-base.h>

namespace LinBox 
{

template<class Field, class Vector>
class BlackboxContainer : public BlackboxContainerBase<Field, Vector> {
public:
	typedef typename Field::RandIter RandIter;

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
			DOTPROD (_value, u, v);  // GV 
			even = 0;
		} else {
			_BB->apply (w, v);  // GV
			DOTPROD (_value, u, w);  // GV
			even = 1;
		}  
	}

	void _wait () {}
};
 
}

#endif // __BLACKBOX_CONTAINER_H

