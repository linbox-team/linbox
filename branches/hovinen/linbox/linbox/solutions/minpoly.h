/* -*- mode: c; style: linux -*- */

/* linbox/solutions/minpoly.h
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

#ifndef __MINPOLY_H
#define __MINPOLY_H

#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"     // massey recurring sequence solver
#include "linbox/solutions/methods.h"

#include "linbox/field/archetype.h"
#include "linbox/blackbox/archetype.h"

namespace LinBox 
{

template <class Field, class Polynomial, class Vector, class Trait>
Polynomial &minpoly (Polynomial                       &P,
		     const Blackbox_archetype<Vector> &A,
		     const Field                      &F,
		     const Trait                      &M = Trait ())
{
	typename Field::RandIter i (F);
	unsigned long            deg;

	BlackboxContainer<Field, Vector> TF (&A, F, i);
	MasseyDomain< Field, BlackboxContainer<Field, Vector> > WD (&TF, M.Early_Term_Threshold ());

	WD.pseudo_minpoly (P, deg);

	return P;
}
 
}

#endif // __MINPOLY_H
