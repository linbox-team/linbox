/* -*- mode: c; style: linux -*- */

/* linbox/solutions/det.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

#ifndef __DET_H
#define __DET_H

#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/algorithms/massey-domain.h"

#include "linbox/vector/vector-traits.h"

#include "linbox/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** Compute the rank of a linear operator A, represented as a black box
	 */

	template <class Field, class Vector>
	typename Field::element &det (typename Field::element          &res,
				      const Blackbox_archetype<Vector> &A,
				      const Field                      &F,
				      const MethodTrait::Wiedemann     &M = MethodTrait::Wiedemann ()) 
	{
		typedef std::vector<typename Field::element> Polynomial;

		linbox_check (A.coldim () == A.rowdim ());

		Polynomial      phi;
		unsigned long   deg;
		Field::RandIter iter (F);

		// Precondition here to separate the eigenvalues, so that
		// minpoly (B) = charpoly (B) with high probability

		Vector d (A.coldim ()), dinv (A.coldim ());
		int i;

		for (i = 0; i < A.coldim (); i++) {
			do iter.random (d[i]) while (F.isZero (d[i]));
			F.inv (dinv[i], d[i]);
		}
	
		Diagonal<Field, Vector> D (d), Dinv (dinv);

		Compose<Vector> B(Compose (Compose (Dinv, A), D1);
		BlackboxContainer<Field, Vector> TF (&B, F, iter);
		MasseyDomain<Field, BlackboxContainer<Field, Vector> > WD (&TF, M.Early_Term_Threshold ());

		WD.pseudo_minpoly (phi, deg);

		F.assign (res, phi[0]);

		return res;
	}
}

#endif // __DET_H
