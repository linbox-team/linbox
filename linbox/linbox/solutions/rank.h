/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/rank.h
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

#ifndef __RANK_H
#define __RANK_H

#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"

#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/methods.h"

#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** Compute the rank of a linear operator A, represented as a black box
	 */

	template <class Field, class Vector>
	unsigned long &rank (unsigned long                    &res,
			     const BlackboxArchetype<Vector> &A,
			     const Field                      &F,
			     const MethodTrait::Wiedemann     &M = MethodTrait::Wiedemann ()) 
	{
		typename Field::RandIter iter (F);

		commentator.start ("Rank", "rank");

		Vector d1, d2;
		size_t i;

		VectorWrapper::ensureDim (d1, A.coldim ());
		VectorWrapper::ensureDim (d2, A.rowdim ());

		for (i = 0; i < A.coldim (); i++)
			do iter.random (d1[i]); while (F.isZero (d1[i]));

		for (i = 0; i < A.rowdim (); i++)
			do iter.random (d2[i]); while (F.isZero (d2[i]));

		Diagonal<Field, Vector> D1 (F, d1), D2 (F, d2);
		Transpose<Vector> AT (&A);

		Compose<Vector> B1 (&D1, &AT);
		Compose<Vector> B2 (&B1, &D2);
		Compose<Vector> B3 (&B2, &A);
		Compose<Vector> B (&B3, &D1);
		BlackboxContainer<Field, Vector> TF (&B, F, iter);
		MasseyDomain<Field, BlackboxContainer<Field, Vector> > WD (&TF, M.earlyTermThreshold ());

		WD.pseudo_rank (res);

		commentator.stop ("done", NULL, "rank");

		return res;
	}
}

#endif // __RANK_H
