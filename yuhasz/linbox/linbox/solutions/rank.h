/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/rank.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * See COPYING for license information.
 */

#ifndef __RANK_H
#define __RANK_H

#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/algorithms/blackbox-container-symmetrize.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/gauss.h"

#include "linbox/vector/vector-traits.h"
#include "linbox/solutions/methods.h"

#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** Compute the rank of a linear operator A, represented as a black box
	 */

	template <class Field, class Blackbox>
	unsigned long &rank (unsigned long                   &res,
			     const Blackbox                  &A,
			     const Field                     &F,
			     const MethodTrait::Wiedemann    &M = MethodTrait::Wiedemann ()) 
	{
		typename Field::RandIter iter (F);

		commentator.start ("Rank", "rank");

		std::vector<typename Field::Element> d1, d2;
		size_t i;

		VectorWrapper::ensureDim (d1, A.coldim ());
		VectorWrapper::ensureDim (d2, A.rowdim ());

		for (i = 0; i < A.coldim (); i++)
			do iter.random (d1[i]); while (F.isZero (d1[i]));

		for (i = 0; i < A.rowdim (); i++)
			do iter.random (d2[i]); while (F.isZero (d2[i]));

		Diagonal<Field> D1 (F, d1), D2 (F, d2);
		Transpose<Blackbox> AT (&A);

		Compose<Diagonal<Field>,Transpose<Blackbox> > B1 (&D1, &AT);
		Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> > B2 (&B1, &D2);
		Compose<Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> >, Blackbox> B3 (&B2, &A);
		Compose<Compose<Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> >, Blackbox>, Diagonal<Field> > B (&B3, &D1);
                    // JGD 22.03.03
// 		BlackboxContainer<Field, Vector> TF (&B, F, iter);
// 		MasseyDomain<Field, BlackboxContainer<Field, Vector> > WD (&TF, M.earlyTermThreshold ());

		typedef Compose<Compose<Compose<Compose<Diagonal<Field>,Transpose<Blackbox> >, Diagonal<Field> >, Blackbox>, Diagonal<Field> > Blackbox1;
		BlackboxContainerSymmetric<Field, Blackbox1> TF (&B, F, iter);
		MasseyDomain<Field, BlackboxContainerSymmetric<Field, Blackbox1> > WD (&TF, M.earlyTermThreshold ());

                    // Here there is an extra diagonal computation
                    // The probability of success is also divided by two, as 
                    // D2^2 contains only squares and squares are half the total elements
// 		BlackboxContainerSymmetrize<Field, Vector> TF (&B2, F, iter);
// 		MasseyDomain<Field, BlackboxContainerSymmetrize<Field, Vector> > WD (&TF, M.earlyTermThreshold ());

		WD.pseudo_rank (res);

		commentator.stop ("done", NULL, "rank");

		return res;
	}

// 	template <class Field, class Blackbox>
// 	unsigned long &rank (unsigned long                   &res,
// 			     const Blackbox                  &A,
// 			     const Field                     &F,
// 			     const MethodTrait::Wiedemann    &M = MethodTrait::Wiedemann ()) 
// 	{
// 		return rank<Field, Blackbox, std::vector<typename Field::Element> > (res, A, F, M);
// 	}


	template <class Field, class Matrix>
	unsigned long &rank (unsigned long                   &res,
			     const Matrix                    &A,
			     const Field                     &F,
			     const MethodTrait::Elimination  &M) 
	{
		commentator.start ("Rank", "rank");

		GaussDomain<Field> GD (F);

		Matrix A1 (A);   // We make a copy as these data will be destroyed

		GD.rankin (res, A1, M.strategy ());
                
		commentator.stop ("done", NULL, "rank");
                
		return res;
	}
    
	template <class Field, class Matrix>
	unsigned long &rankin (unsigned long                   &res,
			             Matrix                    &A,
			       const Field                     &F,
			       const MethodTrait::Elimination  &M) 
	{
		commentator.start ("Rank", "rank");

		GaussDomain<Field> GD (F);

                GD.rankin( res, A, M.strategy ());

		commentator.stop ("done", NULL, "rank");

		return res;
	}
}

#endif // __RANK_H
