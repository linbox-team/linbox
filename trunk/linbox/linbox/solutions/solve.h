/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/solve.h
 * Copyright (C) 2002 Zhendong Wan,
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>
 *
 * ------------------------------------
 * 2002-08-09  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed from solver.h to solve.h, tweak indentation,
 * Add argument SolverTraits for additional parameters,
 * Add header
 * ------------------------------------
 * See COPYING for license information.
 */

#ifndef __SOLVE_H
#define __SOLVE_H

#include <vector>
#include <algorithm>

#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/massey-domain.h" 
#include "linbox/util/debug.h"
#include "linbox/field/vector-domain.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{

template <class Field, class BlackBox, class Vector>
Vector& solve (const BlackBox     &A,
	       Vector             &x,		       
	       const Vector       &y,
	       const Field        &F,
	       const SolverTraits &traits = SolverTraits ())
{
	typedef std::vector<typename Field::Element> Polynomial;

	linbox_check ((x.size () == A.coldim ()) &&
		      (y.size () == A.rowdim ()));

	typename Field::RandIter randiter(F);
	Polynomial               P;
	unsigned long            deg;

	BlackboxContainer<Field, Vector> TF (&A, F, y);
	MasseyDomain< Field, BlackboxContainer<Field, Vector> > WD (&TF);
		
	WD.minpoly (P, deg);
		
	Polynomial::iterator p = P.begin ();

	for (++p; p != P.end (); ++p) {
		F.divin (*p, P[0]);
		F.negin (*p);
	}
		
	std::copy (P.begin () + 1, P.end (), P.begin ());
	P.resize (P.size () - 1);
		
	if (P.empty ())
		return x = y;

	VectorDomain <Field> VD (F);		
	int                  i;

	VD.mul (x, y, P[P.size () - 1]);

	for (i = P.size () - 2; i >= 0; i--) {
		A.applyIn (x);
		VD.axpyin (x, P[i], y);
	}

	return x;
}								       		

}

#endif // __SOLVE_H
