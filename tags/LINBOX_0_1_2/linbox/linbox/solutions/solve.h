/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/solve.h
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-08-09  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed from solver.h to solve.h, tweak indentation,
 * Add argument SolverTraits for additional parameters,
 * Add header
 *
 * Moved solver code proper to linbox/algorithms/wiedemann.h; solve () is now
 * just a wrapper (consequently copyright notices have changed)
 * ------------------------------------
 * See COPYING for license information.
 */

#ifndef __SOLVE_H
#define __SOLVE_H

#include <vector>
#include <algorithm>

#include "linbox/algorithms/wiedemann.h"
#include "linbox/util/debug.h"
#include "linbox/field/vector-domain.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{

/** Solve a system Ax=b over the field F
 */

template <class Field, class BlackBox, class Vector>
Vector &solve (const BlackBox     &A,
	       Vector             &x,		       
	       const Vector       &b,
	       const Field        &F,
	       const SolverTraits &traits = SolverTraits ())
{
	switch (traits.method ()) {
	    case SolverTraits::METHOD_WIEDEMANN:
		return solveWiedemann (A, x, b, F, traits);

	    case SolverTraits::METHOD_LANCZOS:
		throw LinboxError ("Lanczos-based solver not implemented");

	    case SolverTraits::METHOD_ELIMINATION:
		throw LinboxError ("Elimination-based solver not implemented");
	}

	return x;
}

}

#endif // __SOLVE_H
