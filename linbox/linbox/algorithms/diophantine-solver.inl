/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/rational-solver.inl
 * Copyright (C) 2004 David Pritchard
 *
 * Written by David Pritchard <daveagp@mit.edu>
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

#ifndef __LINBOX_DIOPHANTINE_SOLVER_INL
#define __LINBOX_DIOPHANTINE_SOLVER_INL

#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/lambda-sparse.h>
#include <linbox/algorithms/rational-solver.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/debug.h>

#include <linbox-config.h>

namespace LinBox {


	template<class RationalSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<RationalSolver>::solve (Vector1& x,
								     const IMatrix& A,
								     const Vector2& b, 
								     int maxPrimes = DEFAULT_MAXPRIMES) {
		
		SolverReturnStatus status=FAILED;
		status = (SolverReturnStatus)(int)_rationalSolver.solve(x, A, b, false, maxPrimes);
		return status;
	}
	
	template<class RationalSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<RationalSolver>::randomSolve (Vector1& x,
									   const IMatrix& A,
									   const Vector2& b, 
									   int maxPrimes = DEFAULT_MAXPRIMES) {
		
		SolverReturnStatus status=FAILED;
		status = (SolverReturnStatus)(int)_rationalSolver.findRandomSolution(x, A, b, maxPrimes);
		return status;
	}
	
	template<class RationalSolver>
	template<class IMatrix, class Vector1, class Vector2>	
	SolverReturnStatus DiophantineSolver<RationalSolver>::diophantineSolve (Vector1& x,
										const IMatrix& A,
										const Vector2& b, 
										int maxPrimes = DEFAULT_MAXPRIMES) {
		
		SolverReturnStatus status=FAILED;
		return status;
	}
	
} //end of namespace LinBox
#endif
