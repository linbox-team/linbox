/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/diophantine-solver.h
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


#ifndef __LINBOX_DIOPHANTINE_SOLVER_H
#define __LINBOX_DIOPHANTINE_SOLVER_H

#include <linbox/algorithms/rational-solver.h>
#include <linbox/solutions/methods.h>
#include <linbox/blackbox/archetype.h>
#include <linbox/blackbox/lambda-sparse.h>
#include <linbox/blackbox/compose.h>

namespace LinBox {
	
	const char* solverReturnString[] 
	= {"OK", "FAILED", "SINGULAR", "INCONSISTENT", "BAD_PRECONDITIONER", "BAD_PRIME"};

	/** 
	 * @memo DiophantineSolver<QSolver> creates a diophantine solver using a QSolver to generate rational solutions
	 * @doc  Methods solve, randomSolve just expose functions from underlying rational solver.
	 *       Method diophantineSolve creates a solution with minimal denominator, and can also create
	 *       a certificate of minimality (described in 'Certified Dense Linear System Solving' by Mulders+Storjohann)
	 *       which will be left in the public field lastCertificate.
	 */
	template<class QSolver>
	class DiophantineSolver {
		
	protected:

		typedef typename QSolver::RingType    Ring;
		typedef typename Ring::Element        Integer;
		QSolver                               _rationalSolver;
		Ring                                  _R;
		Integer                               _rone;
		
	public:
		VectorFraction<Ring>                  lastCertificate;

		/* Constructor from a rationalSolver
		 * @param rs  , a rationalSolver
		 */
		DiophantineSolver (QSolver rs) :
			_rationalSolver(rs), _R(rs.getRing()), lastCertificate(_R, 0) {
			_R.init(_rone, 1);
		};

		/** Find a solution of the linear system Ax=b over quotient field of a ring.
		 * Calls solve from QSolver so it is deterministic.
		 *
		 * @param A   , Matrix of linear system
		 * @param x   , Vector in which to store solution
		 * @param b   , Right-hand side of system
		 *
		 * @return status of solution - OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER
		 */	
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solve(Vector1& x, const IMatrix& A, const Vector2& b, int maxPrimes = DEFAULT_MAXPRIMES);				
		/** Find a random solution of the linear system Ax=b over quotient field of a ring.
		 *
		 * @param A   , Matrix of linear system
		 * @param x   , Vector in which to store solution
		 * @param b   , Right-hand side of system
		 *
		 * @return status of solution - OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER
		 */	
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus randomSolve(Vector1& x, const IMatrix& A, const Vector2& b, int maxPrimes = DEFAULT_MAXPRIMES);
 
		/** Find a diophantine solution of the linear system Ax=b over quotient field of a ring.
		 * If no diophantine solution exists, one with minimal denominator is returned.
		 *
		 * @param A   , Matrix of linear system
		 * @param x   , Vector in which to store solution
		 * @param b   , Right-hand side of system
		 *
		 * @return status of solution - OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER
		 */	
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus diophantineSolve(Vector1& x, const IMatrix& A, const Vector2& b, bool makeCertificate = false,
						    int maxPrimes = DEFAULT_MAXPRIMES);

					
	};

}
#include <linbox/algorithms/diophantine-solver.inl>

#endif
