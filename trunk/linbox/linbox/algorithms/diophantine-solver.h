/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/algorithms/lifting-container.h
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
	
	enum SolverReturnStatus
	{
			OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER, BAD_PRIME, OK_NOT_DIOPHANTINE
	};

	const char* solverReturnString[] 
	= {"OK", "FAILED", "SINGULAR", "INCONSISTENT", "BAD_PRECONDITIONER", "BAD_PRIME", "OK_NOT_DIOPHANTINE"};
		
	/** _Ring integer ring
	 *  _Field, finite field for lifting
	 */

	template<class RationalSolver>
	class DiophantineSolver {
	protected:
		RationalSolver _rationalSolver;
		
	public:
		/* Constructor from a rationalSolver
		 * @param rs  , a rationalSolver
		 */
		DiophantineSolver (const RationalSolver& rs = RationalSolver()):
			_rationalSolver(rs){};

		
		/** Find a solution of the linear system Ax=b over quotient field of a ring.
		 * Calls solve from RationalSolver so it is deterministic.
		 * Should remove ability to return 'bad_preconditioner'
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
		 * Should remove ability to return 'bad_preconditioner'
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
		 * If no diophantine solution was found
		 * Should remove ability to return 'bad_preconditioner'
		 *
		 * @param A   , Matrix of linear system
		 * @param x   , Vector in which to store solution
		 * @param b   , Right-hand side of system
		 *
		 * @return status of solution - OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER, OK_NOT_DIOPHANTINE
		 */	
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus diophantineSolve(Vector1& x, const IMatrix& A, const Vector2& b, int maxPrimes = DEFAULT_MAXPRIMES);
				
	};

}
#include <linbox/algorithms/diophantine-solver.inl>

#endif
