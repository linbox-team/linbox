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


#ifndef __LINBOX_diophantine_solver_H
#define __LINBOX_diophantine_solver_H

#include <linbox/algorithms/rational-solver.h>
#include <linbox/solutions/methods.h>
#include <linbox/blackbox/archetype.h>
#include <linbox/blackbox/lambda-sparse.h>
#include <linbox/blackbox/compose.h>

namespace LinBox 
{
	
	extern const char* solverReturnString[6] ;

	/** 
	 * \brief DiophantineSolver<QSolver> creates a diophantine solver using a QSolver to generate rational solutions

	 * Methods solve, randomSolve just expose functions from underlying rational solver.
	 *       Method diophantineSolve creates a solution with minimal denominator, and can also create
	 *       a certificate of minimality (described in 'Certified Dense Linear System Solving' by Mulders+Storjohann)
	 *       which will be left in the public field lastCertificate.
	 */
	template<class QSolver>
	class DiophantineSolver {
		
	protected:

		typedef typename QSolver::RingType    Ring;
		typedef typename Ring::Element        Integer;
		QSolver&                              _rationalSolver;
		Ring                                  _R;
		Integer                               _rone;

	public:
		// information for last diophantine solve
		mutable int numSolutionsNeeded;
		mutable int numFailedCallsToSolver;
		mutable int numRevelantSolutions;

		VectorFraction<Ring>                  lastCertificate;

		/* Constructor from a rationalSolver
		 * @param rs  , a rationalSolver
		 */
		DiophantineSolver (QSolver& rs) :
			_rationalSolver(rs), _R(rs.getRing()), lastCertificate(_R, 0) {
			_R.init(_rone, 1);
		};

		/** Solve a linear system Ax=b over quotient field of a ring
		 * 
		 * @param A        , Matrix of linear system
		 * @param x        , Vector in which to store solution
		 * @param b        , Right-hand side of system
		 * @param maxPrimes, maximum number of moduli to try
		 * @param level    , level of certification to be used
		 *
		 * @return status of solution. if (return != SS_FAILED), and (level >= SL_LASVEGAS), solution is guaranteed correct.
		 *   SS_FAILED - all primes used were bad
		 *   SS_OK - solution found. 
		 *   SS_INCONSISTENT - system appreared inconsistent. certificate is in lastCertificate if (level >= SL_CERTIFIED)
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solve(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES, 
					 const SolverLevel level = SL_DEFAULT);

		/** Find a random solution of the general linear system Ax=b over quotient field of a ring.
		 * 
		 * @param A        , Matrix of linear system
		 * @param x        , Vector in which to store solution
		 * @param b        , Right-hand side of system
		 * @param maxPrimes, maximum number of moduli to try
		 * @param level    , level of certification to be used
		 *
		 * @return status of solution. if (return != SS_FAILED), and (level >= SL_LASVEGAS), solution is guaranteed correct.
		 *   SS_FAILED - all primes used were bad
		 *   SS_OK - solution found. 
		 *   SS_INCONSISTENT - system appreared inconsistent. certificate is in lastCertificate if (level >= SL_CERTIFIED)
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus randomSolve(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES, 
					       const SolverLevel level = SL_DEFAULT);
 
		/** 
		 * Find a solution of the linear system Ax=b whose denominator (when written as an integer vector over a single denom) is minimal.
		 *
		 * @param A        , Matrix of linear system
		 * @param x        , Vector in which to store solution
		 * @param b        , Right-hand side of system
		 * @param maxPrimes, maximum number of moduli to try
		 * @param level    , level of certification to be used
		 *
		 * @return status of solution. if (return != SS_FAILED) and (level >= SL_LASVEGAS), solution is guaranteed correct
		 *                             if (return == SS_OK) and (level >= SL_LASVEGAS), solution is guaranteed minimal.
		 *   SS_FAILED - all primes used were bad
		 *   SS_OK - solution found. certificate of minimality is in lastCertificate if (level >= SL_CERTIFIED)
		 *   SS_INCONSISTENT - system appreared inconsistent. certificate of inconsistency is in lastCertificate if (level >= SL_CERTIFIED)
		 *
		 * @return status of solution - OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER
		 */	
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus diophantineSolve(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES, 
						    const SolverLevel level = SL_DEFAULT);
	};

}
#include <linbox/algorithms/diophantine-solver.inl>

#endif //__LINBOX_diophantine_solver_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
