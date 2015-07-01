/* linbox/algorithms/diophantine-solver.h
 * Copyright (C) 2004 David Pritchard
 *
 * Written by David Pritchard <daveagp@mit.edu>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
  * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#ifndef __LINBOX_diophantine_solver_H
#define __LINBOX_diophantine_solver_H

#include "linbox/algorithms/rational-solver.h"
#include "linbox/solutions/methods.h"
#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/lambda-sparse.h"
#include "linbox/blackbox/compose.h"

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
		Ring                                  _ring;

	public:
		// information for last diophantine solve
		mutable int numSolutionsNeeded;
		mutable int numFailedCallsToSolver;
		mutable int numRevelantSolutions;

		VectorFraction<Ring>                  lastCertificate;

		/*! Constructor from a rationalSolver
		 * @param rs  a rationalSolver
		 */
		DiophantineSolver (QSolver& rs) :
			_rationalSolver(rs), _ring(rs.getRing()), lastCertificate(_ring, 0)
		{ }

		/** Solve a linear system \c Ax=b over quotient field of a ring.
		 *
		 * @param A        Matrix of linear system
		 * @param x        Vector in which to store solution
		 * @param b        Right-hand side of system
		 * @param maxPrimes maximum number of moduli to try
		 * @param level    level of certification to be used
		 * @param den
		 *
		 * @return status of solution. if \c (return != SS_FAILED), and \c (level >= SL_LASVEGAS), solution is guaranteed correct.
		 *   \c SS_FAILED - all primes used were bad
		 *   \c SS_OK - solution found.
		 *   \c SS_INCONSISTENT - system appreared inconsistent. certificate is in \p lastCertificate if \c (level >= SL_CERTIFIED)
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus solve(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES,
					 const SolverLevel level = SL_DEFAULT);

		/** Find a random solution of the general linear system \c Ax=b over quotient field of a ring.
		 *
		 * @param A        Matrix of linear system
		 * @param x        Vector in which to store solution
		 * @param b        Right-hand side of system
		 * @param maxPrimes maximum number of moduli to try
		 * @param level    level of certification to be used
		 * @param den
		 *
		 * @return status of solution. if \c (return != SS_FAILED), and \c (level >= SL_LASVEGAS), solution is guaranteed correct.
		 *  \c  SS_FAILED - all primes used were bad
		 *  \c  SS_OK - solution found.
		 *  \c  SS_INCONSISTENT - system appreared inconsistent. certificate is in \p lastCertificate if \c (level >= SL_CERTIFIED)
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus randomSolve(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES,
					       const SolverLevel level = SL_DEFAULT);

		/**
		 * Find a solution of the linear system \c Ax=b whose denominator (when written as an integer vector over a single denom) is minimal.
		 *
		 * @param A        Matrix of linear system
		 * @param x        Vector in which to store solution
		 * @param b        Right-hand side of system
		 * @param maxPrimes maximum number of moduli to try
		 * @param level     level of certification to be used
		 * @param den
		 *
		 * @return status of solution. if \c (return != SS_FAILED) and \c  (level >= SL_LASVEGAS), solution is guaranteed correct
		 *                             if \c (return == SS_OK) and \c (level >= SL_LASVEGAS), solution is guaranteed minimal.
		 *   \c SS_FAILED - all primes used were bad
		 *   \c SS_OK - solution found. certificate of minimality is in lastCertificate if \c (level >= SL_CERTIFIED)
		 *  \c SS_INCONSISTENT - system appreared inconsistent. certificate of inconsistency is in \p lastCertificate if \c (level >= SL_CERTIFIED)
		 *
		 * @return status of solution - OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER
		 */
		template<class IMatrix, class Vector1, class Vector2>
		SolverReturnStatus diophantineSolve(Vector1& x, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes = DEFAULT_MAXPRIMES,
						    const SolverLevel level = SL_DEFAULT);
	};

}
#include "linbox/algorithms/diophantine-solver.inl"

#endif //__LINBOX_diophantine_solver_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

