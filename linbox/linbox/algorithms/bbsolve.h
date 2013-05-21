/* Copyright (c) LinBox
 * linbox/algorithms/bbsolve.h
 * (was linbox/solutions/solve.h)
 * written
 *  by Bradford Hovinen <hovinen@cis.udel.edu>
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_bbsolve_H
#define __LINBOX_bbsolve_H

#include <algorithm>

// must fix this list...
#include "linbox/algorithms/wiedemann.h"
#include "linbox/algorithms/lanczos.h"
#include "linbox/algorithms/mg-block-lanczos.h"
#include "linbox/util/debug.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/solutions/methods.h"

namespace LinBox
{

#if 0
	/** @name Solvers
	 * @brief Solving linear system Ax = b over the field F.
	 */
	@{
#endif

#if 0
	// for specialization on method.
	template <class Field, class Vector, class Blackbox, class MyMethod>
	Vector& solve (const Blackbox&A,
		       Vector                          &x,
		       const Vector                    &b,
		       const Field                     &F,
		       const MyMethod & m );
#endif

	/** Solve Ax=b over field F using Wiedemann's method, with inconsistency certificate.
	 *
	 * This is a general interface for the linear system solving
	 * capabilities of LinBox. If the system is nonsingular, it returns
	 * the unique solution, storing it in the vector x. If the system is
	 * consistent and singular, it returns a random solution. Repeated
	 * calls to this function can give a complete description of the
	 * solution manifold. If the system is inconsistent and the
	 * \ref SolverTraits structure supplied requests certification of
	 * inconsistency, it fills in the certificate of
	 * inconsistency. Otherwise, it runs through the number of iterations
	 * specified in \c traits and throws a \ref SolveFailed exception
	 * if it cannot find a solution.
	 *
	 * This specialization uses Wiedemann's algorithm and is the default.
	 *
	 * @param A Black box matrix of the system
	 * @param x Place to store solution vector
	 * @param b Right-hand side
	 * @param u Vector in which to store certificate of inconsistency, if required
	 * @param F Field over which to perform computations
	 * @param traits \ref SolverTraits structure with user-specified parameters
	 * @return Reference to solution vector
	 */

	template <class Field, class Vector, class Blackbox>
	typename WiedemannSolver<Field>::ReturnStatus
	solve (const Blackbox                  &A,
	       Vector                          &x,
	       const Vector                    &b,
	       Vector                          &u,
	       const Field                     &F,
	       const WiedemannTraits &traits = WiedemannTraits ())
	{
		WiedemannSolver<Field> solver (F, traits);
		return solver.solve (A, x, b, u);
	}

	/** Solve Ax=b over field F using the Wiedemann method.
	 *
	 * This version differs from the one above in that there is no extra
	 * parameter for the certificate of inconsistency, and it throws
	 * exceptions if the solution fails. It also returns a reference to
	 * the solution vector.
	 */

	template <class Field, class Vector, class Blackbox>
	Vector &solve (const Blackbox                  &A,
		       Vector                          &x,
		       const Vector                    &b,
		       const Field                     &F,
		       const WiedemannTraits &traits = WiedemannTraits ())
	{
		Vector u(A.field());
		WiedemannSolver<Field> solver (F, traits);

		VectorWrapper::ensureDim (u, A.rowdim ());

		switch (solver.solve (A, x, b, u)) {
		case WiedemannSolver<Field>::OK:
			return x;

		case WiedemannSolver<Field>::FAILED:
			throw SolveFailed ();

		case WiedemannSolver<Field>::SINGULAR:
			throw SolveFailed ();

		case WiedemannSolver<Field>::INCONSISTENT:
			throw InconsistentSystem<Vector> (u);

		default:
			throw LinboxError ("Bad return value from solve");
		}
	}

	/** Solve Ax=b over field F using the Lanczos method.
	 *
	 * This is a general interface for the linear system solving capabilities of
	 * LinBox. If the system is nonsingular, it returns the unique solution, storing
	 * it in the vector x. If the system is consistent and singular, it returns a
	 * random solution. Repeated calls to this function can give a complete
	 * description of the solution manifold. If the system is inconsistent and the
	 * \ref SolverTraits  structure has result checking turned on, it runs through
	 * the number of iterations specified in \c traits and throws a
	 * \ref SolveFailed  exception if it cannot find a solution.
	 *
	 * This specialization uses the Lanczos algorithm.
	 *
	 * @param A Black box matrix of the system
	 * @param x Place to store solution vector
	 * @param b Right-hand side
	 * @param F Field over which to perform computations
	 * @param traits \ref SolverTraits  structure with user-specified parameters
	 * @return Reference to solution vector
	 */

	template <class Field, class Vector, class Blackbox>
	Vector &solve (const Blackbox                  &A,
		       Vector                          &x,
		       const Vector                    &b,
		       const Field                     &F,
		       const LanczosTraits        &traits)
	{
		LanczosSolver<Field, Vector> solver (F, traits);
		return solver.solve (A, x, b);
	}

	/** Solve Ax=b over field F using the block Lanczos method.
	 *
	 * This is a general interface for the linear system solving capabilities of
	 * LinBox. If the system is nonsingular, it returns the unique solution, storing
	 * it in the vector x. If the system is consistent and singular, it returns a
	 * random solution. Repeated calls to this function can give a complete
	 * description of the solution manifold. If the system is inconsistent and the
	 * \ref SolverTraits structure has result checking turned on, it runs through
	 * the number of iterations specified in \c traits and throws a
	 * \ref SolveFailed exception if it cannot find a solution.
	 *
	 * This specialization uses the block Lanczos algorithm.
	 *
	 * @param A Black box matrix of the system
	 * @param x Place to store solution vector
	 * @param b Right-hand side
	 * @param F Field over which to perform computations
	 * @param traits \ref SolverTraits structure with user-specified parameters
	 * @return Reference to solution vector
	 */

	template <class Field, class Vector, class Blackbox>
	Vector &solve (const Blackbox                &A,
		       Vector                        &x,
		       const Vector                  &b,
		       const Field                   &F,
		       const BlockLanczosTraits &traits)
	{
		MGBlockLanczosSolver<Field> solver (F, traits);
		solver.solve (A, x, b);
		return x;
	}

	/** Solve Ax=b over field F using Gaussian elimination.
	 *
	 * This is a general interface for the linear system solving capabilities of
	 * LinBox. If the system is nonsingular, it returns the unique solution, storing
	 * it in the vector x. If the system is consistent and singular, it returns a
	 * random solution. Repeated calls to this function can give a complete
	 * description of the solution manifold. If the system is inconsistent and the
	 * \ref SolverTraits  structure supplied requests certification of
	 * inconsistency, it throws an \ref InconsistentSystem exception, which
	 * includes a certificate of inconsistency. Otherwise, it runs through the
	 * number of iterations specified in \p traits and throws a
	 * \ref SolveFailed  exception if it cannot find a solution.
	 *
	 * @param A Black box matrix of the system
	 * @param x Place to store solution vector
	 * @param b Right-hand side
	 * @param F Field over which to perform computations
	 * @param traits \ref SolverTraits  structure with user-specified parameters
	 * @return Reference to solution vector
	 */

	template <class Field, class Matrix, class Vector>
	Vector &solve (const Matrix                     &A,
		       Vector                           &x,
		       const Vector                     &b,
		       const Field                      &F,
		       const BlasEliminationTraits &traits)
	{
		// N.B. This is a place holder; I am intending to fix this very shortly
		throw LinboxError ("Elimination-based solver not implemented");
		return x;
	}

	/** Enumeration for results of next solver.
	 *
	 * SOLVE_SUCCESSFUL - System solution was succesful, \c x holds the solution
	 * vector
	 * SOLVE_INCONSISTENT - System was inconsistent, \c u holds the certificate
	 * of inconsistency and \c x is untouched
	 * SOLVE_FAILED - Neither a system solution nor a certificate of inconsistency
	 * could be obtained before the maximum number of trials elapsed. Both \c x
	 * and \c u are untouched.
	 */

	enum SolveResult {
		SOLVE_SUCCESSFUL, SOLVE_INCONSISTENT, SOLVE_FAILED
	};

	/** Solve Ax=b over field F, returning consistency indicator
	 *
	 * This is a variant of \c solve that does not throw any exceptions unless
	 * the user makes an error. It returns a \ref SolveResult  enum indicating
	 * whether the solve operation was successful, the system was inconsistent, or
	 * the solve operation failed. The certificate of inconsistency, if requested,
	 * is stored in a reference parameter supplied to this variant.
	 *
	 * @param A Black box matrix of the system
	 * @param x Place to store solution vector
	 * @param b Right-hand side
	 * @param u Place to store certificate of inconsistency
	 * @param F Field over which to perform computations
	 * @param traits \ref SolverTraits  structure with user-specified parameters
	 * @return \ref SolveResult indicating whether the solve operation was
	 * successful
	 */

	template <class Field, class Blackbox, class Vector, class MethodTraits>
	SolveResult solve (const Blackbox     &A,
			   Vector             &x,
			   const Vector       &b,
			   const Field        &F,
			   Vector             &u,
			   const MethodTraits &traits = MethodTraits ())
	{
		try {
			solve (A, x, b, F, traits);
		}
		catch (SolveFailed) {
			return SOLVE_FAILED;
		}
		catch (InconsistentSystem<Vector> e) {
			VectorDomain<Field> VD (F);
			F.copy (u, e.certificate ());
			return SOLVE_INCONSISTENT;
		}

		return SOLVE_SUCCESSFUL;
	}
#if 0
	@}
#endif

}

#endif // __LINBOX_bbsolve_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
