/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/wiedemann.h
 * Copyright (C) 2002 Zhendong Wan
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-10-02  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Refactoring:
 * Put everything inside a WiedemannSolver class, with the following
 * interface:
 *    solveNonsingular - Solve a nonsingular system
 *    solveSingular - Solve a general singular system
 *    findRandomSolution - Find a random solution to a singular preconditioned
 *                         problem
 *    findNullSpaceElement - Find an element of the right nullspace
 *    certifyInconsistency - Find a certificate of inconsistency for a
 *                           linear system
 *    precondition - Form a preconditioner and return it
 * ------------------------------------
 * 2002-08-09  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Move the Wiedemann stuff over to this file
 *
 * Create a singular and nonsingular version that is a bit intelligent about
 * which one to use in different circumstances
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __WIEDEMANN_H
#define __WIEDEMANN_H

#include <vector>
#include <algorithm>

#include "linbox/blackbox/archetype.h"
#include "linbox/util/debug.h"
#include "linbox/field/vector-domain.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{

/** Wiedemann system solver class
 *
 * This class encapsulates all of the functionality for linear system
 * solving with Wiedemann's algorithm. It includes the random solution and
 * random nullspace element of Kaltofen and Saunders (1991), as well as the
 * certificate of inconsistency of Giesbrecht, Lobo, and Saunders (1998).
 */

template <class Field, class Vector>
class WiedemannSolver 
{
    public:

	/** Constructor
	 * @param F Field over which to operate
	 * @param traits @ref{SolverTraits} structure describing user
	 *               options for the solver 
	 */
	WiedemannSolver (const Field &F, const SolverTraits &traits)
		: _traits (traits), _F (F), _randiter (F), _VD (F)
	{}

	/** Constructor with a random iterator
	 * @param F Field over which to operate
	 * @param traits @ref{SolverTraits} structure describing user
	 *               options for the solver 
	 * @param r Random iterator to use for randomization
	 */
	WiedemannSolver (const Field &F, const SolverTraits &traits, typename Field::RandIter r)
		: _traits (traits), _F (F), _randiter (r), _VD (F)
	{}

	// @name Solvers

	//@{

	/** Solve a system Ax=b, giving a random solution if the system is
	 * singular and consistent, and a certificate of inconsistency (if
	 * requested) otherwise.
	 *
	 * @param A Black box of linear system
	 * @param x Vector in which to store solution
	 * @param b Right-hand side of system
	 * @return Reference to solution vector
	 */
	Vector &solve (const BlackboxArchetype<Vector> &A, Vector &x, const Vector &b);

	/** Solve a nonsingular system Ax=b. Throws @ref{SingularSystem} if the
	 * system turns out to be singular, and @ref{SolveFailed} if the system
	 * solution computed is not actually a solution.
	 *
	 * This is a "Las Vegas" method, which makes use of randomization. It
	 * attempts to certify that the system solution is correct. It will only
	 * make one attempt to solve the system before giving up.
	 *
	 * @param A Black box of linear system
	 * @param x Vector in which to store solution
	 * @param b Right-hand side of system
	 * @param useRandIter true if solveNonsingular should use a random
	 *                    iterator for the Krylov sequence computation;
	 *                    false if it should use the right-hand side
	 * @return Reference to solution vector
	 */
	Vector &solveNonsingular (const BlackboxArchetype<Vector> &A, Vector &x, const Vector &b, bool useRandIter = false);

	/** Solve a general singular linear system. Throws
	 * @ref{InconsistentSystem} if the system turns out to be inconsistent,
	 * producing a certificate of inconsistency if requested. Throws
	 * BadPreconditioner if the preconditioner formed to give A a generic
	 * rank profile does not work.
	 *
	 * @param A Black box of linear system
	 * @param x Vector in which to store solution
	 * @param b Right-hand side of system
	 * @param r Rank of A
	 * @return Reference to solution vector
	 */
	Vector &solveSingular (const BlackboxArchetype<Vector> &A, Vector &x, const Vector &b, unsigned long r);

	/** Get a random solution to a singular system Ax=b of rank r with
	 * generic rank profile. Throws @ref{BadPreconditioner} if the black box
	 * A does not have generic rank profile.
	 *
	 * @param A Black box of linear system
	 * @param x Vector in which to store solution
	 * @param b Right-hand side of system
	 * @param r Rank of A
	 * @param P Left preconditioner (NULL if none needed)
	 * @param Q Right preconditioner (NULL if none needed)
	 * @return Reference to solution vector
	 */
	Vector &findRandomSolution (const BlackboxArchetype<Vector> &A,
				    Vector                          &x,
				    const Vector                    &b,
				    size_t                           r,
				    const BlackboxArchetype<Vector> *P = NULL,
				    const BlackboxArchetype<Vector> *Q = NULL);

	/** Get a random element of the right nullspace of A.
	 *
	 * @param x Vector in which to store nullspace element
	 * @param A Black box of which to find nullspace element
	 * @param r Rank of A
	 * @param P Left preconditioner, if applicable
	 * @param Q Right preconditioner, if applicable
	 */
	Vector &findNullspaceElement (Vector                          &x,
				      const BlackboxArchetype<Vector> &A);

	/** Get a certificate u that the given system Ax=b is
	 * inconsistent, if one can be found.
	 *
	 * @param u Vector in which to store certificate
	 * @param A Black box for the linear system
	 * @param b Right-hand side for the linear system
	 * @param r Rank of A
	 * @param P Left preconditioner, if applicable
	 * @return true if a certificate can be found in one iteration; u
	 *         is filled in with that certificate; and false otherwise
	 */
	bool certifyInconsistency (Vector                          &u,
				   const BlackboxArchetype<Vector> &A,
				   const Vector                    &b);

	//@}

	// @name Preconditioners

	//@{

	/** Given a blackbox archetype A, construct preconditioners P and Q of
	 * the type requested in the constructor and return the preconditioned
	 * matrix PAQ. Returns a pointer to A and leaves P and Q untouched if no
	 * preconditioning was requested.
	 * @param A Black box to precondition
	 * @param PAQ Pointer to preconditioned matrix
	 * @param P Pointer to left preconditioner
	 * @param Q Pointer to right preconditioner
	 */
	const BlackboxArchetype<Vector> *precondition (const BlackboxArchetype<Vector>  &A,
						       BlackboxArchetype<Vector>       *&PAQ,
						       BlackboxArchetype<Vector>       *&P,
						       BlackboxArchetype<Vector>       *&Q);

	//@}

	// @name Exceptions

	//@{

	/** Exception thrown when the computed solution vector is not a true
	 * solution to the system, but none of the problems cited below exist
	 */
	class SolveFailed {};

	/** Exception thrown when a singular system is passed to
	 * @ref{solveNonsingular}
	 */
	class SingularSystem {};

	/** Exception thrown when the system to be solved is
	 * inconsistent. Contains a certificate of inconsistency.
	 */
	class InconsistentSystem 
	{
	    public:
		InconsistentSystem (Vector &u)
			: _cert (true), _u (u)
		{}

		InconsistentSystem ()
			: _cert (false)
		{}

		const Vector &u () const { return _u; }
		bool certified () const { return _cert; }

	    private:

		bool _cert;
		Vector _u;
	};

	/** Exception thrown when the system is not properly conditioned,
	 * i.e. does not have generic rank profile
	 */
	class BadPreconditioner {};

	//@}

    private:

	// Make an m x m lambda-sparse matrix, c.f. Mulders (2000)
	SparseMatrix0<Field, Vector> *makeLambdaSparseMatrix (size_t m);

	const SolverTraits       &_traits;
	const Field              &_F;
	typename Field::RandIter  _randiter;
	VectorDomain<Field>       _VD;
};

}

#include "linbox/algorithms/wiedemann.inl"

#endif // __WIEDEMANN_H
