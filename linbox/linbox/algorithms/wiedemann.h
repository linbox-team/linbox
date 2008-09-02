/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/wiedemann.h
 * Copyright (C) 2002 Zhendong Wan
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Zhendong Wan <wan@mail.eecis.udel.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ----------------------------------------------------
 * 2003-02-05  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Ripped out all the exception code. Exceptions decided one day to
 * just stop working on my compiler, and they were controversal
 * anyway. Now all the solve functions return a status. There are most
 * likely still bugs in this code, though.
 * ----------------------------------------------------
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
#include "linbox/blackbox/sparse.h"
#include "linbox/util/debug.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{

/** \brief Linear system solvers based on Wiedemann's method.
 * 
 * This class encapsulates all of the functionality for linear system
 * solving with Wiedemann's algorithm. It includes the random solution and
 * random nullspace element of Kaltofen and Saunders (1991), as well as the
 * certificate of inconsistency of Giesbrecht, Lobo, and Saunders (1998).
 */

template <class Field>
class WiedemannSolver 
{
    public:

	/// { OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER }
	enum ReturnStatus {
		OK, FAILED, SINGULAR, INCONSISTENT, BAD_PRECONDITIONER
	};

	/** Constructor
	 * @param F Field over which to operate
	 * @param traits @ref{SolverTraits} structure describing user
	 *               options for the solver 
	 */
	WiedemannSolver (const Field &F, const WiedemannTraits &traits)
		: _traits (traits), _F (F), _randiter (F), _VD (F)
	{}

	/** Constructor with a random iterator
	 * @param F Field over which to operate
	 * @param traits @ref{SolverTraits} structure describing user
	 *               options for the solver 
	 * @param r Random iterator to use for randomization
	 */
	WiedemannSolver (const Field &F,
			 const WiedemannTraits &traits,
			 typename Field::RandIter r)
		: _traits (traits), _F (F), _randiter (r), _VD (F)
	{}

	// @name Solvers
	// try to make the idea work doxy
	/// \defgroup Solvers

	//@{

	/** Solve a system Ax=b, giving a random solution if the system is
	 * singular and consistent, and a certificate of inconsistency (if
	 * specified in traits parameter at construction time) otherwise.
	 *
	 * @param A Black box of linear system
	 * @param x Vector in which to store solution
	 * @param b Right-hand side of system
	 * @param u Vector in which to store certificate of inconsistency
	 * @return Reference to solution vector
	 */
	template<class Blackbox, class Vector>	
	ReturnStatus solve (const Blackbox&A, Vector &x, const Vector &b, Vector &u);

	/** Solve a nonsingular system Ax=b.
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
	template<class Blackbox, class Vector>
	ReturnStatus solveNonsingular (const Blackbox&A,
				       Vector &x,
				       const Vector &b,
				       bool useRandIter = false);

	/** Solve a general singular linear system.
	 *
	 * @param A Black box of linear system
	 * @param x Vector in which to store solution
	 * @param b Right-hand side of system
	 * @param u Vector into which certificate of inconsistency will be stored
	 * @param r Rank of A
	 * @return Return status
	 */
	template<class Blackbox, class Vector>
	ReturnStatus solveSingular (const Blackbox&A,
				    Vector &x,
				    const Vector &b,
				    Vector &u,
				    unsigned long r);

	/** Get a random solution to a singular system Ax=b of rank r with
	 * generic rank profile.
	 *
	 * @param A Black box of linear system
	 * @param x Vector in which to store solution
	 * @param b Right-hand side of system
	 * @param r Rank of A
	 * @param P Left preconditioner (NULL if none needed)
	 * @param Q Right preconditioner (NULL if none needed)
	 * @return Return status
	 */
	template<class Blackbox, class Vector, class Prec1, class Prec2>
	ReturnStatus findRandomSolution (const Blackbox          &A,
					 Vector                  &x,
					 const Vector            &b,
					 size_t                   r,
					 const Prec1             *P,
					 const Prec2             *Q);

	/** Get a random element of the right nullspace of A.
	 *
	 * @param x Vector in which to store nullspace element
	 * @param A Black box of which to find nullspace element
	 */
	template<class Blackbox, class Vector>
	ReturnStatus findNullspaceElement (Vector                &x,
					   const Blackbox        &A);

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
	template<class Blackbox, class Vector>
	bool certifyInconsistency (Vector                          &u,
				   const Blackbox                  &A,
				   const Vector                    &b);

	//@}


    private:

	// Make an m x m lambda-sparse matrix, c.f. Mulders (2000)
	SparseMatrix<Field> *makeLambdaSparseMatrix (size_t m);

	WiedemannTraits _traits;
	const Field                         &_F;
	typename Field::RandIter             _randiter;
	VectorDomain<Field>                  _VD;
};

}

#include "linbox/algorithms/wiedemann.inl"

#endif // __WIEDEMANN_H
