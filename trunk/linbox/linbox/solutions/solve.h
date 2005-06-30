/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/solve.h
 *  Evolved from an earlier one by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * See COPYING for license information.
 */

#ifndef __SOLVE_H
#define __SOLVE_H

#include <vector>
#include <algorithm>

// must fix this list...
#include "linbox/algorithms/wiedemann.h"
#include "linbox/algorithms/lanczos.h"
#include "linbox/algorithms/block-lanczos.h"
#include "linbox/blackbox/dense.h"
#include "linbox/util/debug.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/solutions/methods.h"

namespace LinBox 
{

	// for specialization with respect to the DomainCategory
    template< class Vector, class Blackbox, class SolveMethod, class DomainCategory>
    Vector & solve (
		Vector & 				x,
        const Blackbox &        A,
		const Vector &			b,
        const DomainCategory &  tag,
        const SolveMethod &     M
		SolveStatus * 			s = 0);

	/** \brief Solve Ax = b, for x.
	 *
	 * Vector x such that Ax = b is returned.  
	 In the case of a singular matrix A, if the system is consistent, a random
	 solution is returned by default.  The method parameter may contain
	 ??? indicating that an arbitrary element of the solution space may be
	 returned (can be faster).  
	 If the system is inconsistent the zero vector is returned and the 
	 status, if non-null, is set to indicate inconsistency.
	 
         \ingroup solutions
        */
    template< class Vector, class Blackbox, class SolveMethod>
    Vector & solve (
		Vector & 				x,
        const Blackbox &        A,
		const Vector &			b,
        const SolveMethod &     M
		SolveStatus * 			s = 0)
	{ return solve(x, A, b, 
			typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

	// the solve with default method
    template< class Vector, class Blackbox>
	Vector& solve(Vector& x, const Blackbox& A, const Vector& b)
	{ return solve(x, A, b, Method::Hybrid()); }

// in methods.h FoobarMethod and Method::Foobar are the same class.
// in methods.h template<BB> bool useBB(const BB& A) is defined.
//   rowDim > 500 or colDim > 500 might be it.

	// specialize this on blackboxes which have local methods
	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
					const Method::Hybrid& m)
	{	
		if (useBB(A)) return solve(x, A, b, Method::Blackbox(m)); 
		else return solve(x, A, b, Method::Elimination(m));
	}

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
					const Method::Blackbox& m)
	{ 
	// chosen because it is best and/or most reliable currently available choice
		return solve(x, A, b, Method::Wiedemann(m));
	// future:
	//	return solve(x, A, b, Method::BlockLanzos(m));
	}

// temporary
#define inBlasRange(p) true

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
					const Method::Elimination& m)
	{ 
		integer c, p;
		A.field().cardinality(c);
		A.field().characteristic(p);
		if ( p == 0 || (c == p && inBlasRange(p)) )
			return solve(x, A, b, 
					FieldTraits<typename BB::Field>::categoryTag(), 
					Method::BlasElimination(m)); 
  		else 
			return solve(x, A, b, 
					FieldTraits<typename BB::Field>::categoryTag(), 
					Method::NonBlasElimination(m)); 
	}

// BlasElimination section ///////////////////

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
           			const RingCategories::ModularTag tag, 
					const Method::BlasElimination& m)
	{ //Pascal puts in call to base case of dense rational solver (which copies BB to blasmat...)
		return x;
	}

	// specialization when no need to copy matrix
	template <class Vector, class Field> 
	Vector& solve(Vector& x, const DenseMatrix<Field>& A, const Vector& b, 
           			const RingCategories::ModularTag tag, 
					const Method::BlasElimination& m)
	{ //Pascal puts in call to base case of dense rational solver 
		return x;
	}

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
           			const RingCategories::IntegerTag tag, 
					const Method::BlasElimination& m)
	{ 
		DenseMatrix<typename BB::Field> B(A); // copy A into B
		return solve(x, B, b, tag, m);
	} 

	// specialization when no need to copy matrix
	template <class Vector, class Field> 
	Vector& solve(Vector& x, const DenseMatrix<Field>& A, const Vector& b, 
   			        const RingCategories::IntegerTag tag, 
					const Method::BlasElimination& m)
	{ // Pascal puts in call to dense rational solver (which choses prime...)
		return x;
	} 

/*
	struct BlasEliminationCRASpecifier;
	// Extra case put in (1) for timing comparison or (2) for parallelism or 
	// (3) as an example of how we might leave an abandoned choice around in a 
	// callable state for future reference 
	template <class Vector, class Field> 
	Vector& solve(Vector& x, const DenseMatrix<Field>& A, const Vector& b, 
   			        const RingCategories::IntegerTag tag, 
					const BlasEliminationCRASpecifier & m)
	{ // (low priority) J-G puts in code using CRA object CRA and solve(x, A, b, ModularTag, Method::BlasElimination) 
		return x;
	} 
*/

// NonBlasElimination section ////////////////

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
   					const RingCategories::ModularTag tag, 
					const Method::NonBlasElimination& m)
	{	DenseMatrix<typename BB::Field> B(A); // copy
		return solve(x, B, b, tag, m);
	}

	// specialization when no need to copy
	template <class Vector, class Field> 
	Vector& solve(Vector& x, const DenseMatrix<Field>& A, const Vector& b, 
   			        const RingCategories::ModularTag tag, 
					const Method::NonBlasElimination& m)
	{ //Dave finds a call in original solve.h maybe.
		return x;
	}

// note: no need for NonBlasElimination when RingCategory is integer

// WiedemannElimination section ////////////////

	// may throw SolverFailed or InconsistentSystem
	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
   			        const RingCategories::ModularTag tag, 
					const Method::Wiedemann& m)
	{
		// adapt to earlier signature of wiedemann solver
		return solve (A, x, b, A.field(), m);
	}

	template <class Vector, class BB> 
	Vector& solve(Vector& x, const BB& A, const Vector& b, 
   			        const RingCategories::IntegerTag tag, 
					const Method::Wiedemann& m)
	{ 	// Must put in cra loop
		/*struct solver 
		{ Vector& operator ()(Vector& x, const Modular<double>& F) { 
			// make modular Am bm from A, b, make sm, then
		MatrixHom::mod(
		solve (xm, Am, bm, m, sm)
		return xm;
		}
		*/
		return x;
	}

/* remark 1.  I used copy constructors when switching method types.
But if the method types are (empty) child classes of a common  parent class containing
all the information, then casts can be used in place of copies.

remark 2. Stopped here.  It remains to put some of the below stuff in algorithms
and call it in the appropriate places above.

*/ 
#if 0
/** @name Solvers
 * @memo Solving linear system Ax = b over the field F.
 */
//@{
/** Solve Ax=b over field F using Wiedemann's method, with inconsistency certificate.
 * 
 * This is a general interface for the linear system solving
 * capabilities of LinBox. If the system is nonsingular, it returns
 * the unique solution, storing it in the vector x. If the system is
 * consistent and singular, it returns a random solution. Repeated
 * calls to this function can give a complete description of the
 * solution manifold. If the system is inconsistent and the
 * \Ref{SolverTraits} structure supplied requests certification of
 * inconsistency, it fills in the certificate of
 * inconsistency. Otherwise, it runs through the number of iterations
 * specified in @code{traits} and throws a \Ref{SolveFailed} exception
 * if it cannot find a solution.
 *
 * This specialization uses Wiedemann's algorithm and is the default.
 *
 * @param A Black box matrix of the system
 * @param x Place to store solution vector
 * @param b Right-hand side
 * @param u Vector in which to store certificate of inconsistency, if required
 * @param F Field over which to perform computations
 * @param traits \Ref{SolverTraits} structure with user-specified parameters
 * @return Reference to solution vector
 */

template <class Field, class Vector>
typename WiedemannSolver<Field, Vector>::ReturnStatus 
solve (const BlackboxArchetype&A,
       Vector                          &x,		       
       const Vector                    &b,
       Vector                          &u,
       const Field                     &F,
       const WiedemannTraits &traits = WiedemannTraits ())
{
	WiedemannSolver<Field, Vector> solver (F, traits);
	return solver.solve (A, x, b, u);
}

/** Solve Ax=b over field F using the Wiedemann method.
 *
 * This version differs from the one above in that there is no extra
 * parameter for the certificate of inconsistency, and it throws
 * exceptions if the solution fails. It also returns a reference to
 * the solution vector.
 */

template <class Field, class Vector>
Vector &solve (const BlackboxArchetype&A,
	       Vector                          &x,		       
	       const Vector                    &b,
	       const Field                     &F,
	       const WiedemannTraits &traits = WiedemannTraits ())
{
	Vector u;
	WiedemannSolver<Field, Vector> solver (F, traits);

	VectorWrapper::ensureDim (u, A.rowdim ());

	switch (solver.solve (A, x, b, u)) {
	    case WiedemannSolver<Field, Vector>::OK:
		return x;

	    case WiedemannSolver<Field, Vector>::FAILED:
		throw SolveFailed ();

	    case WiedemannSolver<Field, Vector>::SINGULAR:
		throw SolveFailed ();

	    case WiedemannSolver<Field, Vector>::INCONSISTENT:
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
 * \Ref{SolverTraits} structure has result checking turned on, it runs through
 * the number of iterations specified in @code{traits} and throws a
 * \Ref{SolveFailed} exception if it cannot find a solution.
 *
 * This specialization uses the Lanczos algorithm.
 *
 * @param A Black box matrix of the system
 * @param x Place to store solution vector
 * @param b Right-hand side
 * @param F Field over which to perform computations
 * @param traits \Ref{SolverTraits} structure with user-specified parameters
 * @return Reference to solution vector
 */

template <class Field, class Vector>
Vector &solve (const BlackboxArchetype&A,
	       Vector                          &x,		       
	       const Vector                    &b,
	       const Field                     &F,
	       const LanczosTraits &traits)
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
 * \Ref{SolverTraits} structure has result checking turned on, it runs through
 * the number of iterations specified in @code{traits} and throws a
 * \Ref{SolveFailed} exception if it cannot find a solution.
 *
 * This specialization uses the block Lanczos algorithm.
 *
 * @param A Black box matrix of the system
 * @param x Place to store solution vector
 * @param b Right-hand side
 * @param F Field over which to perform computations
 * @param traits \Ref{SolverTraits} structure with user-specified parameters
 * @return Reference to solution vector
 */

template <class Field, class Vector, class Blackbox>
Vector &solve (const Blackbox &A,
	       Vector         &x,		       
	       const Vector   &b,
	       const Field    &F,
	       const BlockLanczosTraits &traits)
{
	BlockLanczosSolver<Field> solver (F, traits);
	return solver.solve (A, x, b);
}

/** Solve Ax=b over field F using Gaussian elimination.
 * 
 * This is a general interface for the linear system solving capabilities of
 * LinBox. If the system is nonsingular, it returns the unique solution, storing
 * it in the vector x. If the system is consistent and singular, it returns a
 * random solution. Repeated calls to this function can give a complete
 * description of the solution manifold. If the system is inconsistent and the
 * \Ref{SolverTraits} structure supplied requests certification of
 * inconsistency, it throws an \Ref{InconsistentSystem} exception, which
 * includes a certificate of inconsistency. Otherwise, it runs through the
 * number of iterations specified in @code{traits} and throws a
 * \Ref{SolveFailed} exception if it cannot find a solution.
 *
 * @param A Black box matrix of the system
 * @param x Place to store solution vector
 * @param b Right-hand side
 * @param F Field over which to perform computations
 * @param traits \Ref{SolverTraits} structure with user-specified parameters
 * @return Reference to solution vector
 */

template <class Field, class Matrix, class Vector>
Vector &solve (const Matrix &A,
	       Vector       &x,		       
	       const Vector &b,
	       const Field  &F,
	       const BlasEliminationTraits &traits)
{
	// N.B. This is a place holder; I am intending to fix this very shortly
	throw LinboxError ("Elimination-based solver not implemented");
	return x;
}

/** Enumeration for results of next solver.
 *
 * SOLVE_SUCCESSFUL - System solution was succesful, @code{x} holds the solution
 * vector 
 * SOLVE_INCONSISTENT - System was inconsistent, @code{u} holds the certificate
 * of inconsistency and @code{x} is untouched
 * SOLVE_FAILED - Neither a system solution nor a certificate of inconsistency
 * could be obtained before the maximum number of trials elapsed. Both @code{x}
 * and @code{u} are untouched.
 */

enum SolveResult {
	SOLVE_SUCCESSFUL, SOLVE_INCONSISTENT, SOLVE_FAILED
};

/** Solve Ax=b over field F, returning consistency indicator
 *
 * This is a variant of @code{solve} that does not throw any exceptions unless
 * the user makes an error. It returns a \Ref{SolveResult} enum indicating
 * whether the solve operation was successful, the system was inconsistent, or
 * the solve operation failed. The certificate of inconsistency, if requested,
 * is stored in a reference parameter supplied to this variant.
 *
 * @param A Black box matrix of the system
 * @param x Place to store solution vector
 * @param b Right-hand side
 * @param u Place to store certificate of inconsistency
 * @param F Field over which to perform computations
 * @param traits \Ref{SolverTraits} structure with user-specified parameters
 * @return \Ref{SolveResult} indicating whether the solve operation was
 * successful
 */

template <class Field, class Blackbox, class Vector, class MethodTraits>
SolveResult solve (const Blackbox     &A,
		   Vector             &x,		       
		   const Vector       &b,
		   const Field        &F,
		   Vector             &u,
		   const SolverTraits<MethodTraits> &traits = SolverTraits<MethodTraits> ())
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
//@}
#endif

}

#endif // __SOLVE_H
