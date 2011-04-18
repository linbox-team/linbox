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

#ifndef __LINBOX_wiedemann_H
#define __LINBOX_wiedemann_H

#include <vector>
#include <algorithm>

#include "linbox/blackbox/archetype.h"
#include "linbox/blackbox/squarize.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/util/debug.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/solutions/methods.h"

#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"

// massey recurring sequence solver
#include "linbox/algorithms/massey-domain.h"     

namespace LinBox 
{


	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     RingCategories::ModularTag tag,
			     const Method::Wiedemann& M = Method::Wiedemann ())
	{
		typedef typename Blackbox::Field Field;
		typename Field::RandIter i (A.field());
		unsigned long            deg;

		commentator.start ("Wiedemann Minimal polynomial", "minpoly");

                if (A.coldim() != A.rowdim()) {
                    commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "Virtually squarize matrix" << std::endl;
                    
                    Squarize<Blackbox> B(&A);
                    BlackboxContainer<Field, Squarize<Blackbox> > TF (&B, A.field(), i);
                    MasseyDomain< Field, BlackboxContainer<Field, Squarize<Blackbox> > > WD (&TF, M.earlyTermThreshold ());
                    
                    WD.minpoly (P, deg);                    
                } else if (M.symmetric ()) {
                    typedef BlackboxContainerSymmetric<Field, Blackbox> BBContainerSym;
                    BBContainerSym TF (&A, A.field(), i);
                    MasseyDomain< Field, BBContainerSym > WD (&TF, M.earlyTermThreshold ());
                    
                    WD.minpoly (P, deg);
		} else {
                    typedef BlackboxContainer<Field, Blackbox> BBContainer;
                    BBContainer TF (&A, A.field(), i);
                    MasseyDomain< Field, BBContainer > WD (&TF, M.earlyTermThreshold ());
                    
                    WD.minpoly (P, deg);
#ifdef INCLUDE_TIMING
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for applies:      " << TF.applyTime () << std::endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for dot products: " << TF.dotTime () << std::endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for discrepency:  " << WD.discrepencyTime () << std::endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for LSR fix:      " << WD.fixTime () << std::endl;
#endif // INCLUDE_TIMING
                }



		commentator.stop ("done", NULL, "minpoly");

		return P;
	}
}

#ifdef __LINBOX_HAVE_GIVARO
#ifndef LINBOX_EXTENSION_DEGREE_MAX
#define LINBOX_EXTENSION_DEGREE_MAX 19
#endif

#include "linbox/blackbox/sparse.h"
#include "linbox/field/modular.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/field/givaro-extension.h"
#include "linbox/field/map.h"

namespace LinBox 
{  
	// The minpoly with BlackBox Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::ExtensionWiedemann& M)
	{
            typedef typename Blackbox::Field Field;
            const Field& F = A.field();
            integer a,c; F.cardinality(a); F.characteristic(c);
            if (a != c) {
                unsigned long extend = (unsigned long)FF_EXPONENT_MAX(a,(integer)LINBOX_EXTENSION_DEGREE_MAX);
                if (extend > 1) {
                    commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << "Extension of degree " << extend << std::endl;
                    GivaroExtension<Field> EF( F, extend);
                    
                    typedef typename Blackbox::template rebind< GivaroExtension<Field>  >::other FBlackbox;
                    
                    FBlackbox Ap(A, EF);
                    
                    std::vector< typename GivaroExtension<Field>::Element > eP;
                    minpoly(eP, Ap, tag, Method::Wiedemann(M));
                    
                    return PreMap<Field, GivaroExtension<Field> >(F,EF)(P, eP);
                } else
                    return minpoly(P, A, tag, Method::Wiedemann(M)); 
                
            } else {
                unsigned long extend = (unsigned long)FF_EXPONENT_MAX(c,(integer)LINBOX_EXTENSION_DEGREE_MAX);
                if (extend > 1) {
                    commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << "Word size extension : " << extend << std::endl;
                    GivaroGfq EF( (unsigned long)c, extend);                    
                    typedef typename Blackbox::template rebind< GivaroGfq >::other FBlackbox;
                    FBlackbox Ap(A, EF);                    
                    std::vector< typename GivaroGfq::Element > eP;
                    minpoly(eP, Ap, tag, Method::Wiedemann(M));
                    return PreMap<Field, GivaroGfq >(F,EF)(P, eP);
                    
                } else
                    return minpoly(P, A, tag, Method::Wiedemann(M)); 
            }
        }
}
#else
namespace LinBox 
{
	// The minpoly with BlackBox Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::ExtensionWiedemann& M)
	{
            commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << " WARNING, no extension available, returning only a factor of the minpoly\n";
            return minpoly(P, A, tag, Method::Wiedemann (M));
	}
}
#endif

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

#endif //  __LINBOX_wiedemann_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
