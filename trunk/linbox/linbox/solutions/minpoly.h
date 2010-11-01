/* linbox/solutions/minpoly.h
 * Copyright (C) 1999, 2001 Jean-Guillaume Dumas
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
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

#ifndef __LINBOX_minpoly_H
#define __LINBOX_minpoly_H

#include "linbox/solutions/methods.h"
#include "linbox/util/commentator.h"

#include "linbox/blackbox/squarize.h"

#include "linbox/algorithms/blas-domain.h"
#include "linbox/algorithms/wiedemann.h"

#ifdef __LINBOX_HAVE_MPI
#include "linbox/util/mpicpp.h"
#include "linbox/algorithms/cra-mpi.h"
#endif

#include "linbox/algorithms/minpoly-integer.h"

namespace LinBox 
{
	
	/*- @brief Minimal polynomial of a blackbox linear operator A.
	 * The resulting polynomial is a vector of coefficients.
	 * Somewhere we should document our handling of polys.
	 */
	template < class Blackbox, class Polynomial, class DomainCategory, class MyMethod>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     const DomainCategory& tag,
			     const MyMethod& M);
/*
	template < class Blackbox, class Polynomial, class MyMethod>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     const RingCategories::RationalTag& tag,
			     const MyMethod& M)
	{
		throw LinboxError("LinBox ERROR: minpoly is not yet defined over a rational domain");
	}
*/

        /** \brief  ...using an optional Method parameter
	    \param P - the output minimal polynomial.  If the polynomial is
	    of degree d, this random access container has size d+1, the 0-th entry is 
	    the constant coefficient and the d-th is 1 since the minpoly is monic.
	    \param A - a blackbox matrix
	    Optional \param M - the method object.  Generally, the default
	    object suffices and the algorithm used is determined by the class of M.
	    Basic methods are Method::Blackbox, Method::Elimination, and Method::Hybrid
	    (the default).
	    See methods.h for more options.
	    \return a reference to P.
	*/
	template < class Blackbox, class Polynomial, class MyMethod>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     const MyMethod& M){
		return minpoly (P, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

        /// \brief  ...using default Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (Polynomial &P, 
			     const Blackbox &A)    
	{        return minpoly (P, A, Method::Hybrid());    }



	// The minpoly with Hybrid Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::Hybrid& M)
	{
		// not yet a hybrid
		return minpoly(P, A, tag, Method::Blackbox(M));
	}

	// The minpoly with Hybrid Method on DenseMatrix
	template<class Polynomial, class Field>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const DenseMatrix<Field> 			&A,
			     const RingCategories::ModularTag          &tag,
			     const Method::Hybrid& M)
	{
		return minpoly(P, A, tag, Method::Elimination(M));
	}
	
	// The minpoly with Hybrid Method on BlasBlackbox
	template<class Polynomial, class Field>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const BlasBlackbox<Field> 			&A,
			     const RingCategories::ModularTag          &tag,
			     const Method::Hybrid& M)
	{
		return minpoly(P, A, tag, Method::Elimination(M));
	}

	// The minpoly with Elimination Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::Elimination& M)
	{
		return minpoly(P, A, tag, Method::BlasElimination(M));
	}

	// The minpoly with BlasElimination Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::BlasElimination& M)
	{
            commentator.start ("Convertion to BLAS Minimal polynomial", "blasconvert");
            
            if (A.coldim() != A.rowdim()) {
                commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "Squarize matrix" << std::endl;
                Squarize<Blackbox> B(&A);              
                BlasBlackbox< typename Blackbox::Field > BBB (B);
                BlasMatrixDomain< typename Blackbox::Field > BMD (BBB.field());
                commentator.stop ("done", NULL, "blasconvert");
       
                return BMD.minpoly (P, static_cast<const BlasMatrix<typename Blackbox::Field::Element>& >(BBB));
            } else {
                BlasBlackbox< typename Blackbox::Field > BBB (A);
                BlasMatrixDomain< typename Blackbox::Field > BMD (BBB.field());
                commentator.stop ("done", NULL, "blasconvert");
                return BMD.minpoly (P, static_cast<const BlasMatrix<typename Blackbox::Field::Element>& >(BBB));
            }
	}
        
	// The minpoly with BlackBox Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::Blackbox& M)
	{
            if (M.certificate()) {
                    // Will make a word size extension 
                    // when field size is too small
                return minpoly(P, A, tag, Method::ExtensionWiedemann (M));
            } else
		return minpoly(P, A, tag, Method::Wiedemann (M));
	}


        
}

// ---------------------------------------------------------
// Chinese Remaindering generic wrappers for integer minpoly
// ---------------------------------------------------------

#include "linbox/field/modular-double.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"

#include "linbox/algorithms/rational-cra2.h"
#include "linbox/algorithms/varprec-cra-early-multip.h"
#include "linbox/algorithms/minpoly-rational.h"

namespace LinBox {
   
	template <class Blackbox, class MyMethod>
	struct IntegerModularMinpoly {       
		const Blackbox &A;
		const MyMethod &M;

		IntegerModularMinpoly(const Blackbox& b, const MyMethod& n) 
			: A(b), M(n) {}
        
        
		template<typename Polynomial, typename Field>
		Polynomial& operator()(Polynomial& P, const Field& F) const {
                    typedef typename Blackbox::template rebind<Field>::other FBlackbox;
                    FBlackbox Ap(A, F);
                    return minpoly( P, Ap, typename FieldTraits<Field>::categoryTag(), M);
		}            
	};

	template <class Polynomial, class Blackbox, class MyMethod>
	Polynomial &minpoly (Polynomial 			&P, 
                             const Blackbox                     &A,
                             const RingCategories::IntegerTag   &tag,
                             const MyMethod                     &M)
	{
//             if (A.coldim() != A.rowdim())
//                     throw LinboxError("LinBox ERROR: matrix must be square for minimal polynomial computation\n");

#ifdef __LINBOX_HAVE_MPI
		Communicator *c = M.communicatorp();
		if(!c || c->rank() == 0) 
			commentator.start ("Integer Minpoly", "Iminpoly");
		else{
			//commentator.setMaxDepth(0);
			//commentator.setMaxDetailLevel(0);
			//commentator.setPrintParameters(0, 0, 0);
		}
#else
		commentator.start ("Integer Minpoly", "Iminpoly");
#endif
		// 0.7213475205 is an upper approximation of 1/(2log(2))
		RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 
		IntegerModularMinpoly<Blackbox,MyMethod> iteration(A, M);
		std::vector<integer> PP; // use of integer du to non genericity of cra. PG 2005-08-04
#ifdef __LINBOX_HAVE_MPI
		MPIChineseRemainder< EarlyMultipCRA<Modular<double> > > cra(3UL, c);
		cra(PP, iteration, genprime);
#else
		ChineseRemainder< EarlyMultipCRA<Modular<double> > > cra(3UL);
		cra(PP, iteration, genprime);
#endif
		size_t i =0;
		P.resize(PP.size());
		for (typename Polynomial::iterator it= P.begin(); it != P.end(); ++it, ++i)
			A.field().init(*it, PP[i]);

#ifdef __LINBOX_HAVE_MPI
		if(c || c->rank() == 0) 
#endif
			commentator.stop ("done", NULL, "Iminpoly");
		return P;
	}

	template < class Blackbox, class Polynomial, class MyMethod>
        Polynomial &minpoly (Polynomial& P, const Blackbox& A,
                               const RingCategories::RationalTag& tag, const MyMethod& M)
        {
	        commentator.start ("Rational Minpoly", "Rminpoly");

		RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205));
		RationalRemainder2< VarPrecEarlyMultipCRA<Modular<double> > > rra(3UL);
		IntegerModularMinpoly<Blackbox,MyMethod> iteration(A, M);

		std::vector<Integer> PP; // use of integer due to non genericity of cra. PG 2005-08-04
		Integer den;
		rra(PP,den, iteration, genprime);
		size_t i =0;
		P.resize(PP.size());
		for (typename Polynomial::iterator it= P.begin(); it != P.end(); ++it, ++i)
			A.field().init(*it, PP[i],den);

		commentator.stop ("done", NULL, "Rminpoly");

		return P;
	}

	template < class Field, template<class> class Polynomial, class MyMethod>
	Polynomial<typename Field::Element> &minpoly (Polynomial<typename Field::Element>& P, const DenseMatrix<Field>& A,
				const RingCategories::RationalTag& tag, const MyMethod& M)
	{
		commentator.start ("Dense Rational Minpoly", "Rminpoly");

		rational_minpoly(P,A,M);
		return P;
	}
} // end of LinBox namespace
#endif // __LINBOX_minpoly_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
