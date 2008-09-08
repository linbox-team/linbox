/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#ifndef __MINPOLY_H
#define __MINPOLY_H

#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/blackbox/squarize.h"
#include "linbox/algorithms/massey-domain.h"     // massey recurring sequence solver
#include "linbox/algorithms/blas-domain.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/commentator.h"
#ifdef __LINBOX_HAVE_MPI
#include "linbox/util/mpicpp.h"
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

	//error handler for rational domain
	template < class Blackbox, class Polynomial, class MyMethod>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     const RingCategories::RationalTag& tag,
			     const MyMethod& M)
	{
		throw LinboxError("LinBox ERROR: minpoly is not yet defined over a rational domain");
	}
 

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
        
}

#ifdef __LINBOX_HAVE_GIVARO
#define LINBOX_EXTENSION_DEGREE_MAX 20

#include "linbox/blackbox/sparse.h"
#include "linbox/field/modular.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/field/givaro-extension.h"
#include "linbox/field/map.h"

namespace LinBox {  
	// The minpoly with BlackBox Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::Blackbox& M)
	{
            if (M.certificate()) {
                typedef typename Blackbox::Field Field;
                const Field& F = A.field();
                integer a,c; F.cardinality(a); F.characteristic(c);
                if (a != c) {
                    unsigned long extend = (unsigned long)FF_EXPONENT_MAX(a,(integer)LINBOX_EXTENSION_DEGREE_MAX);
                    if (extend > 1) {
                        commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << "Extension of degree " << extend << std::endl;
                        GivaroExtension<Field> EF( F, extend);
                        typedef typename Blackbox::template rebind< GivaroExtension<Field>  >::other FBlackbox;
                        FBlackbox * Ap;
                        MatrixHom::map(Ap, A, EF );
                        std::vector< typename GivaroExtension<Field>::Element > eP;

                        minpoly(eP, *Ap, tag, Method::Wiedemann(M));

                        return PreMap<Field, GivaroExtension<Field> >(F,EF)(P, eP);
                    } else
                        return minpoly(P, A, tag, Method::Wiedemann(M)); 
                    
                } else {
                    unsigned long extend = (unsigned long)FF_EXPONENT_MAX(c,(integer)LINBOX_EXTENSION_DEGREE_MAX);
                    if (extend > 1) {
                        commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << "Word size extension : " << extend << std::endl;
                        GivaroGfq EF( (unsigned long)c, extend);                    
                        typedef typename Blackbox::template rebind< GivaroGfq >::other FBlackbox;
                        FBlackbox * Ap;
                        MatrixHom::map(Ap, A, EF );
                        std::vector< typename GivaroGfq::Element > eP;
                        minpoly(eP, *Ap, tag, Method::Wiedemann(M));

                        return PreMap<Field, GivaroGfq >(F,EF)(P, eP);
                        
                    } else
                        return minpoly(P, A, tag, Method::Wiedemann(M)); 
                }
            } else
		return minpoly(P, A, tag, Method::Wiedemann (M));
	}
}
#else
namespace LinBox {
	// The minpoly with BlackBox Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::Blackbox& M)
	{
            commentator.report (Commentator::LEVEL_ALWAYS,INTERNAL_WARNING) << " WARNING, no extension available, returning only a factor of the minpoly\n";
            return minpoly(P, A, tag, Method::Wiedemann (M));
	}
}
#endif
namespace LinBox {
	/*
	  template<class Polynomial, class Blackbox>
	  Polynomial &minpoly (Polynomial& P,
	  const Blackbox& A,
	  const RingCategories::ModularTag          &tag,
	  const Method::Wiedemann& M = Method::Wiedemann ());
	
	  return minpoly (P, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	  }

	  template<class Polynomial, class Blackbox>
	  Polynomial &minpoly (Polynomial& P,
	  const Blackbox& A,
	  RingCategories::IntegerTag tag,
	  const Method::Wiedemann& M = Method::Wiedemann ())
	  {	
	  typedef Modular<double> ModularField;
	  MinPoly<typename Blackbox::Field::Element, ModularField>::minPoly(P, A);

	  return P;
	  }
	*/

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
                }


#ifdef INCLUDE_TIMING
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for applies:      " << TF.applyTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for dot products: " << TF.dotTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for discrepency:  " << WD.discrepencyTime () << endl;
		commentator.report (Commentator::LEVEL_IMPORTANT, TIMING_MEASURE)
			<< "Time required for LSR fix:      " << WD.fixTime () << endl;
#endif // INCLUDE_TIMING

		commentator.stop ("done", NULL, "minpoly");

		return P;
	}
}

#include "linbox/field/modular-double.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"

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
                    FBlackbox * Ap;
                    MatrixHom::map(Ap, A, F);
                    minpoly( P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
                    delete Ap;
                    return P;
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
		if(!M.communicatorp() || (M.communicatorp())->rank() == 0) 
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
		ChineseRemainder< EarlyMultipCRA<Modular<double> > > cra(3UL);
		IntegerModularMinpoly<Blackbox,MyMethod> iteration(A, M);
		std::vector<integer> PP; // use of integer du to non genericity of cra. PG 2005-08-04
#ifdef __LINBOX_HAVE_MPI
		cra(PP, iteration, genprime, M.communicatorp());
#else
		cra(PP, iteration, genprime);
#endif
		size_t i =0;
		P.resize(PP.size());
		for (typename Polynomial::iterator it= P.begin(); it != P.end(); ++it, ++i)
			A.field().init(*it, PP[i]);

#ifdef __LINBOX_HAVE_MPI
		if(!M.communicatorp() || (M.communicatorp())->rank() == 0) 
#endif
			commentator.stop ("done", NULL, "Iminpoly");
		return P;
	}
	
} // end of LinBox namespace
#endif // __MINPOLY_H
