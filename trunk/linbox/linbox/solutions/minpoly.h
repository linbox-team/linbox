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
#include "linbox/algorithms/massey-domain.h"     // massey recurring sequence solver
#include "linbox/algorithms/blas-domain.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/commentator.h"
#ifdef __LINBOX_HAVE_MPI
#include "linbox/util/mpicpp.h"
#endif
#include <linbox/algorithms/minpoly-integer.h>

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
		throw LinboxError("LinBox ERROR: minpoly is not yet define over a rational domain");
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
		return minpoly(P, A, tag, Method::Wiedemann(M));
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
	    	if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for minimal polynomial computation\n");

		BlasBlackbox< typename Blackbox::Field > BBB (A);
		BlasMatrixDomain< typename Blackbox::Field > BMD (BBB.field());
		return BMD.minpoly (P, static_cast<const BlasMatrix<typename Blackbox::Field::Element>& >(BBB));
	}

	// The minpoly with BlackBox Method 
	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (
			     Polynomial         &P, 
			     const Blackbox                            &A,
			     const RingCategories::ModularTag          &tag,
			     const Method::Blackbox& M)
	{
		return minpoly(P, A, tag, Method::Wiedemann (M));
	}

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
	/*

	template < class Blackbox, class Polynomial, class FieldCategoryTag>
	Polynomial &minpolySymmetric (Polynomial& P,
	const Blackbox& A,
	FieldCategoryTag tag,
	const Method::Wiedemann& M = Method::Wiedemann ());

	template < class Blackbox, class Polynomial>
	Polynomial &minpolySymmetric (Polynomial& P,
	const Blackbox& A,
	const RingCategories::ModularTag          &tag,
	const Method::Wiedemann& M = Method::Wiedemann ()) 
	{

	minpolySymmetric(P, A,  typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	return P;
	}
				 
	template < class Blackbox, class Polynomial>
	Polynomial &minpolySymmetric (Polynomial& P,
	const Blackbox& A,
	RingCategories::IntegerTag tag,
	const Method::Wiedemann& M = Method::Wiedemann ())
	{	
	typedef typename Blackbox::Field::Element Integer;
	typedef Modular<double> ModularField;
	MinPoly<Integer, ModularField>::minPoly(P, A);

	return P;
	}
	*/

	template<class Polynomial, class Blackbox>
	Polynomial &minpoly (Polynomial& P,
			     const Blackbox& A,
			     RingCategories::ModularTag tag,
			     const Method::Wiedemann& M = Method::Wiedemann ())
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for minimal polynomial computation\n");

		typedef typename Blackbox::Field Field;
		typename Field::RandIter i (A.field());
		unsigned long            deg;

		commentator.start ("Wiedemann Minimal polynomial", "minpoly");

		BlackboxContainer<Field, Blackbox> TF (&A, A.field(), i);
		MasseyDomain< Field, BlackboxContainer<Field, Blackbox> > WD (&TF, M.earlyTermThreshold ());

		WD.minpoly (P, deg);

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

	/*
	  template < class Blackbox, class Polynomial>
	  Polynomial &minpolySymmetric (Polynomial& P,
	  const Blackbox& A,
	  RingCategories::ModularTag tag,
	  const Method::Wiedemann& M = Method::Wiedemann ())
	  {
	  typedef typename Blackbox::Field Field;
	  typename Field::RandIter i (A.field());
	  unsigned long            deg;

	  commentator.start ("Minimal polynomial", "minpoly");

	  BlackboxContainerSymmetric<Field, Blackbox> TF (&A, A.field(), i);
	  MasseyDomain< Field, BlackboxContainerSymmetric<Field, Blackbox> > WD (&TF, M.earlyTermThreshold ());

	  WD.minpoly (P, deg);

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
	*/
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
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for minimal polynomial computation\n");

#ifdef __LINBOX_HAVE_MPI
		if(!M.communicatorp() || (M.communicatorp())->rank() == 0) 
			commentator.start ("Integer Minpoly", "Iminpoly");
		else{
			commentator.setMaxDepth(0);
			commentator.setMaxDetailLevel(0);
			commentator.setPrintParameters(0, 0, 0);
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
