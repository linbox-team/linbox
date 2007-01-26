/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox/solutions/det.h
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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
#ifndef __DET_H
#define __DET_H

#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/compose.h"
#include "linbox/solutions/methods.h"
#include "linbox/blackbox/dense.h"

#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blackbox-container.h"
#include "linbox/algorithms/blackbox-container-symmetric.h"
#include "linbox/algorithms/massey-domain.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/prime-stream.h"
#include "linbox/util/debug.h"
#include "linbox/util/mpicpp.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	// for specialization with respect to the DomainCategory
	template< class Blackbox, class DetMethod, class DomainCategory> 
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element          &d, 
						const Blackbox                             &A,
						const DomainCategory                     &tag,
						const DetMethod                           &M);

	
	//error handler for rational domain
	template< class Blackbox, class DetMethod> 
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d, 
						const Blackbox                            &A,
						const RingCategories::RationalTag       &tag,
						const DetMethod                          &M)
	{
		throw LinboxError("LinBox ERROR: determinant is not yet define over a rational domain");
	}


	/** \brief Compute the determinant of A
	 *
	 * The determinant of a linear operator A, represented as a
	 * black box, is computed over the ring or field of A.
	 *
	 * @param d - Field element into which to store the result
	 * @param A - Black box of which to compute the determinant
	 * @param M - optional method.  The default is Method::Hybrid(), Other options
	 include Blackbox, Elimination, Wiedemann, and BlasElimination.  Sometimes it helps to
	 indicate properties of the matrix in the method object (for instance symmetry).
	 See class Method for details.
         \ingroup solutions
        */
	template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d, 
						const Blackbox                              &A,
						const MyMethod                           &M) 
	{
		return det(d, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

	// The det with default Method 
	template<class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d, 
						const Blackbox                               &A)
	{
		return det(d, A, Method::Hybrid());
	}

	// The det with Hybrid Method 
	template<class Blackbox>
	typename Blackbox::Field::Element &det (
						typename Blackbox::Field::Element         &d, 
						const Blackbox                            &A,
						const RingCategories::ModularTag          &tag,
						const Method::Hybrid& M)
	{
		// not yet a hybrid
		/*
		  if (useBB(A)) 
		  return det(d, A, tag, Method::Blackbox(M));
		  else
		*/
		return det(d, A, tag, Method::Elimination(M));
	}

	// The det with Hybrid Method on DenseMatrix
	template<class Field>
	typename Field::Element &det (
				      typename Field::Element         &d, 
				      const DenseMatrix<Field> 			&A,
				      const RingCategories::ModularTag          &tag,
				      const Method::Hybrid& M)
	{
		return det(d, A, tag, Method::Elimination(M));
	}

	// Forward declaration saves us from including blackbox/toeplitz.h
	template<class A, class B> class Toeplitz;
       
	// Toeplitz determinant
	template<class CField, class PField >
	typename CField::Element& det(
				      typename CField::Element& res,
				      const Toeplitz<CField,PField>& A )
	{ 
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");
		return A.det(res);
	}

	// The det with Elimination Method 
	template<class Blackbox>
	typename Blackbox::Field::Element &det (
						typename Blackbox::Field::Element         &d, 
						const Blackbox                            &A,
						const RingCategories::ModularTag          &tag,
						const Method::Elimination& M)
	{
		return det(d, A, tag, Method::BlasElimination(M));
	}

	// The det with BlackBox Method 
	template<class Blackbox>
	typename Blackbox::Field::Element &det (
						typename Blackbox::Field::Element         &d, 
						const Blackbox                            &A,
						const RingCategories::ModularTag          &tag,
						const Method::Blackbox& M)
	{
		return det(d, A, tag, Method::Wiedemann(M));
	}


	// The det with Wiedemann, finite field.
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d, 
						const Blackbox                              &A,
						const RingCategories::ModularTag          &tag,
						const Method::Wiedemann                &M) 
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");
		
		typedef typename Blackbox::Field Field;
		typedef std::vector<typename Field::Element> Polynomial;
		Field F = A.field();
		
		if(M.symmetric()) {
			commentator.start ("Symmetric Determinant", "det");
			linbox_check (A.coldim () == A.rowdim ());       
			Polynomial               phi;
			unsigned long            deg;
			typename Field::RandIter iter (F);
                    
			// Precondition here to separate the eigenvalues, so that
			// minpoly (B) = charpoly (B) with high probability
			// Here there is an extra diagonal computation
			// The probability of success is also divided by two, as 
			// diag^2 contains only squares and squares are half the total elements
			std::vector<typename Field::Element> diag (A.coldim ());
                    
			typename Field::Element pi;
			size_t i;
			size_t iternum = 1;
			do {
				F.init (pi, 1);
				for (i = 0; i < A.coldim (); i++) {
					do iter.random (diag[i]); while (F.isZero (diag[i]));
					F.mulin (pi, diag[i]);
				}
                        
				Diagonal<Field> D (F, diag);
				Compose<Blackbox,Diagonal<Field> > B0 (&A, &D);                        
				typedef Compose<Diagonal<Field>,Compose<Blackbox,Diagonal<Field> > > Blackbox1;
				Blackbox1 B(&D, &B0);

				BlackboxContainerSymmetric<Field, Blackbox1> TF (&B, F, iter);
			
				MasseyDomain<Field, BlackboxContainerSymmetric<Field, Blackbox1> > WD (&TF, M.earlyTermThreshold ());

				WD.minpoly (phi, deg);
				//                         std::cout << "\tdet: iteration # " << iternum << "\tMinpoly deg= " 
				//                                   << phi.size() << "\n" ;
				//                         std::cout << "[" ;
				//                         for(typename Polynomial::const_iterator refs =  phi.begin();
				// 			        refs != phi.end() ;
				// 				      ++refs )
				// 		          std::cout << (*refs) << " " ;
				//                         std::cout << "]" << std::endl;
			
				++iternum;
			} while ( (phi.size () < A.coldim () + 1) && ( !F.isZero (phi[0]) ) );
                    
                   
			// Divided twice since multiplied twice by the diagonal matrix
			F.div (d, phi[0], pi);
			F.divin (d, pi);
                   
			if ( (deg & 1) == 1)
				F.negin (d);
                 
			commentator.stop ("done", NULL, "det");
                    
			return d;                   
		} else {
			commentator.start ("Determinant", "det");
			linbox_check (A.coldim () == A.rowdim ());
                    
			Polynomial               phi;
			unsigned long            deg;
			typename Field::RandIter iter (F);
                    
			// Precondition here to separate the eigenvalues, so that
			// minpoly (B) = charpoly (B) with high probability
			std::vector<typename Field::Element> diag (A.coldim ());
                    
			typename Field::Element pi;
			size_t i;
			size_t iternum = 1;
			do {
				F.init (pi, 1);
				for (i = 0; i < A.coldim (); i++) {
					do iter.random (diag[i]); while (F.isZero (diag[i]));
					F.mulin (pi, diag[i]);
				}
                        
				Diagonal<Field> D (F, diag);
                        
				Compose<Blackbox,Diagonal<Field> > B (&A, &D);
                        
				typedef Compose<Blackbox,Diagonal<Field> > Blackbox1;

				BlackboxContainer<Field, Blackbox1> TF (&B, F, iter);
			
				MasseyDomain<Field, BlackboxContainer<Field, Blackbox1> > WD (&TF, M.earlyTermThreshold ());

				WD.minpoly (phi, deg);
			
				++iternum;
			} while ( (phi.size () < A.coldim () + 1) && ( !F.isZero (phi[0]) ) );
                    
			F.div (d, phi[0], pi);
                    
			if (deg & 1 == 1)
				F.negin (d);
                    
                  
			commentator.stop ("done", NULL, "det");
                    
			return d;
		}
	}



	// the det with Blas, finite field.
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d,
						const Blackbox                              &A,
						const RingCategories::ModularTag          &tag,
						const Method::BlasElimination           &M) 
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");

		typedef typename Blackbox::Field Field;
		Field F = A.field();
		
		commentator.start ("Blas Determinant", "det");

		linbox_check (A.coldim () == A.rowdim ());

		BlasMatrix<typename Field::Element> B(A);
		BlasMatrixDomain<Field> BMD(F);
		d= BMD.det(B);
		commentator.stop ("done", NULL, "det");

		return d;
	}

	
	// This should work for a DenseMatrix too ?
	/** \brief A will be modified.

	\returns d determinant of A.
	\param A this BlasBlackbox matrix will be modified in place in the process.
	\ingroup solutions
 	*/
	template <class Field>
	typename Field::Element &detin (typename Field::Element             &d,
					BlasBlackbox<Field>                   &A)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");

		Field F = A.field();
		
		commentator.start ("Determinant", "det");
		linbox_check (A.coldim () == A.rowdim ());

		BlasMatrixDomain<Field> BMD(F);
		d= BMD.detin(static_cast<BlasMatrix<typename Field::Element>& > (A));
		commentator.stop ("done", NULL, "det");

		return d;
	}
} // end of LinBox namespace 

#include "linbox/field/modular.h"
//#include "linbox/field/givaro-zpz.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"

namespace LinBox {
    
	template <class Blackbox, class MyMethod>
	struct IntegerModularDet {       
		const Blackbox &A;
		const MyMethod &M;
        
		IntegerModularDet(const Blackbox& b, const MyMethod& n) 
			: A(b), M(n) {}
        
        
		template<typename Field>
		typename Field::Element& operator()(typename Field::Element& d, const Field& F) const {
			typedef typename Blackbox::template rebind<Field>::other FBlackbox;
			FBlackbox * Ap;
			MatrixHom::map(Ap, A, F);
			det( d, *Ap, M);
			delete Ap;
			return d;
		}
	};
    

	template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &cra_det (typename Blackbox::Field::Element         &d,
						    const Blackbox                            &A,
						    const RingCategories::IntegerTag          &tag,
						    const MyMethod                            &M
#ifdef __LINBOX_HAVE_MPI
							 ,Communicator 										 *C = NULL
#endif
																								)
	{
	   //  if no parallelism or if this is the parent process
		//  begin the verbose output
#ifdef __LINBOX_HAVE_MPI
	   if(!C || C->rank() == 0)
#endif
			commentator.start ("Integer Determinant", "det");
		// 0.7213475205 is an upper approximation of 1/(2log(2))
		RandomPrime genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 
		ChineseRemainder< Modular<double> > cra(3UL);
		IntegerModularDet<Blackbox,MyMethod> iteration(A, M);
		integer dd; // use of integer due to non genericity of cra. PG 2005-08-04
		//Communicator* C = M.communicatorp();
		/*
		//  if no parallelism or if parent process
		if(!C || !C->rank()){
			//  if parallel, report size of parallel world
			if(C) std::cout << "C->size()... " << C->size() << std::endl;
			std::cout << "using cra_det With C = " << int(C) 
						 << " and dd = " << dd << std::endl;
		}
		*/
		//  will call regular cra if C=0
#ifdef __LINBOX_HAVE_MPI
		cra(dd, iteration, genprime, C);
#else
		cra(dd, iteration, genprime);
#endif

#ifdef __LINBOX_HAVE_MPI
		if(!C || C->rank() == 0){
#endif
			A.field().init(d, dd); // convert the result from integer to original type
			commentator.stop ("done", NULL, "det");
#ifdef __LINBOX_HAVE_MPI
		}
#endif
		return d;
	}

} // end of LinBox namespace

#ifdef __LINBOX_HAVE_NTL
# include "linbox/algorithms/hybrid-det.h"
# define CRA_DET lif_cra_det
#else
# define CRA_DET cra_det
#endif

namespace LinBox {
    
	template <class Blackbox, class MyMethod>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d,
						const Blackbox                            &A,
						const RingCategories::IntegerTag          &tag,
						const MyMethod                            &M)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for determinant computation\n");
		return CRA_DET(d, A, tag, M);
	}    
} // end of LinBox namespace

#ifdef __LINBOX_HAVE_MPI
namespace LinBox {
    
	template <class Blackbox>
	typename Blackbox::Field::Element &det (typename Blackbox::Field::Element         &d,
						const Blackbox                            &A,
						/*const*/ Communicator						&C)
	{
		return det(d, A, Method::Hybrid(C));
	}
}
#endif //__LINBOX_HAVE_MPI

#endif // __DET_H
