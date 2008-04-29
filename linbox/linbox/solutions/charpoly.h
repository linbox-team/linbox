/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/solutions/charpoly.h
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <clement.pernet@imag.fr>
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

#ifndef __CHARPOLY_H
#define __CHARPOLY_H



#include "linbox/solutions/methods.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"
#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"

#ifdef __LINBOX_HAVE_GIVARO
// BBcharpoly without givaropolynomials is not yet implemented
#include "linbox/algorithms/bbcharpoly.h"
#endif
// Namespace in which all LinBox library code resides

namespace LinBox
{

	// for specialization with respect to the DomainCategory
 	template< class Blackbox, class Polynomial, class MyMethod, class DomainCategory>
 	Polynomial &charpoly ( Polynomial            &P, 
 			       const Blackbox        &A,
 			       const DomainCategory  &tag,
 			       const MyMethod        &M);
	

	//error handler for rational domain
	template <class Blackbox, class Polynomial>
	Polynomial &charpoly (Polynomial& P,
			      const Blackbox& A,
			      const RingCategories::RationalTag& tag,
			      const Method::Hybrid& M)
	{
		throw LinboxError("LinBox ERROR: charpoly is not yet define over a rational domain");
	}

        /** \brief  ...using an optional Method parameter
	    \param P - the output characteristic polynomial.  If the polynomial 
	    is of degree d, this random access container has size d+1, the 0-th 
	    entry is the constant coefficient and the d-th is 1 since the charpoly 
	    is monic.
	    \param A - a blackbox matrix
	    Optional \param M - the method object.  Generally, the default
	    object suffices and the algorithm used is determined by the class of M.
	    Basic methods are Method::Blackbox, Method::Elimination, and 
	    Method::Hybrid (the default).
	    See methods.h for more options.
	    \return a reference to P.
	*/
	template <class Blackbox, class Polynomial, class MyMethod>
	Polynomial &charpoly (Polynomial         & P, 
			      const Blackbox     & A,
			      const MyMethod     & M){
		return charpoly( P, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}
}


namespace LinBox
{

	/// \brief ...using default method
	template<class Blackbox, class Polynomial>
	Polynomial &charpoly (Polynomial        & P, 
			      const Blackbox    & A)
	{
		return charpoly (P, A, Method::Hybrid());
	}

	// The charpoly with Hybrid Method 
	template<class Polynomial, class Blackbox, class DomainCategory>
	Polynomial &charpoly (Polynomial            &P, 
			      const Blackbox        &A,
			      const DomainCategory  &tag,
			      const Method::Hybrid  &M)
	{
		// not yet a hybrid
		//return charpoly(P, A, tag, Method::Blackbox(M));
		return charpoly(P, A, tag, Method::BlasElimination(M));
	}

	// The charpoly with Hybrid Method 
	template<class Polynomial, class Domain, class DomainCategory>
	Polynomial &charpoly (Polynomial            &P, 
			      const SparseMatrix<Domain>  &A,
			      const DomainCategory  &tag,
			      const Method::Hybrid  &M)
	{
		// not yet a hybrid
		return charpoly(P, A, tag, Method::Blackbox(M));
	}
	// The charpoly with Hybrid Method 
	template<class Polynomial, class Domain, class DomainCategory>
	Polynomial &charpoly (Polynomial            &P, 
			      const DenseMatrix<Domain> &A,
			      const DomainCategory  &tag,
			      const Method::Hybrid  &M)
	{
		// not yet a hybrid
		return charpoly(P, A, tag, Method::BlasElimination(M));
	}
	
	// The charpoly with Elimination Method 
	template<class Polynomial, class Blackbox, class DomainCategory>
	Polynomial &charpoly (Polynomial                &P, 
			      const Blackbox            &A,
			      const DomainCategory      &tag,
			      const Method::Elimination &M)
	{
		return charpoly(P, A, tag, Method::BlasElimination(M));
	}


	/** @brief Compute the characteristic polynomial over <bold>Z</bold><sub>p</sub>
	 *
	 * Compute the characteristic polynomial of a matrix using dense 
	 * elimination methods

	 * @param P Polynomial where to store the result
	 * @param A Blackbox representing the matrix
	 */
	template < class Polynomial, class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::ModularTag & tag,
			      const Method::BlasElimination    & M) 
	{ 	
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for characteristic polynomial computation\n");
		
		BlasBlackbox< typename Blackbox::Field > BBB (A);
		BlasMatrixDomain< typename Blackbox::Field > BMD (BBB.field());
		return BMD.charpoly (P, static_cast<BlasMatrix<typename Blackbox::Field::Element> >(BBB));
	}

}

//	#if 0
#if defined(__LINBOX_HAVE_NTL) && defined(__LINBOX_HAVE_GIVARO)

#include "linbox/algorithms/cia.h"
namespace LinBox {
	/** @brief Compute the characteristic polynomial over {\bf Z}
	 *
	 * Compute the characteristic polynomial of a matrix using dense 
	 * elimination methods

	 * @param P Polynomial where to store the result
	 * @param A \ref{Blacbox} representing the matrix
	 */


	template < class Polynomial, class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::IntegerTag & tag,
			      const Method::BlasElimination    & M) 
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for characteristic polynomial computation\n");
		GivPolynomial<typename Blackbox::Field::Element> Pg;
		return P = cia (Pg, A, M);
	}

	/** Compute the characteristic polynomial over {\bf Z}
	 *
	 * Compute the characteristic polynomial of a matrix, represented via 
	 * a blackBox.
	 * 
	 * @param P Polynomial where to store the result
	 * @param A \ref{Blacbox} representing the matrix
	 */
	template < class Polynomial, class Blackbox/*, class Categorytag*/ >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::IntegerTag & tag,
			      const Method::Blackbox           & M)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for characteristic polynomial computation\n");
		GivPolynomial<typename Blackbox::Field::Element> Pg;
		return P = blackboxcharpoly (Pg, A, tag, M);
	}


#else


#include "linbox/field/modular.h"
#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/cra-full-multip.h"
#include "linbox/algorithms/cra-early-multip.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"

namespace LinBox {
	template <class Blackbox, class MyMethod>
	struct IntegerModularCharpoly {       
		const Blackbox &A;
		const MyMethod &M;
		
		IntegerModularCharpoly(const Blackbox& b, const MyMethod& n) 
			: A(b), M(n) {}
		
		template<typename Polynomial, typename Field>
		Polynomial& operator()(Polynomial& P, const Field& F) const {
			typedef typename Blackbox::template rebind<Field>::other FBlackbox;
			FBlackbox * Ap;
			MatrixHom::map(Ap, A, F);
			charpoly( P, *Ap, typename FieldTraits<Field>::categoryTag(), M);
			integer p;
			F.characteristic(p);
			//std::cerr<<"Charpoly(A) mod "<<p<<" = "<<P;

			delete Ap;
			return P;
		}            
	};
       
	template < class Polynomial,class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::IntegerTag & tag,
			      const Method::Blackbox           & M)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for characteristic polynomial computation\n");

		commentator.start ("Integer BlackBox Charpoly : No NTL installation -> chinese remaindering", "IbbCharpoly");
		
		RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 

		// typename Blackbox::ConstRawIterator it = A.rawBegin();
// 		typename Blackbox::ConstRawIterator it_end = A.rawEnd();
// 		integer max = 1,min=0;
// 		while( it != it_end ){
// 			//      cerr<<"it="<<(*it)<<endl;
// 			if (max < (*it))
// 				max = *it; 
// 			if ( min > (*it))
// 				min = *it;
// 			it++;
// 		}
// 		if (max<-min) 
// 			max=-min;
// 		size_t n=A.coldim();		
// 		double hadamarcp = n/2.0*(log(double(n))+2*log(double(max))+0.21163275)/log(2.0);

// 		ChineseRemainder< FullMultipCRA<Modular<double> > > cra(hadamarcp);
 		ChineseRemainder< EarlyMultipCRA<Modular<double> > > cra(3UL);		

		IntegerModularCharpoly<Blackbox,Method::Blackbox> iteration(A, M);
		cra.operator() (P, iteration, genprime);
		commentator.stop ("done", NULL, "IbbCharpoly");
		return P;
	}
	template < class Polynomial,class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::IntegerTag & tag,
			      const Method::BlasElimination    & M)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for characteristic polynomial computation\n");

		commentator.start ("Integer Dense Charpoly : No NTL installation -> chinese remaindering", "IbbCharpoly");
		
		RandomPrimeIterator genprime( 26-(int)ceil(log((double)A.rowdim())*0.7213475205)); 


// 		typename Blackbox::ConstRawIterator it = A.rawBegin();
// 		typename Blackbox::ConstRawIterator it_end = A.rawEnd();
// 		integer max = 1,min=0;
// 		while( it != it_end ){
// 			//      cerr<<"it="<<(*it)<<endl;
// 			if (max < (*it))
// 				max = *it; 
// 			if ( min > (*it))
// 				min = *it;
// 			it++;
// 		}
// 		if (max<-min) 
// 			max=-min;
// 		size_t n=A.coldim();
// 		double hadamarcp = n/2.0*(log(double(n))+2*log(double(max))+0.21163275)/log(2.0);

		
//		ChineseRemainder< FullMultipCRA<Modular<double> > > cra(hadamarcp);

		ChineseRemainder< EarlyMultipCRA<Modular<double> > > cra(3UL);
		IntegerModularCharpoly<Blackbox,Method::BlasElimination> iteration(A, M);
		cra(P, iteration, genprime);
		commentator.stop ("done", NULL, "IbbCharpoly");
		return P;
	}

#endif	

	
	/** Compute the characteristic polynomial over <bold>Z</bold><sub>p</sub>.
	 *
	 * Compute the characteristic polynomial of a matrix, represented via 
	 * a blackBox.
	 * 
	 * @param P Polynomial where to store the result
	 * @param A Blackbox representing the matrix
	 */
	template < class Polynomial, class Blackbox/*, class Categorytag*/ >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::ModularTag & tag,
			      const Method::Blackbox           & M)
	{
		if (A.coldim() != A.rowdim())
			throw LinboxError("LinBox ERROR: matrix must be square for characteristic polynomial computation\n");

#ifdef __LINBOX_HAVE_GIVARO
		GivPolynomial<typename Blackbox::Field::Element> Pg;
		return P = blackboxcharpoly (Pg, A, tag, M);
#else
		return charpoly(P, A, tag, Method::BlasElimination());
#endif
	}
	

}

#endif // __CHARPOLY_H
