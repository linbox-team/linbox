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


#include "linbox/field/modular.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/field/field-traits.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	// for specialization with respect to the DomainCategory
	template< class Blackbox, class Polynomial, class MyMethod, class DomainCategory>
	Polynomial &charpoly ( Polynomial            &P, 
			       const Blackbox        &A,
			       const DomainCategory  &tag,
			       const MyMethod        &M);

	/** Computes the characteristic polynomial of A.
	 * The characteristic polynomial of a linear operator A, represented as a
	 * black box, is computed over the ring or field of A.
	 *
	 * This implementation is essentially direct, in that it does not
	 * perform any modular reduction and reconstruction. Thus, it is not
	 * recommended that one use this function to compute the characteristic polynomial of
	 * an integer or rational matrix. One should instead use the version
	 * indicated below that uses \ref{BlackboxFactory}.
	 *
	 * @param res Field element into which to store the result
	 * @param A Black box of which to compute the characteristic polynomial
	 * @param F Field over which to compute the characteristic polynomial
	 * @param M Method traits
	 */
	template <class Blackbox, class Polynomial, class MyMethod>
	Polynomial &charpoly (Polynomial         & P, 
			      const Blackbox     & A,
			      const MyMethod     & M){
		return charpoly( P, A, typename FieldTraits<typename Blackbox::Field>::categoryTag(), M);
	}

	// Charpoly with the default method
	template<class Blackbox, class Polynomial>
	Polynomial &charpoly ( Polynomial        & P, 
			       const Blackbox    & A ){
		return charpoly( P, A, FieldTraits<typename Blackbox::Field>::categoryTag(), Method::BlasElimination());
	}


	// Instantiation for the BlasElimination Method over a finite field
	template < class Polynomial, class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::ModularTag & tag,
			      const Method::BlasElimination    & M) { 
		
		BlasBlackbox< typename Blackbox::Field > BBB (A);
		BlasMatrixDomain< typename Blackbox::Field > BMD (BBB.field());
		return BMD.charpoly (P, BBB);
	}

	// Instantiation for the BlasElimination Method over the integers
	template < class Polynomial, class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::IntegerTag & tag,
			      const Method::BlasElimination    & M) { 
		
		typedef Modular<double> Field;
		typedef typename Blackbox::template rebind<Field>::other FBlackbox;
		typedef std::vector<Element> FCharPoly;

		Polynomial IntMinPoly, IntCharPoly;
		minpoly ( IntMinPoly, A, tag, M);
		/* Factorization of the minimal polynomial over Z */
		NTL::ZZXFac_InitNumPrimes = 1;
		NTL::ZZX f;
		for (size_t i = 0; i < IntMinPoly.size(); ++i)
			NTL::SetCoeff (f, i, NTL::to_ZZ((std::string( IntMinPoly[i] )).c_str()) );
		NTL::vec_pair_ZZX_long factors;
		NTL::ZZ c;
		NTL::factor (c, factors, f);
		
		/* One modular characteristic polynomial computation */
		RandomPrime primeg (22);
		integer p;
		primeg.randomPrime (p);
		Field F(p);
		FBlackbox * fbb;
		MatrixMod::mod (fbb, A, F);
		BlasBlackbox< typename Blackbox::Field > fbbb (*fbb);
		charpoly (FCharPoly, fbbb, M);
		
		/* Determine the multiplicities */
		it_f=Factor_l_field.begin();
		FieldPolynomial currPol=P_f;
		FieldPolynomial currFact;
		FieldPolynomial q(P_f); 
		FieldPolynomial r,tmp;
		vector<int> multip(nb_factor);
		for ( int i=0; i<nb_factor; ++i ){
			//		cerr<<"Facteur "<<i<<" : "<<(*it_f)<<endl;
			currFact = (*it_f++);
			r.clear();
			int m=0;
			q=currPol;
			do{
				currPol = q;
				PolDom.divmod( q, r, currPol, currFact);
				// 			cerr<<"Apres q,r,currPol,currFact= "
				// 			    <<q<<" "<<r<<" "<<currPol<<" "<<currFact;
				m++;
			} while ( PolDom.iszero( r ) );
			multip[i] = m-1;
		}


	
	/** Compute the characteristic polynomial over {\bf Z}
	 *
	 * Compute the characteristic polynomial of a matrix, represented via 
	 * a blackBox.
	 * Perform the necessary modular reductions and
	 * reconstruct the result via Chinese remaindering or rational number
	 * reconstruction.
	 *
	 * @param P Polynomial into which to store the result
	 * @param A \ref{Blacbox} that represents the matrix
	 */

	// FIXME: Right now we only support doing this over Modular<uint32> --
	// that's probably a bad idea. There needs to be a way to get from the
	// field some idea of where a "good" choice of moduli is.
	// Dan Roche 8-6-04 Fixed using FieldTraits

}

#endif // __CHARPOLY_H
