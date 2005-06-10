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


#include "linbox/matrix/blas-matrix.h"
#include "linbox/algorithms/matrix-mod.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/util/debug.h"
#include <NTL/ZZXFactoring.h>
#include "linbox/field/ntl-ZZ.h"
#include "linbox/field/modular.h"
#include "linbox/field/field-traits.h"
#include <givaro/givpoly1.h>

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
		
		typename Blackbox::Field  IntRing = A.field();
		typedef Modular<double> Field;
		typedef typename Blackbox::template rebind<Field>::other FBlackbox;
		typedef vector<typename Field::Element> FieldPolynomial;
		vector<Polynomial> intFactors;    
		vector<FieldPolynomial> fieldFactors;

		Polynomial intMinPoly, intCharPoly;
		FieldPolynomial fieldCharPoly;
		minpoly ( intMinPoly, A, tag, M);
		
		/* One modular characteristic polynomial computation */
		RandomPrime primeg (22);
		integer p;
		primeg.randomPrime (p);
		Field F(p);
		FBlackbox * fbb;
		MatrixMod::mod (fbb, A, F);
		BlasBlackbox< Field > fbbb (*fbb);
		charpoly ( fieldCharPoly, fbbb, M);
		
		/* Factorization of the minimal polynomial over Z */
		NTL::ZZXFac_InitNumPrimes = 1;
		NTL::ZZX f;
		for (size_t i = 0; i < intMinPoly.size(); ++i)
			NTL::SetCoeff (f, i, NTL::to_ZZ((std::string( intMinPoly[i] )).c_str()) );
		NTL::vec_pair_ZZX_long factors;
		NTL::ZZ c;
		NTL::factor (c, factors, f);
	
		
		NTL::ZZ t; 
		NTL_ZZ NTLIntDom;
		for (size_t i= 0; i<factors.length(); ++i) {
			intFactors[i].resize( deg(factors[i].a)+1 );
			for(unsigned long j = 0; j <= deg(factors[i].a); ++j) {
				NTL::GetCoeff(t,factors[i].a,j);
				NTLIntDom.convert( intFactors[i][j], t );
				F.init ( fieldFactors[i][j], intFactors[i][j]);
			}
		}
    	
		/* Determine the multiplicities */
		Poly1Dom<Field,Dense> PolDom (F);
		Poly1Dom<typename Blackbox::Field, Dense> IntPolDom (IntRing);
		FieldPolynomial currPol = fieldCharPoly;
		FieldPolynomial currFact;
		FieldPolynomial r,tmp,q;
		vector<int> multip (factors.length());
		for (int i = 0; i < factors.length(); ++i) {
			//cerr<<"Facteur "<<i<<" : "<<(*it_f)<<endl;
			currFact = fieldFactors[i];
			r.clear();
			int m=0;
			q=currPol;
			do{
				currPol = q;
				PolDom.divmod (q, r, currPol, currFact);
				//cerr<<"Apres q,r,currPol,currFact= "
				//    <<q<<" "<<r<<" "<<currPol<<" "<<currFact;
				m++;
			} while (PolDom.iszero (r));
			multip[i] = m-1;
		}
		intCharPoly.resize (A.coldim());
		IntRing.init (intCharPoly[0], 1);
		for (int i = 0; i < factors.lenght(); ++i){
			IntPolDom.pow( P, IntPoly[i], multip[i] );
			IntPolDom.mulin( intCharPoly, P );
		}
		return intCharPoly;
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
