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
#include "linbox/blackbox/polynomial.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/sum.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/det.h"
#include "linbox/util/debug.h"
#include "NTL/ZZXFactoring.h"
#include "linbox/field/ntl-ZZ.h"
#include "linbox/field/modular.h"
#include "linbox/field/field-traits.h"
#include "linbox/element/givaro-polynomial.h"
#include "linbox/ring/givaro-polynomial.h"

// Namespace in which all LinBox library code resides

namespace LinBox
{
	// for specialization with respect to the DomainCategory
	template< class Blackbox, class Polynomial, class MyMethod, class DomainCategory>
	Polynomial &charpoly ( Polynomial            &P, 
			       const Blackbox        &A,
			       const DomainCategory  &tag,
			       const MyMethod        &M);

        /** \brief  ...using an optional Method parameter
	    \parameter P - the output characteristic polynomial.  If the polynomial 
	    is of degree d, this random access container has size d+1, the 0-th 
	    entry is the constant coefficient and the d-th is 1 since the charpoly 
	    is monic.
	    \parameter A - a blackbox matrix
	    Optional \parameter M - the method object.  Generally, the default
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
		return charpoly(P, A, tag, Method::Blackbox(M));
		//		return charpoly(P, A, tag, Method::BlasElimination(M));
	}


	// The charpoly with Hybrid Method on DenseMatrix
	// Forces the elminination method
	template<class Polynomial, class Field, class DomainCategory>
	Polynomial &charpoly (Polynomial                 &P, 
			      const DenseMatrix<Field>   &A,
			      const DomainCategory       &tag,
			      const Method::Hybrid       &M)
	{
		return charpoly(P, A, tag, Method::Elimination(M));
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


	// Instantiation for the BlasElimination Method over a finite field
	template < class Polynomial, class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::ModularTag & tag,
			      const Method::BlasElimination    & M) 
	{ 
		BlasBlackbox< typename Blackbox::Field > BBB (A);
		BlasMatrixDomain< typename Blackbox::Field > BMD (BBB.field());
		return BMD.charpoly (P, static_cast<BlasMatrix<typename Blackbox::Field::Element> >(BBB));
	}

	// Instantiation for the BlasElimination Method over the integers
	template < class Polynomial, class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::IntegerTag & tag,
			      const Method::BlasElimination    & M) { 
		
		typename Blackbox::Field intRing = A.field();
		typedef Modular<double> Field;
		typedef BlasBlackbox<Field> FBlackbox;
		typedef GivPolynomialRing<typename Blackbox::Field, Dense> IntPolyDom;
		typedef GivPolynomialRing<Field, Dense> FieldPolyDom;
		typedef typename GivPolynomialRing<typename Blackbox::Field, Dense>::Element IntPoly;
		typedef typename GivPolynomialRing<Field, Dense>::Element FieldPoly;

		IntPolyDom IPD(intRing);

		
		FieldPoly fieldCharPoly;
		/* Computation of the integer minimal polynomial */
		IntPoly intMinPoly;
		minpoly (intMinPoly, A, tag, M);
		
		/* Factorization over the integers */
		vector<IntPoly> intFactors;    
		IPD.factor (intFactors, intMinPoly);
		size_t nf = intFactors.size();

		/* One modular characteristic polynomial computation */
		RandomPrime primeg (22);
		integer p;
		primeg.randomPrime (p);
		Field F(p);
		FBlackbox * fbb;
		MatrixHom::map<Field,Blackbox> (fbb, A, F);
		charpoly ( fieldCharPoly, *fbb, M);
		
    	
		/* Determination of the multiplicities */
		FieldPolyDom FPD (F);
		vector<FieldPoly> fieldFactors (nf);
		for (size_t i = 0; i < nf; ++i)
			for (size_t j = 0; j < intFactors[i].size(); ++j)
				F.init (fieldFactors[i][j], intFactors[i][j]);
		
		FieldPoly currPol = fieldCharPoly;
		FieldPoly r,tmp,q;
		vector<int> multip (nf);
		for (size_t i = 0; i < nf; ++i) {
			//cerr<<"Facteur "<<i<<" : "<<(*it_f)<<endl;
			FieldPoly& currFact = fieldFactors[i];
			r.clear();
			int m=0;
			q=currPol;
			do{
				currPol = q;
				FPD.divmod (q, r, currPol, currFact);
				//cerr<<"Apres q,r,currPol,currFact= "
				//    <<q<<" "<<r<<" "<<currPol<<" "<<currFact;
				m++;
			} while (FPD.isZero (r));
			multip[i] = m-1;
		}
		
		IntPoly intCharPoly (A.coldim());
		intRing.init (intCharPoly[0], 1);
		for (size_t i = 0; i < nf; ++i){
			IPD.pow( P, intFactors[i], multip[i] );
			IPD.mulin( intCharPoly, P );
		}
		return P = intCharPoly;
	}

	// Instantiation for the BlackBox Method over the integers
	template < class Polynomial, class Blackbox >
	Polynomial& charpoly (Polynomial                       & P, 
			      const Blackbox                   & A,
			      const RingCategories::IntegerTag & tag,
			      const Method::Blackbox           & M) { 
		
		typename Blackbox::Field intRing = A.field();
		typedef Modular<uint32> Field;
		typedef typename Blackbox::template rebind<Field>::other FieldBlackbox;
		typedef GivPolynomialRing<typename Blackbox::Field, Dense> IntPolyDom;
		typedef GivPolynomialRing<Field, Dense> FieldPolyDom;
		typedef GivPolynomial<typename Blackbox::Field::Element> IntPoly;
		typedef  GivPolynomial<typename Field::Element> FieldPoly;
		// Pair formed by a factor and its multiplicity
		typedef std::pair<IntPoly, unsigned long> IntFactorMult;
		typedef std::pair<FieldPoly, unsigned long> FieldFactorMult;
		// Set of factors-multiplicities pairs sorted by degree
		typedef multimap<unsigned long,FieldFactorMult> FactPoly;
		typedef typename multimap<unsigned long, FieldFactorMult>::iterator FactPolyIterator;
		FactPoly factCharPoly;

		IntPolyDom IPD(intRing);
		
		/* Computation of the integer minimal polynomial */
		IntPoly intMinPoly;
		minpoly (intMinPoly, A, tag, M);
		
		/* Factorization over the integers */
		vector<IntPoly> intFactors;    
		IPD.factor (intFactors, intMinPoly);
		for (size_t i = 0; i < intFactors.size(); ++i) {
			FieldPoly * Pp = new FieldPoly(intFactors[i].size());
			//IntPoly::rebind (Pp, intFactors[i]);
			factCharPoly.insert( pair<size_t, FieldFactorMult>(intFactors[i].size()-1, FieldFactorMult(*Pp,1)) );
		}
		
		/* Choose a modular prime field */
		RandomPrime primeg (31);
		integer p;
		primeg.randomPrime (p);
		Field F(p);
		
		FieldBlackbox * Ap;
		MatrixHom::map(Ap, A, F);

		/* Rank for the linear factors */
		int factnum = intFactors.size();
		FactPolyIterator itf = factCharPoly.begin();
		size_t n = A.coldim();
		int goal = n;
				
		while ( ( factnum > 5 ) && ( itf->first == 1) ){
			/* The matrix Pi (A) */
			
			PolynomialBB<FieldBlackbox, FieldPoly > PA (Ap, &itf->second.first);
			long unsigned int r;
			rank( r, PA, M ) ;
			itf->second.second = n-r;
			goal -= (n-r);
			factnum--;
			itf++;
		}
	
		int maxIter = 5;//MAX( 5, sqrt(IspBB.mnz() ) );
		// Rank for the other factorspair
		while ( factnum > maxIter ){
			PolynomialBB<FieldBlackbox, FieldPoly > PA (Ap, &itf->second.first);
			long unsigned int r;
			rank( r, PA, M ) ;
			itf->second.second = (n-r)/(itf->second.first.size()-1);
			goal -= (n-r);
			factnum--;
			itf++;
		}
		FactPolyIterator firstUnknowFactor = itf;
		// Recursive search if feasible
		if ( factnum <= 5 ){
			std::vector<FieldFactorMult> unknownFact (factnum);
			for (int i = 0; i < factnum; ++i, itf++)
				unknownFact[i] = *itf;
			std::list<vector<FieldFactorMult> > sols;
			trials (sols, goal,unknownFact, 0);

			if (sols.size()>1){
				// Evaluation of P in random gamma
				Field::Element d, gamma, mgamma, d2;
				Field::RandIter g(F);
				do
					g.random(gamma);
				while (F.isZero(gamma));
				

				//Building the matrix A + gamma.Id mod p
				F.neg( mgamma, gamma );
				ScalarMatrix<Field> gammaId( F, n, gamma ); 
				Sum<FieldBlackbox,ScalarMatrix<Field> > Agamma(*Ap, gammaId);

				// Compute det (A+gamma.Id)
				det( d, Agamma, M );
				if (A.rowdim()%2)
					F.negin(d);
				
				Field::Element tmp, e;
				F.init (d2,1);
				FactPolyIterator it_f=factCharPoly.begin();
				FieldPolyDom FPD (F);
				for (size_t i = 0; i < factCharPoly.size()-factnum; ++i, it_f++){
					FPD.eval (tmp, it_f->second.first, mgamma);
					for (size_t j=0; j < it_f->second.second; ++j)
						F.mulin (d2, tmp);
				}
				
				typename list<vector<FieldFactorMult> >::iterator uf_it = sols.begin();
				while ( uf_it != sols.end() ){
					F.init (e,1);
					for (size_t i = 0; i < uf_it->size(); ++i){
						FPD.eval( tmp, (*uf_it)[i].first, mgamma );
						for (size_t  j=0; j < (*uf_it)[i].second; ++j)
							F.mulin( e, tmp );
					}
					F.mulin( e, d2);
					if (F.areEqual(e,d)){
						cerr<<"Trouvé la bonne multiplicité"<<endl;
						FactPolyIterator it_f = firstUnknowFactor;
						typename std::vector<FieldFactorMult>::iterator it_fm = (*uf_it).begin();
						for (; it_f != factCharPoly.end(); it_f++, it_fm++)
							it_f->second.second = it_fm->second;
						break;
					}
					uf_it++;
				}
				
			}
		}
		
		IntPoly intCharPoly (A.coldim());
		IntPoly tmpP;
		intRing.init (intCharPoly[0], 1);
		for (FactPolyIterator it_f = factCharPoly.begin(); it_f != factCharPoly.end(); it_f++){
			IPD.pow (tmpP, it_f->first, it_f->second);
			IPD.mulin (intCharPoly, tmpP);
		}
		return P = intCharPoly;
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

template < class FactorMult>
void trials( list<vector<FactorMult> >& sols, const int goal,
	     vector<FactorMult>& ufv, const int i0 ){
	//	cerr<<"trials goal="<<goal<<endl;
	if ( !goal ){
		//cerr<<"Une solution trouvee"<<endl;
		sols.push_back( ufv);
	}
	else if ( goal > 0 ){
		for (int i=i0; i<ufv.size(); ++i){
			ufv[i].second++;
			trials( sols, goal - ufv[i].first.size()+1, ufv, i );
			ufv[i].second--;
		}
	}
}
#endif // __CHARPOLY_H
