/* linbox/algorithms/bbchapoly.h
 * Copyright(C) LinBox
 * Written 
 *  by Clement Pernet <clement.pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_bbcharpoly_H
#define __LINBOX_bbcharpoly_H

#define __MAXITER 5
#include <vector>
#include <map>

#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/sum.h"
#include "linbox/ring/givaro-polynomial.h"
#include "linbox/field/modular.h"
#include "linbox/field/field-traits.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/blackbox/polynomial.h"

namespace LinBox 
{	
	template<class FieldPoly, class IntPoly=FieldPoly>
	class FactorMult ;

	template < class Blackbox, class Polynomial, class Categorytag >
	Polynomial&
	blackboxcharpoly (Polynomial & P, 
			  const Blackbox                   & A,
			  const Categorytag                & tag,
			  const Method::Blackbox           & M);
	
	/* Algorithm computing the integer characteristic polynomial
	 * of a blackbox.
	 */
	template < class Blackbox >
	typename GivPolynomialRing<typename Blackbox::Field>::Element&
	blackboxcharpoly (typename GivPolynomialRing<typename Blackbox::Field>::Element & P, 
			  const Blackbox                   & A,
			  const RingCategories::IntegerTag & tag,
			  const Method::Blackbox           & M)
	{
		commentator.start ("Integer Blackbox Charpoly ", "IbbCharpoly");

 		typename Blackbox::Field intRing = A.field();
 		typedef Modular<uint32> Field;
 		typedef typename Blackbox::template rebind<Field>::other FieldBlackbox;
 		typedef GivPolynomialRing<typename Blackbox::Field, Dense> IntPolyDom;
 		typedef typename IntPolyDom::Element IntPoly;
 		typedef GivPolynomialRing<Field>::Element FieldPoly;
 		// Set of factors-multiplicities sorted by degree
 		typedef FactorMult<FieldPoly,IntPoly> FM;
 		typedef multimap<unsigned long,FM*> FactPoly;
		typedef typename FactPoly::iterator FactPolyIterator;
		multimap<FM*,bool> leadingBlocks;
		//typename multimap<FM*,bool>::iterator lead_it; 
 		FactPoly factCharPoly;
 		size_t n = A.coldim();

 		IntPolyDom IPD(intRing);
		
 		/* Computation of the integer minimal polynomial */
 		IntPoly intMinPoly;
		minpoly (intMinPoly, A, M);
		if (intMinPoly.size() == n+1){
			commentator.stop ("done", NULL, "IbbCharpoly");
 			return P = intMinPoly;
		}
 		/* Factorization over the integers */
 		vector<IntPoly*> intFactors;    
 		vector<unsigned long> exp;
 		IPD.factor (intFactors, exp, intMinPoly);
		size_t factnum = intFactors.size();

 		/* Choose a modular prime field */
 		RandomPrimeIterator primeg (28);
		++primeg;
 		Field F(*primeg);

		/* Building the structure of factors */
		int goal = n;

		for (size_t i = 0; i < intFactors.size(); ++i) {
			unsigned long deg =  (intFactors[i]->size()-1);
			FactorMult<FieldPoly,IntPoly>* FFM=NULL;
			if (exp[i] > 1) {
				IntPoly *tmp = new IntPoly(*intFactors[i]);
				FM* depend = NULL;
				for (size_t j = 1; j <= exp[i]; ++j){
					IntPoly * tmp2 = new IntPoly(*tmp);
					FieldPoly *tmp2p = new FieldPoly(tmp2->size());
					typename IntPoly::template rebind<Field>() (*tmp2p, *tmp2, F);

					FFM = new FM (tmp2p, tmp2, 0, depend);
					factCharPoly.insert (pair<size_t, FM*> (deg, FFM));
					factnum++;
					depend = FFM;
					deg += intFactors[i]->size()-1;
					if (j < exp[i])
						IPD.mul (*tmp, *tmp2, *intFactors[i]);
				}
				delete tmp;
				factnum--;
				FFM->multiplicity = 1; // The last factor is present in minpoly
				goal -= deg-intFactors[i]->size()+1;
				leadingBlocks.insert (pair<FM*,bool>(FFM,false));
			} else {
				FieldPoly* fp=new FieldPoly(intFactors[i]->size());
				typename IntPoly::template rebind<Field>() (*fp, *(intFactors[i]), F);
				FFM = new FM (fp,intFactors[i],1,NULL);
				factCharPoly.insert (pair<size_t, FM* > (intFactors[i]->size()-1, FFM));
				leadingBlocks.insert (pair<FM*,bool>(FFM,false));
				goal -= deg;
			}
		}
		 
 		FieldBlackbox Ap(A, F);
		
		findMultiplicities (Ap, factCharPoly, leadingBlocks, goal, M);

 		// Building the integer charpoly
		IntPoly intCharPoly (n+1);
		IntPoly tmpP;
		intRing.init (intCharPoly[0], 1);
		for (FactPolyIterator it_f = factCharPoly.begin(); it_f != factCharPoly.end(); it_f++){
			IPD.pow (tmpP, *it_f->second->intP, it_f->second->multiplicity);
			IPD.mulin (intCharPoly, tmpP);
			delete it_f->second->intP;
			delete it_f->second->fieldP;
			delete it_f->second;
		}
		commentator.stop ("done", NULL, "IbbCharpoly");

		return P = intCharPoly;
	}
	/* Algorithm computing the  characteristic polynomial
	 * of a blackbox over a prime field.
	 */
	template < class Blackbox >
	typename GivPolynomialRing<typename Blackbox::Field>::Element& 
	blackboxcharpoly (typename GivPolynomialRing<typename Blackbox::Field>::Element & P, 
			  const Blackbox                                   & A,
			  const RingCategories::ModularTag                 & tag,
			  const Method::Blackbox           & M)
	{
		commentator.start ("Modular Blackbox Charpoly ", "MbbCharpoly");
		typedef typename Blackbox::Field Field;
		typedef GivPolynomialRing<Field, Dense> PolyDom;
		typedef typename PolyDom::Element Polynomial;
		// Set of factors-multiplicities sorted by degree
		typedef std::multimap<unsigned long,FactorMult<Polynomial>* > FactPoly;
		typedef typename FactPoly::iterator FactPolyIterator;
		multimap<FactorMult<Polynomial>*,bool> leadingBlocks;
		//typename multimap<FactorMult<Polynomial>*,bool>::iterator lead_it; 

		Field F = A.field();
		PolyDom PD (F);
		FactPoly factCharPoly;
		size_t n = A.coldim();

		/* Computation of the minimal polynomial */
		Polynomial minPoly;
		minpoly (minPoly, A, M);
		//std::cerr<<"Minpoly = "<<minPoly;
		if (minPoly.size() == n+1){
			commentator.stop ("done", NULL, "MbbCharpoly");
			return P = minPoly;
		}
		 
		/* Factorization over the field */
		std::vector<Polynomial*> factors;    
		std::vector<unsigned long> exp;
		PD.factor (factors, exp, minPoly);
		size_t factnum = factors.size();

		/* Building the structure of factors */
		int goal = n;

		for (size_t i = 0; i < factors.size(); ++i) {
			unsigned long deg =  (factors[i]->size()-1);
			FactorMult<Polynomial>* FFM=NULL;
			if (exp[i] > 1) {
				Polynomial* tmp = new Polynomial(*factors[i]);
				FactorMult<Polynomial>* depend = NULL;
				for (size_t j = 1; j <= exp[i]; ++j){
					Polynomial * tmp2 = new Polynomial(*tmp);
					FFM = new FactorMult<Polynomial> (tmp2, tmp2, 0, depend);
					//	std::cerr<<"Inserting new factor (exp>1) : "<<(*tmp2)<<std::endl;

					factCharPoly.insert (pair<size_t, FactorMult<Polynomial>*> (deg, FFM));
					factnum++;
					depend = FFM;
					deg += factors[i]->size()-1;
					if (j < exp[i])
						PD.mul (*tmp, *tmp2, *factors[i]);
				}
				delete tmp;
				factnum--;
				FFM->multiplicity = 1; // The last factor is present in minpoly
				goal -= deg-factors[i]->size()+1;
				leadingBlocks.insert (pair<FactorMult<Polynomial>*,bool>(FFM,false));
			} else {
				FFM = new FactorMult<Polynomial> (factors[i],factors[i],1,NULL);
				//std::cerr<<"Inserting new factor : "<<*factors[i]<<std::endl;
				factCharPoly.insert (pair<size_t, FactorMult<Polynomial>* > (factors[i]->size()-1, FFM));
				leadingBlocks.insert (pair<FactorMult<Polynomial>*,bool>(FFM,false));
				goal -= deg;
			}
		}
		 
		findMultiplicities ( A, factCharPoly, leadingBlocks, goal, M);
		 
		// Building the product 
		Polynomial charPoly (n+1);
		Polynomial tmpP;
		F.init (charPoly[0], 1);
		for (FactPolyIterator it_f = factCharPoly.begin(); it_f != factCharPoly.end(); it_f++){
			PD.pow (tmpP, *it_f->second->fieldP, it_f->second->multiplicity);
			PD.mulin (charPoly, tmpP);
			delete it_f->second->fieldP;
			delete it_f->second;

		}


		commentator.stop ("done", NULL, "MbbCharpoly");

		return P = charPoly;
	}
	
	template<class FieldPoly, class IntPoly>
	class FactorMult {
	public:
		FactorMult():multiplicity(0),dep(NULL){}
		FactorMult( FieldPoly* FP, IntPoly* IP, unsigned long m, FactorMult<FieldPoly,IntPoly>*d)
			: fieldP(FP), intP(IP), multiplicity(m), dep(d) {}
		FactorMult (const FactorMult<FieldPoly>& FM) 
			: fieldP(FM.fieldP), intP(FM.intP), multiplicity(FM.multiplicity), dep(FM.dep) {}
	
		int update (const size_t n, int * goal)
		{
			if (dep->dep != NULL){ 
				FactorMult<FieldPoly,IntPoly>*curr = dep;
				int k = dep->update (n,goal)+1;
				int d = (dep->fieldP->size()-1)/k;
				int tmp = (n-dep->multiplicity) / d;
				int i = k-1;
				while (curr->dep!=NULL){
					curr = curr->dep;
					tmp-=i*curr->multiplicity;
					i--;
				}
				tmp = tmp/k + (multiplicity - dep->multiplicity) / d;
				dep->multiplicity = tmp ;
				//std::cerr<<"Updating "<<*dep->fieldP<<" --> mul = "<<tmp<<std::endl;

				*goal -= tmp * (dep->fieldP->size()-1);
				return k;
			}
			else{
				int tmp =  (n - 2 * dep->multiplicity + multiplicity) / (dep->fieldP->size()-1);
				*goal -= tmp * (dep->fieldP->size()-1);
				//std::cerr<<"Updating (leaf)"<<*dep->fieldP<<" --> mul = "<<tmp<<std::endl;
				dep->multiplicity = tmp;
				return 1;
			}
		}


		FieldPoly                       *fieldP;
		IntPoly                         *intP;
		unsigned long                   multiplicity;
		FactorMult<FieldPoly, IntPoly>  *dep;

		std::ostream& write(std::ostream& os){
			return os<<"  FieldPoly --> "<<fieldP
				 <<"  IntPoly --> "<<intP
				 <<"  multiplicity --> "<<multiplicity<<std::endl
				 <<"  dep --> "<<dep<<std::endl;
		
		}
	};

	template < class FieldPoly,class IntPoly>
	void trials( list<vector<FactorMult<FieldPoly,IntPoly> > >& sols, const int goal,
		     vector<FactorMult<FieldPoly,IntPoly> >& ufv, const int i0 )
	{
		if ( !goal ){
			sols.push_back( ufv);
		}
		else if ( goal > 0 ){
			for (size_t i=i0; i<ufv.size(); ++i){
				ufv[i].multiplicity++;
				trials( sols, goal - ufv[i].fieldP->size()+1, ufv, i );
				ufv[i].multiplicity--;
			}
		}
	}	
	template <class Blackbox, class FieldPoly, class IntPoly>
	void findMultiplicities( const Blackbox& A, 
				 std::multimap<unsigned long, FactorMult<FieldPoly,IntPoly>* >& factCharPoly,
				 std::multimap<FactorMult<FieldPoly,IntPoly>*,bool>& leadingBlocks,
				 int goal,
				 const Method::Blackbox &M)
	{
		typedef multimap<unsigned long, FactorMult<FieldPoly,IntPoly>* > FactPoly;
		typedef typename Blackbox::Field Field;
		typedef GivPolynomialRing<Field, Dense> PolyDom;
		typename FactPoly::iterator itf = factCharPoly.begin();
		typename multimap<FactorMult<FieldPoly,IntPoly>*,bool>::iterator lead_it; 
		Field F = A.field();
		PolyDom PD(F);
		size_t factnum = factCharPoly.size();
		size_t n = A.coldim();
		
		/* Rank for the linear factors */
		while ( ( factnum > 1 ) && ( itf->first == 1) ){
			
			lead_it = leadingBlocks.find(itf->second);
			if ( lead_it != leadingBlocks.end())
				lead_it->second = true;
			long unsigned int r;
			
			/* The matrix Pi (A) */
			if (F.isZero (itf->second->fieldP->operator[](0))){
				rank (r, A, M) ;
			} else {
				PolynomialBB<Blackbox, FieldPoly > PA (A, *itf->second->fieldP);
				rank (r, PA,  M) ;
			}
			itf->second->multiplicity = r;

			//std::cerr<<"Rank 1 : "<<*itf->second->fieldP<<" --> "<<r<<std::endl;
			factnum--;
			itf++;
		}
		
		size_t maxIter = __MAXITER;//MAX( __MAXITER, sqrt(IspBB.mnz() ) );
		// Rank for the other factorspair
		while ( factnum > maxIter ){
			lead_it = leadingBlocks.find (itf->second);
			if ( lead_it != leadingBlocks.end())
				lead_it->second = true;
			
			PolynomialBB<Blackbox, FieldPoly > PA (A, *itf->second->fieldP);
			long unsigned int r;
			rank (r, PA,  M);
						
			itf->second->multiplicity =r;
			//std::cerr<<"Rank > 1 : "<<*itf->second->intP<<" --> "<<r<<std::endl;

			factnum--;
			itf++;
		}
		
		// update the multiplicities
		for (lead_it = leadingBlocks.begin(); lead_it != leadingBlocks.end(); lead_it++){
			
			FactorMult<FieldPoly,IntPoly>* currFFM = lead_it->first;
			//std::cerr<<"Updating multiplicities of "<<*lead_it->first->fieldP<<std::endl;
						
			if (!lead_it->second){ // the leading block has not been computed
				
				// go to the last computed multiplicity of the sequence
				while (currFFM->dep!=NULL){
					if (currFFM->dep->multiplicity != 0)
						break;
					currFFM = currFFM->dep;
				}
				if (currFFM->dep != NULL){
					
					// Need one more computation:
					PolynomialBB<Blackbox, FieldPoly > PA (A, *currFFM->fieldP);
					long unsigned int r;
					rank (r, PA, M) ;
					//std::cerr<<"extra factor : "<<*currFFM->fieldP<<" --> "<<r<<std::endl;

					int tmp = currFFM->multiplicity;
					currFFM->multiplicity = r;
					currFFM->update (n,&goal);
					currFFM->multiplicity = tmp;
				}
			} else {
				int lbm;
				if (currFFM->dep != NULL){
					
					int k = currFFM->update (n,&goal)+1;
					int d = (lead_it->first->fieldP->size()-1) / k;
					
					lbm = (n-lead_it->first->multiplicity) / d;
					currFFM = currFFM->dep;
					do{
						lbm -= currFFM->multiplicity * (currFFM->fieldP->size()-1); 
						currFFM = currFFM->dep;
					} while (currFFM!=NULL);
					lbm /= k;
					goal -= (lbm-1)*(lead_it->first->fieldP->size()-1);
				}else {
					lbm = (n-lead_it->first->multiplicity) / (lead_it->first->fieldP->size()-1);
					goal -=  (lbm-1)*(lead_it->first->fieldP->size()-1); 
				}
				lead_it->first->multiplicity = lbm;
			}
		}
		
		// Recursive search 
		typename FactPoly::iterator firstUnknowFactor = itf;
		if ( factnum <= __MAXITER ){
			std::vector<FactorMult<FieldPoly,IntPoly> > unknownFact (factnum);
			for (size_t  i = 0; i < factnum; ++i, itf++){
				unknownFact[i] = *itf->second;
			}
			std::list<vector<FactorMult<FieldPoly,IntPoly> > > sols;
			trials (sols, goal,unknownFact, 0);
			typename list<vector<FactorMult<FieldPoly,IntPoly> > >::iterator uf_it = sols.begin();

			if (sols.size()>1){
				// Evaluation of P in random gamma
				typename Field::Element d, gamma, mgamma, d2;
				typename Field::RandIter g(F);
				do
					g.random(gamma);
				while (F.isZero(gamma));
				

				//Building the matrix A + gamma.Id mod p
				F.neg( mgamma, gamma );
				ScalarMatrix<Field> gammaId( F, n, gamma ); 
				Sum<Blackbox,ScalarMatrix<Field> > Agamma(A, gammaId);

				// Compute det (A+gamma.Id)
				det (d, Agamma, M);
				if (A.rowdim()%2)
					F.negin(d);
				
				// Compute Prod(Pi(-gamma)^mi)
				typename Field::Element tmp, e;
				F.init (d2,1);
				typename FactPoly::iterator it_f=factCharPoly.begin();
				PolyDom PD (F);
				for (size_t i = 0; i < factCharPoly.size()-factnum; ++i, it_f++){
					PD.eval (tmp, *it_f->second->fieldP, mgamma);
					for (size_t j=0; j < it_f->second->multiplicity; ++j)
						F.mulin (d2, tmp);
				}
				while ( uf_it != sols.end() ){
					F.init (e,1);
					for (size_t i = 0; i < uf_it->size(); ++i){
						PD.eval( tmp, *(*uf_it)[i].fieldP, mgamma );
						for (size_t  j=0; j < (*uf_it)[i].multiplicity; ++j)
							F.mulin( e, tmp );
					}
					F.mulin( e, d2);
					if (F.areEqual(e,d))
						break;
					uf_it++;
				}
				if (uf_it == sols.end())
					std::cerr<<"FAIL:: No solutions found in recursive seach"<<endl;
			} // At this point, uf_it points on the good solution
			// update with the good multiplicities
			typename FactPoly::iterator it_f = firstUnknowFactor;
			typename std::vector<FactorMult<FieldPoly,IntPoly> >::iterator it_fm = (*uf_it).begin();
			for (; it_f != factCharPoly.end(); it_f++, it_fm++)
				it_f->second->multiplicity = it_fm->multiplicity;
		}
	}

}

#undef __MAXITER
#endif // __BBCHARPOLY_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
