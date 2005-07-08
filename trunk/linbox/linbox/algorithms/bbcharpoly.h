/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/algorithms/bbchapoly.h
 * 
 *  by Clement Pernet <clement.pernet@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __BBCHARPOLY_H
#define __BBCHARPOLY_H

#include "linbox/ring/givaro-polynomial.h"
#include "linbox/field/modular.h"
#include "linbox/field/field-traits.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/det.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/blackbox/polynomial.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/sum.h"

namespace LinBox 
{
	
	template<class>
	class FieldFactorMult ;

	template<class,class>
	class FactorMult ;

	template < class Blackbox, class Polynomial, class Categorytag >
	Polynomial&
	blackboxcharpoly (Polynomial & P, 
			  const Blackbox                   & A,
			  const Categorytag                & tag);

	/* Algorithm computing the integer characteristic polynomial
	 * of a blackbox.
	 */
	template < class Blackbox >
	GivPolynomial<typename Blackbox::Field::Element>&
	blackboxcharpoly (GivPolynomial<typename Blackbox::Field::Element> & P, 
			  const Blackbox                   & A,
			  const RingCategories::IntegerTag & tag)
	{
		typename Blackbox::Field intRing = A.field();
		typedef Modular<uint32> Field;
		typedef typename Blackbox::template rebind<Field>::other FieldBlackbox;
		typedef GivPolynomialRing<typename Blackbox::Field, Dense> IntPolyDom;
		typedef GivPolynomialRing<Field, Dense> FieldPolyDom;
		typedef GivPolynomial<typename Blackbox::Field::Element> IntPoly;
		typedef GivPolynomial<typename Field::Element> FieldPoly;
		// Pair formed by a factor and its multiplicity
		//		typedef std::pair<FieldPoly, unsigned long> FactorMult;
		// Set of factors-multiplicities pairs sorted by degree
		typedef FactorMult<FieldPoly,IntPoly> FM;
		typedef multimap<unsigned long,FM> FactPoly;
		typedef typename multimap<unsigned long, FM>::iterator FactPolyIterator;
		FactPoly factCharPoly;
		size_t n = A.coldim();

		IntPolyDom IPD(intRing);
		
		/* Computation of the integer minimal polynomial */
		IntPoly intMinPoly;
		minpoly (intMinPoly, A, tag, M);
		
		if (intMinPoly.size() == n+1)
			return P = intMinPoly;
		/* Factorization over the integers */
		vector<IntPoly> intFactors;    
		vector<unsigned long> mult;

		IPD.factor (intFactors, mult, intMinPoly);

		/* Choose a modular prime field */
		RandomPrime primeg (31);
		integer p;
		primeg.randomPrime (p);
		Field F(p);

		int goal = n;

		for (size_t i = 0; i < intFactors.size(); ++i) {
			FieldPoly * Pp = new FieldPoly(intFactors[i].size());
			typename IntPoly::template rebind<Field>() (Pp, intFactors[i], F);
			goal -= intFactors[i].size() - 1;
			factCharPoly.insert( pair<size_t, FM>(intFactors[i].size()-1, FM(*Pp,intFactors[i],1)) );
		}
		
		
		FieldBlackbox * Ap;
		MatrixHom::map(Ap, A, F);

		/* Rank for the linear factors */
		int factnum = intFactors.size();
		FactPolyIterator itf = factCharPoly.begin();
				
		while ( ( factnum > 1 ) && ( itf->first == 1) ){
			/* The matrix Pi (A) */
			
			PolynomialBB<FieldBlackbox, FieldPoly > PA (Ap, &itf->second.fieldP);
			long unsigned int r;
			rank( r, PA, M ) ;
			itf->second.multiplicity = n-r;
			goal -= (n-r-1);
			//std::cerr<<"deg 1 factor --> goal="<<goal<<std::endl;

			factnum--;
			itf++;
		}
	
		int maxIter = 5;//MAX( 5, sqrt(IspBB.mnz() ) );
		// Rank for the other factorspair
		while ( factnum > maxIter ){
			PolynomialBB<FieldBlackbox, FieldPoly > PA (Ap, &itf->second.fieldP);
			long unsigned int r;
			rank( r, PA, M ) ;
			itf->second.multiplicity = (n-r)/(itf->second.fieldP.size()-1);
			goal -= (n-r - itf->second.fieldP.size() + 1);
			//std::cerr<<"deg >1 factor --> goal="<<goal<<std::endl;
			factnum--;
			itf++;
		}
		FactPolyIterator firstUnknowFactor = itf;
		// Recursive search if feasible
		if ( factnum <= 5 ){
			//std::std::cerr<<"trials with "<<factnum<<" unknown factors goal="<<goal<<std::endl;
			std::vector<FM> unknownFact (factnum);
			for (int i = 0; i < factnum; ++i, itf++){
				//itf->second.write(std::cerr<<"UnknownFactor["<<i<<"] = "<<std::endl);
				unknownFact[i] = itf->second;
			}
			std::list<vector<FM> > sols;
			trials (sols, goal,unknownFact, 0);

			typename list<vector<FM> >::iterator uf_it = sols.begin();
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
				
				// Compute Prod(Pi(-gamma)^mi)
				Field::Element tmp, e;
				F.init (d2,1);
				FactPolyIterator it_f=factCharPoly.begin();
				FieldPolyDom FPD (F);
				for (size_t i = 0; i < factCharPoly.size()-factnum; ++i, it_f++){
					FPD.eval (tmp, it_f->second.fieldP, mgamma);
					for (size_t j=0; j < it_f->second.multiplicity; ++j)
						F.mulin (d2, tmp);
				}
				
				while ( uf_it != sols.end() ){
					F.init (e,1);
					for (size_t i = 0; i < uf_it->size(); ++i){
						FPD.eval( tmp, (*uf_it)[i].fieldP, mgamma );
						for (size_t  j=0; j < (*uf_it)[i].multiplicity; ++j)
							F.mulin( e, tmp );
					}
					F.mulin( e, d2);
					if (F.areEqual(e,d))
						break;
					uf_it++;
				}
				
			}
			// update with the good multiplicities
			FactPolyIterator it_f = firstUnknowFactor;
			typename std::vector<FM>::iterator it_fm = (*uf_it).begin();
			for (; it_f != factCharPoly.end(); it_f++, it_fm++)
				it_f->second.multiplicity = it_fm->multiplicity;
		}

		// Building the integer charpoly
		IntPoly intCharPoly (n+1);
		IntPoly tmpP;
		intRing.init (intCharPoly[0], 1);
		for (FactPolyIterator it_f = factCharPoly.begin(); it_f != factCharPoly.end(); it_f++){
			IPD.pow (tmpP, it_f->second.intP, it_f->second.multiplicity);
			IPD.mulin (intCharPoly, tmpP);
		}
		return P = intCharPoly;
	}

	 template < class Blackbox >
	 GivPolynomial<typename Blackbox::Field::Element>& 
	 blackboxcharpoly (GivPolynomial<typename Blackbox::Field::Element> & P, 
			   const Blackbox                                   & A,
			   const RingCategories::ModularTag                 & tag,
			   const Method::Blackbox                           & M) 
	 {
		 typedef typename Blackbox::Field Field;
		 Field F = A.field();
		 typedef GivPolynomialRing<Field, Dense> PolyDom;
		 typedef GivPolynomial<typename Field::Element> Polynomial;
		 // Pair formed by a factor and its multiplicity

		 // Set of factors-multiplicities pairs sorted by degree
		 typedef multimap<unsigned long,FieldFactorMult<Polynomial>* > FactPoly;
		 typedef typename multimap<unsigned long, FieldFactorMult<Polynomial>* >::iterator FactPolyIterator;
		 FactPoly factCharPoly;

		 PolyDom PD (F);

		 /* Computation of the minimal polynomial */
		 Polynomial minPoly;
		 minpoly (minPoly, A, tag, M);
		 size_t n = A.coldim();

		 if (minPoly.size() == n+1)
			 return P = minPoly;
		 /* Factorization over the field */
		 vector<Polynomial> factors;    
		 vector<unsigned long> exp;
		 PD.factor (factors, exp, minPoly);

		 IntFactorDom<> IFD;
		 size_t factnum = factors.size();
		 //		cerr<<"factors.size()="<<factnum<<std::endl;
		 int goal = n;
		 multimap<FieldFactorMult<Polynomial>*,bool> leadingBlocks;
		 for (size_t i = 0; i < factors.size(); ++i) {
			 std::cerr<<"Traitement du facteur : "<<factors[i]<<" de multiplicite : "<<exp[i]<<std::endl;
			 unsigned long deg =  (factors[i].size()-1);
			 FieldFactorMult<Polynomial>*FFM;
			 if (exp[i] > 1) {
				 
				 Polynomial tmp = factors[i];
				 FieldFactorMult<Polynomial>*depend=NULL;
				 for (size_t j = 1; j <= exp[i]; ++j){
					 Polynomial * tmp2 = new Polynomial(tmp);
					 FFM = new FieldFactorMult<Polynomial> (*tmp2,0,depend);
					 factCharPoly.insert (pair<size_t, FieldFactorMult<Polynomial>*> (deg, FFM));
					 factnum++;
					 if (depend!=NULL)
						 std::cerr<<"insertion d'un FMdeg qui pointe sur "<<&depend<<" "<<depend->fieldP<<std::endl;
					 else
						 std::cerr<<"insertion d'un FMdeg qui pointe sur NULL"<<std::endl;
					 depend = FFM;
					 //std::cerr<<"initialisation de depend2 avec "<<FMdeg->second.fieldP<<std::endl;
					 deg += factors[i].size()-1;
					 PD.mul(tmp, *tmp2, factors[i]);
					 //std::cerr<<"apres le mul depend2 -> "<<depend2->dep->fieldP<<std::endl;
				 }
				 factnum--;
				 FFM->multiplicity = 1; // The last factor is present in minpoly
				 goal -= deg-factors[i].size()+1;
				 std::cerr<<"goal -= "<<exp[i]*deg/2<<" --> goal = "<<goal<<std::endl;
				 leadingBlocks.insert (pair<FieldFactorMult<Polynomial>*,bool>(FFM,false));
			 } else {
				 FFM = new FieldFactorMult<Polynomial> (factors[i],1,NULL);
				 factCharPoly.insert( pair<size_t, FieldFactorMult<Polynomial>* > (factors[i].size()-1, FFM));

				 leadingBlocks.insert (pair<FieldFactorMult<Polynomial>*,bool>(FFM,false));
				 goal -= deg;
				 std::cerr<<"goal -= "<<deg<<" --> goal = "<<goal<<std::endl;
			 }
		 }

		 /* Rank for the linear factors */
		 FactPolyIterator itf = factCharPoly.begin();
		 typename multimap<FieldFactorMult<Polynomial>*,bool>::iterator lead_it; 

		 while ( ( factnum > 1 ) && ( itf->first == 1) ){
	
			 lead_it = leadingBlocks.find(itf->second);
			 if ( lead_it != leadingBlocks.end())
				 lead_it->second = true;
			 long unsigned int r;

			 /* The matrix Pi (A) */
			 if (F.isZero (itf->second->fieldP[0])){
				 rank( r, A, M ) ;
			 } else {
				 PolynomialBB<Blackbox, Polynomial > PA (A, itf->second->fieldP);
				 rank( r, PA, M ) ;
			 }
			 itf->second->multiplicity = r;
			 
			 std::cerr<<"Facteur : "<<(itf->second)<<" "<<itf->second->fieldP<<" --> "<<itf->second->multiplicity<<std::endl;
			 
			 std::cerr<<"deg 1 factor --> goal="<<goal<<std::endl;

			 factnum--;
			 itf++;
		 }

		 size_t maxIter = 3;//MAX( 3, sqrt(IspBB.mnz() ) );
		 // Rank for the other factorspair
		 while ( factnum > maxIter ){
			 lead_it = leadingBlocks.find(itf->second);
			 if ( lead_it != leadingBlocks.end())
				 lead_it->second = true;
			 PolynomialBB<Blackbox, Polynomial > PA (A, itf->second->fieldP);
			 long unsigned int r;
			 rank( r, PA, M ) ;
			 std::cerr<<"rank="<<r<<std::endl;;
			 
			 itf->second->multiplicity =r;
			 
			 std::cerr<<"Facteur : "<<itf->second->fieldP<<" --> "<<itf->second->multiplicity<<std::endl;
 			 factnum--;
			 itf++;
		 }

		 // update the multiplicities
		 for (lead_it = leadingBlocks.begin(); lead_it != leadingBlocks.end(); lead_it++){
			 FieldFactorMult<Polynomial>* currFFM = lead_it->first;
		 
			 std::cerr<<"updating the multiplicity of "<<lead_it->first->fieldP<<std::endl;
			 if (!lead_it->second){ // the leading block has not been computed
				 // go to the last computed multiplicity of the sequence
				 std::cerr<<"Never visited"<<std::endl;
				 if (currFFM->dep != NULL){
					 while (currFFM->dep->multiplicity == 0){
						 std::cerr<<"on descend sur "<<currFFM->dep->fieldP<<std::endl;
						 currFFM = currFFM->dep;
					 }
					 // Need one more computation:
					 PolynomialBB<Blackbox, Polynomial > PA (A, currFFM->fieldP);
					 long unsigned int r;
					 rank( r, PA, M ) ;
					 int tmp = currFFM->multiplicity;
					 currFFM->multiplicity = r;
					 std::cerr<<"Extra Facteur : "<<currFFM->fieldP<<" --> "<<currFFM->multiplicity<<std::endl;
					 updateFactorsMult (currFFM,n,&goal);
					 currFFM->multiplicity = tmp;
				 }
			 } else {
				 int lbm = n-lead_it->first->multiplicity;
				 std::cerr<<"lbm="<<lbm<<std::endl;
				 if (currFFM->dep != NULL){
					 std::cerr<<"goal="<<goal<<std::endl;
					 updateFactorsMult (currFFM,n,&goal);
					 std::cerr<<"goal="<<goal<<std::endl;
					 currFFM = currFFM->dep;
					 do{
						 lbm -= currFFM->multiplicity * (currFFM->fieldP.size()-1); 
						 std::cerr<<"lbm="<<lbm<<std::endl;
						 currFFM = currFFM->dep;
					 } while (currFFM!=NULL);
				 }else
					 goal -=  lbm-(lead_it->first->fieldP.size()-1); 
				 lead_it->first->multiplicity = lbm/(lead_it->first->fieldP.size()-1); 
			 }
		}
		// Recursive search if feasible
		FactPolyIterator firstUnknowFactor = itf;
		if ( factnum <= 3 ){
			std::cerr<<"trials with "<<factnum<<" unknown factors goal="<<goal<<std::endl;
			std::vector<FieldFactorMult<Polynomial> > unknownFact (factnum);
			for (size_t  i = 0; i < factnum; ++i, itf++){
				std::cerr<<"Facteur unknown : "<<itf->second->fieldP<<" --> "<<itf->second->multiplicity<<std::endl;
				unknownFact[i] = *itf->second;
			}
			std::list<vector<FieldFactorMult<Polynomial> > > sols;
			trials (sols, goal,unknownFact, 0);
			std::cerr<<"trials..done"<<std::endl
			    <<sols.size()<<" solutions found"<<std::endl;
			typename list<vector<FieldFactorMult<Polynomial> > >::iterator uf_it = sols.begin();

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
				det( d, Agamma, M );
				if (A.rowdim()%2)
					F.negin(d);
				
				// Compute Prod(Pi(-gamma)^mi)
				typename Field::Element tmp, e;
				F.init (d2,1);
				FactPolyIterator it_f=factCharPoly.begin();
				PolyDom PD (F);
				for (size_t i = 0; i < factCharPoly.size()-factnum; ++i, it_f++){
					PD.eval (tmp, it_f->second->fieldP, mgamma);
					for (size_t j=0; j < it_f->second->multiplicity; ++j)
						F.mulin (d2, tmp);
				}
				std::cerr<<"searching for the good solution among "<<sols.size()<<std::endl;
				while ( uf_it != sols.end() ){
					F.init (e,1);
					for (size_t i = 0; i < uf_it->size(); ++i){
						PD.eval( tmp, (*uf_it)[i].fieldP, mgamma );
						for (size_t  j=0; j < (*uf_it)[i].multiplicity; ++j)
							F.mulin( e, tmp );
					}
					F.mulin( e, d2);
					if (F.areEqual(e,d))
						break;
					uf_it++;
				}
				
			} // At this point, uf_it points on the good solution
			// update with the good multiplicities
			FactPolyIterator it_f = firstUnknowFactor;
			typename std::vector<FieldFactorMult<Polynomial> >::iterator it_fm = (*uf_it).begin();
			for (; it_f != factCharPoly.end(); it_f++, it_fm++)
				it_f->second->multiplicity = it_fm->multiplicity;
		}
		
		// Building the product 
		Polynomial charPoly (n+1);
		Polynomial tmpP;
		F.init (charPoly[0], 1);
		for (FactPolyIterator it_f = factCharPoly.begin(); it_f != factCharPoly.end(); it_f++){
			std::cerr<<"computing "<< it_f->second->fieldP<<" power "<< it_f->second->multiplicity<<std::endl;
			PD.pow (tmpP, it_f->second->fieldP, it_f->second->multiplicity);
			PD.mulin (charPoly, tmpP);
		}
		return P = charPoly;
	}

	template<class FieldPoly, class IntPoly>
	class FactorMult {
	public:
		FactorMult( FieldPoly& FP, IntPoly& IP, unsigned long multiplicity): fieldP (FP), intP (IP), multiplicity (multiplicity) {}
		FactorMult(){};
		FactorMult (const FactorMult& FM) : fieldP(FM.fieldP), intP(FM.intP), multiplicity(FM.multiplicity){}
		FieldPoly fieldP;
		IntPoly intP;
		unsigned long multiplicity;
		std::ostream& write(std::ostream& os){
			return os<<"IntPoly --> "<<intP
				 <<"FieldPoly --> "<<fieldP
				 <<"multiplicity --> "<<multiplicity<<std::endl;
		}
	} ;

	template<class FieldPoly>
	class FieldFactorMult {
	public:
		FieldFactorMult():multiplicity(0),dep(NULL){}
		FieldFactorMult( FieldPoly FP, unsigned long m, FieldFactorMult<FieldPoly>*d)
		{multiplicity=m;dep=d; fieldP = FP;}
		FieldFactorMult (const FieldFactorMult<FieldPoly>& FM) 
			: fieldP(FM.fieldP){multiplicity=FM.multiplicity;dep = FM.dep;}
		FieldPoly fieldP;
		unsigned long multiplicity;
		FieldFactorMult<FieldPoly>* dep;

		std::ostream& write(std::ostream& os){
			return os<<"  FieldPoly --> "<<fieldP
				 <<"  multiplicity --> "<<multiplicity<<std::endl
				 <<"  dep --> "<<dep<<std::endl;
		
		}
	};
	template<class Polynomial>
	int updateFactorsMult (FieldFactorMult<Polynomial>* ffm, 
			       const size_t                 n, 
			       int                        * goal)
	{
		if (ffm->dep->dep != NULL){ 
			FieldFactorMult<Polynomial>*curr = ffm->dep;
			int k = updateFactorsMult (ffm->dep,n,goal)+1;
			int tmp = n-ffm->dep->multiplicity;
			// 			std::cerr<<"initialisation tmp ="<<tmp;
			int i = k-1;
			while (curr->dep!=NULL){
				curr = curr->dep;
				tmp-=i*curr->multiplicity;
				// 				std::cerr<<"tmp-="<<i<<"*"<<curr->multiplicity<<" --> tmp = "<<tmp<<std::endl;
				i--;
			}
			// 			std::cerr<<"tmp = tmp / "<<k<<" - "<<ffm->dep->multiplicity<<" + "<<ffm->multiplicity<<" --> tmp="<<tmp<<std::endl;
			tmp = tmp/k + ffm->multiplicity - ffm->dep->multiplicity ;
			ffm->dep->multiplicity = tmp ;
			*goal -= tmp * (ffm->dep->fieldP.size()-1);
			// 			std::cerr<<"goal -= "<<tmp<<" * "<<(ffm->dep->fieldP.size()-1)<<" --> goal = "<<*goal<<std::endl;
			// 			std::cerr<<"Set mult ("<<ffm->dep->fieldP<<") = "<<ffm->dep->multiplicity<<std::endl;
			return k+1;
		}
		else{
			int tmp =  n - 2 * ffm->dep->multiplicity + ffm->multiplicity;
			*goal -= tmp * (ffm->dep->fieldP.size()-1);
			// 			std::cerr<<"goal -= "<<tmp<<" * "<<(ffm->dep->fieldP.size()-1)<<"goal = "<<*goal<<std::endl;;
			ffm->dep->multiplicity = tmp;
			// 			std::cerr<<"Set mult ("<<ffm->dep->fieldP<<") = "<<ffm->dep->multiplicity<<std::endl;
			return 1;
		}
	}

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
				//ufv[i].write(std::cerr<<"Appel a trials avec goal="<<goal-ufv[i].fieldP.size()+1<<" et ufv[i] modifie ="<<std::endl);
				trials( sols, goal - ufv[i].fieldP.size()+1, ufv, i );
				ufv[i].multiplicity--;
			}
		}
	}
	
	
	template < class Polynomial >
	void trials( list<vector<FieldFactorMult<Polynomial> > >& sols, const int goal,
		     vector<FieldFactorMult<Polynomial> >& ufv, const int i0 )
	{
		if ( !goal ){
			sols.push_back( ufv);
		}
		else if ( goal > 0 ){
			for (size_t i=i0; i<ufv.size(); ++i){
				ufv[i].multiplicity++;
				std::cerr<<"Appel a trials avec goal="<<(goal- (int) ufv[i].fieldP.size()+1)<<" et mul(ufv[i]) modifie ="<<(ufv[i].multiplicity)<<std::endl;
				trials( sols, goal -(int) ufv[i].fieldP.size()+1, ufv, i );
				ufv[i].multiplicity--;
			}
		}
	}

}

#endif // __BBCHARPOLY_H
