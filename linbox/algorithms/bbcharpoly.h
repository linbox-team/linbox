/* linbox/algorithms/bbchapoly.h
 * Copyright(C) LinBox
 * Written
 *  by Clement Pernet <clement.pernet@imag.fr>
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file algorithms/bbcharpoly.h
 * @ingroup algorithms
 * @brief no doc.
 *
 */

#ifndef __LINBOX_bbcharpoly_H
#define __LINBOX_bbcharpoly_H

#define _LB_MAXITER 5

#include <givaro/givpoly1.h>


#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/sum.h"
#include "linbox/ring/modular.h"
#include "linbox/ring/polynomial-ring.h"
#include "linbox/field/field-traits.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/matrix-hom.h"
#include "linbox/polynomial/dense-polynomial.h"
#include "linbox/blackbox/polynomial.h"

#include <vector>
#include <map>

namespace LinBox
{

	/*! @internal
	 * @ingroup charpoly
	 * @brief  BlackBox Characteristic Polynomial
	 * @details NO DOC
	 * @warning rebind comes from Givaro !
	 */
	class BBcharpoly {
		template<class FieldPoly, class IntPoly=FieldPoly>
		class FactorMult ;

		/*! @brief No doc
		*/
		template<class FieldPoly, class IntPoly>
		class FactorMult {
		public:
			FactorMult() :
				fieldP(NULL),intP(NULL)
				,multiplicity(0),dep(NULL)
			{}
			FactorMult( FieldPoly* FP, IntPoly* IP, size_t m, FactorMult<FieldPoly,IntPoly>*d) :
				fieldP(FP), intP(IP), multiplicity(m), dep(d)
			{}

			FactorMult (const FactorMult<FieldPoly>& FM) :
				fieldP(FM.fieldP), intP(FM.intP), multiplicity(FM.multiplicity), dep(FM.dep)
			{}

			int update (const size_t n, int & goal)
			{
// std::clog << "Update(" << n << ',' << goal << ')' << std::endl;
                if (fieldP->field().isZero(fieldP->operator[](0))) {
                        // Factor is power of X
                    int tmp = (n-multiplicity) -(fieldP->size()-1);
                    goal -= tmp;
// std::clog<<"Updating (X) "<< *fieldP <<" --> mul = "<<tmp<< " d:" << (fieldP->size()-1) << ',' << multiplicity << std::endl;
                    if (dep != NULL) {					// Power of X in minpoly
                        auto leaf(dep);
                        while (leaf->dep != NULL) {
                            leaf->multiplicity = 0;
                            leaf = leaf->dep;
                        }
                        leaf->multiplicity = (size_t)tmp;
                        multiplicity = 1;
                        return tmp;
                    } else {
                        multiplicity = n - multiplicity; // Just X in minpoly
                        return multiplicity;
                    }
                }
				if (dep->dep != NULL){
					FactorMult<FieldPoly,IntPoly>*curr = dep;
					int k = dep->update (n,goal)+1;
					int d = ((int)dep->fieldP->size()-1)/k;
					int tmp = (int)(n-dep->multiplicity) / d;
					int i = k-1;
					while (curr->dep!=NULL){
						curr = curr->dep;
						tmp-=i*(int)curr->multiplicity;
						--i;
					}
					tmp = tmp/k + (int)(multiplicity - dep->multiplicity) / d;
					dep->multiplicity =(size_t) tmp ;
// std::clog<<"Updating "<<*dep->fieldP<<" --> mul = "<<tmp<<std::endl;
					goal -= tmp * ((int)dep->fieldP->size()-1);
// std::clog << "Updated(" << n << "), " << goal << ',' << k << std::endl;
					return k;
				}
				else{
                    int tmp =  (int)((n - 2 * dep->multiplicity + multiplicity) / (dep->fieldP->size()-1));
                    goal -= tmp * ((int)dep->fieldP->size()-1);
// std::clog<<"Updating (leaf)"<<*dep->fieldP<<" --> mul = "<<tmp<< " <" << dep->multiplicity << '|' << multiplicity << '>' << std::endl;
                    dep->multiplicity = (size_t)tmp;
					return 1;
				}
			}


			FieldPoly                       *fieldP;
			IntPoly                         *intP;
			size_t                   multiplicity;
			FactorMult<FieldPoly, IntPoly>  *dep;

			std::ostream& write(std::ostream& os)
			{
				return os<<"  FieldPoly(" << fieldP << ") --> "<< *fieldP
				<<"  IntPoly(" << intP << ") --> "<< *intP
				<<"  multiplicity --> "<<multiplicity<<std::endl
				<<"  dep --> "<< dep<<std::endl;

			}
		};

	public :

		BBcharpoly() {} ;
		/** Algorithm computing the characteristic polynomial
		 * of a blackbox.
		 * @warning not implemented in general cases.
		 */
		template < class BlackBox, class Polynomial, class Categorytag >
		static Polynomial&
		blackboxcharpoly (Polynomial & P,
				  const BlackBox                   & A,
				  const Categorytag                & tag,
				  const Method::Blackbox           & M);

		/** Algorithm computing the integer characteristic polynomial
		 * of a blackbox.
		 */
		template < class BlackBox, class Polynomial >
		static Polynomial&
		blackboxcharpoly (Polynomial& P,
				  const BlackBox                   & A,
				  const RingCategories::IntegerTag & tag,
				  const Method::Blackbox           & M)
		{
			commentator().start ("Integer BlackBox Charpoly ", "IbbCharpoly");

			typedef typename BlackBox::Field Ring_t;
			const Ring_t& intRing = A.field();
			typedef Givaro::Modular<uint32_t> Field;
			typedef typename BlackBox::template rebind<Field>::other FieldBlackBox;
			typedef PolynomialRing<Ring_t, Givaro::Dense> IntPolyDom;
			typedef typename IntPolyDom::Element IntPoly;
			typedef typename PolynomialRing<Field,Givaro::Dense>::Element FieldPoly;
			// Set of factors-multiplicities sorted by degree
			typedef FactorMult<FieldPoly,IntPoly> FM;
			typedef std::multimap<size_t,FM*> FactPoly;
			typedef typename FactPoly::iterator FactPolyIterator;
			std::multimap<FM*,bool> leadingBlocks;
			//typename std::multimap<FM*,bool>::iterator lead_it;
			FactPoly factCharPoly;
			size_t n = A.coldim();

			IntPolyDom IPD(intRing,'W');

			/* Computation of the integer minimal polynomial */
			Polynomial intMinPoly(intRing);
			minpoly (intMinPoly, A, M);
			if (intMinPoly.size() == n+1){
				commentator().stop ("done", NULL, "IbbCharpoly");
				return P = intMinPoly;
			}
			/* Factorization over the integers */
			std::vector<IntPoly> intFactors(0, intRing);
			std::vector<uint64_t> exp;
			IPD.factor (intFactors, exp, intMinPoly);
// 			size_t factnum = intFactors.size();

			/* Choose a modular prime field */
			PrimeIterator<IteratorCategories::HeuristicTag> primeg (FieldTraits<Field>::bestBitSize(n));
			Field F(*primeg);
			++primeg;

                        /* Building the structure of factors */
			int goal =(int) n;

			for (size_t i = 0; i < intFactors.size(); ++i) {
				uint64_t deg =  (intFactors[i].size()-1);


				FactorMult<FieldPoly,IntPoly>* FFM=NULL;
				if (exp[i] > 1) {
					IntPoly *tmp = new IntPoly(intFactors[i]);
					FM* depend = NULL;
					for (size_t j = 1; j <= exp[i]; ++j){
						IntPoly * tmp2 = new IntPoly(*tmp);
//						FieldPoly *tmp2p = new FieldPoly(tmp2->size());
//						typename IntPoly::template rebind<Field>() (*tmp2p, *tmp2, F);
						FieldPoly *tmp2p = new FieldPoly(*tmp2, F);

						FFM = new FM (tmp2p, tmp2, 0, depend);
						factCharPoly.insert (std::pair<size_t, FM*> (deg, FFM));
// 						++factnum;
						depend = FFM;
						deg += intFactors[i].size()-1;
						if (j < exp[i])
							IPD.mul (*tmp, *tmp2, intFactors[i]);
					}
					delete tmp;
// 					--factnum;
					FFM->multiplicity = 1; // The last factor is present in minpoly
					goal -= (int)deg-(int)intFactors[i].size()+1;
					leadingBlocks.insert (std::pair<FM*,bool>(FFM,false));
				}
				else {
//					FieldPoly* fp=new FieldPoly(intFactors[i].size());
//					typename IntPoly::template rebind<Field>() (*fp, (intFactors[i]), F);
                    IntPoly * ip=new IntPoly(intFactors[i]);
					FieldPoly* fp=new FieldPoly(intFactors[i], F);

					FFM = new FM (fp,ip,1,NULL);
					factCharPoly.insert (std::pair<size_t, FM* > (intFactors[i].size()-1, FFM));
					leadingBlocks.insert (std::pair<FM*,bool>(FFM,false));
					goal -= (int)deg;
				}
			}

                // CP: const_cast is here only to please clang who would otherwise fail to compile
                // claiming (mistakingly) that there is an ambiguous specialization
			FieldBlackBox Ap(A, const_cast<const Field&>(F));

			findMultiplicities (Ap, factCharPoly, leadingBlocks, goal, M);

			// Building the integer charpoly
			IntPoly intCharPoly (intRing, n+1);
			IntPoly tmpP(intRing);
			intRing.assign(intCharPoly[0], intRing.one);
			for (FactPolyIterator it_f = factCharPoly.begin(); it_f != factCharPoly.end(); ++it_f){
// it_f->second->write(std::clog << "Int Charp Extra factor: ") << std::endl;

				IPD.pow (tmpP, *it_f->second->intP, (long) it_f->second->multiplicity);
				IPD.mulin (intCharPoly, tmpP);
                delete it_f->second->intP;
                delete it_f->second->fieldP;
				delete it_f->second;
			}
			commentator().stop ("done", NULL, "IbbCharpoly");

			return P = Polynomial(A.field(), intCharPoly);
		}

		/** Algorithm computing the  characteristic polynomial
		 * of a blackbox over a prime field.
		 */
		template < class BlackBox, class Polynomial >
		static Polynomial&
		blackboxcharpoly (Polynomial & P,
				  const BlackBox                                                & A,
				  const RingCategories::ModularTag                              & tag,
				  const Method::Blackbox                                        & M)
		{
			commentator().start ("Givaro::Modular BlackBox Charpoly ", "MbbCharpoly");
			typedef typename BlackBox::Field Field;
			typedef PolynomialRing<Field, Givaro::Dense> PolyDom;
			// Set of factors-multiplicities sorted by degree
			typedef std::multimap<size_t,FactorMult<Polynomial>* > FactPoly;
			typedef typename FactPoly::iterator FactPolyIterator;
			std::multimap<FactorMult<Polynomial>*,bool> leadingBlocks;
			//typename std::multimap<FactorMult<Polynomial>*,bool>::iterator lead_it;

			Field F = A.field();
			PolyDom PD (F,'Y');
			FactPoly factCharPoly;
			const size_t n = A.coldim();

			/* Computation of the minimal polynomial */
			Polynomial minPoly(F);
			minpoly (minPoly, A, M);
// PD.write(std::clog<<"Minpoly = ",minPoly) << std::endl;

			if (minPoly.size() == n+1){
				commentator().stop ("done", NULL, "MbbCharpoly");
				return P = minPoly;
			}



			Polynomial charPoly (F);

                /* Factorization over the field */
            std::vector<Polynomial> factors;
            std::vector<uint64_t> exp;

            PD.factor (factors, exp, minPoly);
//             size_t factnum = factors.size();
// std::clog<<"factnum = "<<factnum<<std::endl;

				/* Building the structure of factors */
            int goal = (int)n;

            for (size_t i = 0; i < factors.size(); ++i) {
                uint64_t deg =  (factors[i].size()-1);

                    // Normalize factors
                for(size_t j=0; j<=deg; ++j)
                    F.divin(factors[i][j],factors[i][deg]);


                FactorMult<Polynomial>* FFM=NULL;
                if (exp[i] > 1) {
                            // Power of X, no intermediate
                    if (F.isZero(factors[i][0])) {
                            // tmp2 is deleted after charPoly construction
                        Polynomial* tmp2 = new Polynomial(factors[i]);
                        FFM = new FactorMult<Polynomial> (tmp2,tmp2,1U,NULL);
// PD.write(std::clog<<"Inserting X : ", factors[i]) <<std::endl;
                        factCharPoly.insert (std::pair<size_t, FactorMult<Polynomial>* > (factors[i].size()-1, FFM));
//                         ++factnum;
                        FactorMult<Polynomial>* depend = FFM;
                            // tmp2 is deleted after charPoly construction
                        tmp2 = new Polynomial(factors[i]);
                        PD.pow (*tmp2, factors[i], exp[i]); // X^e
                        FFM = new FactorMult<Polynomial> (tmp2, tmp2, 0, depend);
// std::clog<<"Inserting " << exp[i] << "-th power of X: "<<(*tmp2)<<std::endl;
                        factCharPoly.insert (std::pair<size_t, FactorMult<Polynomial>* > (exp[i], FFM));
                        FFM->multiplicity = 1; // The last factor is present in minpoly
                        goal -= exp[i];
                        leadingBlocks.insert (std::pair<FactorMult<Polynomial>*,bool>(FFM,false));
                    } else {
                        Polynomial* tmp = new Polynomial(factors[i]);
                        FactorMult<Polynomial>* depend = NULL;
                        for (size_t j = 1; j <= exp[i]; ++j){
                                // tmp2 is deleted after charPoly construction
                            Polynomial* tmp2 = new Polynomial(*tmp);
                            FFM = new FactorMult<Polynomial> (tmp2, tmp2, 0, depend);
// std::clog<<"Inserting new factor (exp>1): "<<(*tmp2)<<std::endl;

                            factCharPoly.insert (std::pair<size_t, FactorMult<Polynomial>* > (deg, FFM));
//                             ++factnum;
                            depend = FFM;
                            deg += factors[i].size()-1;
                            if (j < exp[i])
                                PD.mul (*tmp, *tmp2, factors[i]);
                        }
                        delete tmp;
//                         --factnum;
                        FFM->multiplicity = 1; // The last factor is present in minpoly
                        goal -= (int)(deg-factors[i].size())+1;
                        leadingBlocks.insert (std::pair<FactorMult<Polynomial>*,bool>(FFM,false));
                    }
                } else {
                        // tmp2 is deleted after charPoly construction
                    Polynomial* tmp2 = new Polynomial(factors[i]);
                    FFM = new FactorMult<Polynomial> (tmp2,tmp2,1U,NULL);
// PD.write(std::clog<<"Inserting new factor : ", factors[i]) <<std::endl;
                    factCharPoly.insert (std::pair<size_t, FactorMult<Polynomial>* > (factors[i].size()-1, FFM));
                    leadingBlocks.insert (std::pair<FactorMult<Polynomial>*,bool>(FFM,false));
                    goal -= (int)deg;
                }
            }
            findMultiplicities ( A, factCharPoly, leadingBlocks, goal, M);

				// Building the product
            Polynomial tmpP(F);
            PD.assign(charPoly, PD.one);
            for (FactPolyIterator it_f = factCharPoly.begin(); it_f != factCharPoly.end(); ++it_f){

                PD.pow (tmpP, *it_f->second->fieldP,(long) it_f->second->multiplicity);
// PD.write(PD.write(std::clog<<"New multiplicity of ", *it_f->second->fieldP) << " found : ", tmpP) <<std::endl;
                PD.mulin (charPoly, tmpP);

                delete it_f->second->fieldP;
                delete it_f->second;
            }

			commentator().stop ("done", NULL, "MbbCharpoly");

			return P = charPoly;
		}


		template < class FieldPoly,class IntPoly>
		static void trials( std::list<std::vector<FactorMult<FieldPoly,IntPoly> > >& sols,
			     const int goal,
			     std::vector<FactorMult<FieldPoly,IntPoly> >& ufv,
			     const int i0 )
		{
			if ( !goal ){
				sols.push_back( ufv);
			}
			else if ( goal > 0 ){
				for (size_t i=(size_t)i0; i<ufv.size(); ++i){
					ufv[i].multiplicity++;
					trials( sols, goal - (int)ufv[i].fieldP->size()+1, ufv, (int)i );
					ufv[i].multiplicity--;
				}
			}
		}

		template <class BlackBox, class FieldPoly, class IntPoly>
		static void findMultiplicities( const BlackBox& A,
					 std::multimap<size_t, FactorMult<FieldPoly,IntPoly>* >& factCharPoly,
					 std::multimap<FactorMult<FieldPoly,IntPoly>*,bool>& leadingBlocks,
					 int goal,
					 const Method::Blackbox &M)
		{
			typedef std::multimap<size_t, FactorMult<FieldPoly,IntPoly>* > FactPoly;
			typedef typename BlackBox::Field Field;
			typedef PolynomialRing<Field, Givaro::Dense> PolyDom;
			typename FactPoly::iterator itf = factCharPoly.begin();
			typename std::multimap<FactorMult<FieldPoly,IntPoly>*,bool>::iterator lead_it;
			Field F = A.field();
			PolyDom PD(F,'Z');
			size_t factnum = factCharPoly.size();
			size_t n = A.coldim();

			/* Rank for the linear factors */
			while ( ( factnum > 1 ) && ( itf->first == 1) ){

				lead_it = leadingBlocks.find(itf->second);
				if ( lead_it != leadingBlocks.end())
					lead_it->second = true;
				size_t r;

				/* The matrix Pi (A) */
				if (F.isZero (itf->second->fieldP->operator[](0))){
					rank (r, A, M) ;
				}
				else {
					PolynomialBB<BlackBox, FieldPoly > PA (A, *itf->second->fieldP);
					rank (r, PA,  M) ;
				}
				itf->second->multiplicity =(size_t) r;

// std::clog<<"Rank 1 : "<<*itf->second->fieldP<<" --> "<<r<<std::endl;
				--factnum;
				++itf;
			}

			size_t maxIter = _LB_MAXITER;//MAX( _LB_MAXITER, sqrt(IspBB.mnz() ) );
			// Rank for the other factorspair
			while ( factnum > maxIter ){
				lead_it = leadingBlocks.find (itf->second);
				if ( lead_it != leadingBlocks.end())
					lead_it->second = true;

				PolynomialBB<BlackBox, FieldPoly > PA (A, *itf->second->fieldP);
				size_t r;
				rank (r, PA,  M);

				itf->second->multiplicity =(size_t)r;
// std::clog<<"Rank > 1 : "<<*itf->second->intP<<" --> "<<r<<std::endl;

				--factnum;
				++itf;
			}

			// update the multiplicities
			for (lead_it = leadingBlocks.begin(); lead_it != leadingBlocks.end(); ++lead_it){

				FactorMult<FieldPoly,IntPoly>* currFFM = lead_it->first;
// std::clog<<"Updating multiplicities of "<<*lead_it->first->fieldP<< " (goal=" << goal << ')' << std::endl;

				if (!lead_it->second){ // the leading block has not been computed

					// go to the last computed multiplicity of the sequence
					while (currFFM->dep!=NULL){
						if (currFFM->dep->multiplicity != 0)
							break;
						currFFM = currFFM->dep;
					}
					if (currFFM->dep != NULL){

						// Need one more computation:
						PolynomialBB<BlackBox, FieldPoly > PA (A, *currFFM->fieldP);
						size_t r;
						rank (r, PA, M) ;
// std::clog<<"extra factor : "<<*currFFM->fieldP<<" --> "<<r<<std::endl;

						int tmp = (int)currFFM->multiplicity;
						currFFM->multiplicity =(size_t) r;
						currFFM->update (n, goal);
						currFFM->multiplicity = (size_t)tmp;
					}
				}
				else {
                    if (F.isZero(lead_it->first->fieldP->operator[](0))) {
                            // power of X
						currFFM->update (n, goal); // remaining powers of X
                    } else {
                        int lbm;
                        if (currFFM->dep != NULL){

                            int k = currFFM->update (n, goal)+1;
                            int d = (int)(lead_it->first->fieldP->size()-1) / k;

                            lbm = (int)(n-lead_it->first->multiplicity) / d;
                            currFFM = currFFM->dep;
                            do{
                                lbm -= (int)(currFFM->multiplicity * (currFFM->fieldP->size()-1));
                                currFFM = currFFM->dep;
                            } while (currFFM!=NULL);
                            lbm /= k;
                            goal -= (lbm-1)*((int)lead_it->first->fieldP->size()-1);
                        }
                        else {
                            lbm = (int)((int)(n-lead_it->first->multiplicity) / ((int)(lead_it->first->fieldP->size())-1));
                            goal -=  (lbm-1)*((int)lead_it->first->fieldP->size()-1);
                        }
                        lead_it->first->multiplicity = (size_t)lbm;
                    }
                }
// std::clog<<"Updated multiplicities of "<<*lead_it->first->fieldP<< " --> goal: " << goal << std::endl;
			}

			// Recursive search
			typename FactPoly::iterator firstUnknowFactor = itf;
			if ( factnum <= _LB_MAXITER ){
				std::vector<FactorMult<FieldPoly,IntPoly> > unknownFact (factnum);
				for (size_t  i = 0; i < factnum; ++i, ++itf){
					unknownFact[i] = *itf->second;
				}
				std::list<std::vector<FactorMult<FieldPoly,IntPoly> > > sols;
				trials (sols, goal,unknownFact, 0);
				typename std::list<std::vector<FactorMult<FieldPoly,IntPoly> > >::iterator uf_it = sols.begin();

				if (sols.size()>1){
					// Evaluation of P in random gamma
					typename Field::Element d, gamma, mgamma, d2;
					typename Field::RandIter g(F);
					do
						g.random(gamma);
					while (F.isZero(gamma));


					//Building the matrix A + gamma.Id mod p
					F.neg( mgamma, gamma );
					ScalarMatrix<Field> gammaId( F, n, n, gamma );
					Sum<BlackBox,ScalarMatrix<Field> > Agamma(A, gammaId);

					// Compute det (A+gamma.Id)
					det (d, Agamma, M);
					if (A.rowdim()%2)
						F.negin(d);

					// Compute Prod(Pi(-gamma)^mi)
					typename Field::Element tmp, e;
					F.assign(d2,F.one);
					typename FactPoly::iterator it_f=factCharPoly.begin();
					PolyDom PD_f (F,'X');
					for (size_t i = 0; i < factCharPoly.size()-factnum; ++i, ++it_f){
						PD_f.eval (tmp, *it_f->second->fieldP, mgamma);
						for (size_t j=0; j < it_f->second->multiplicity; ++j)
							F.mulin (d2, tmp);
					}
					while ( uf_it != sols.end() ){
						F.assign (e,F.one);
						for (size_t i = 0; i < uf_it->size(); ++i){
							PD_f.eval( tmp, *(*uf_it)[i].fieldP, mgamma );
							for (size_t  j=0; j < (*uf_it)[i].multiplicity; ++j)
								F.mulin( e, tmp );
						}
						F.mulin( e, d2);
						if (F.areEqual(e,d))
							break;
						++uf_it;
					}
					if (uf_it == sols.end())
						std::cerr<<"FAIL:: No solutions found in recursive seach"<<std::endl;
				} // At this point, uf_it points on the good solution
                        // update with the good multiplicities
                    typename FactPoly::iterator it_f = firstUnknowFactor;
                    typename std::vector<FactorMult<FieldPoly,IntPoly> >::iterator it_fm = uf_it->begin();
                    for (; it_f != factCharPoly.end(); ++it_f, ++it_fm)
                        it_f->second->multiplicity = it_fm->multiplicity;
			}
		}

	};
}

#undef _LB_MAXITER
#endif // __BBCHARPOLY_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
