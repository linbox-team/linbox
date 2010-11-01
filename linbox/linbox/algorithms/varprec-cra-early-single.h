/* linbox/blackbox/rational-reconstruction-base.h
 * Copyright (C) 2009 Anna Marszalek
 *
 * Written by Anna Marszalek <aniau@astronet.pl>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_varprec_cra_early_single_H
#define __LINBOX_varprec_cra_early_single_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/field/PID-integer.h"
#include "linbox/field/gmp-rational.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>
#include "linbox/algorithms/cra-early-single.h"
#include "linbox/algorithms/cra-full-multip.h"
#include "linbox/algorithms/lazy-product.h"


namespace LinBox {

/* 
 * Implements early terminated CRA with preconditioning of the result
 * factor | result should be given
 * does not detect bad factors but can run until full termination [for result] in this case
 * factor may be changed by changeFactor function, residues stored in FullMultipCRA are recomputed in this case,
 * FullMulpitCRA should consist of vectors of size 0, but errors happen.
 */


template<class Domain_Type>
struct VarPrecEarlySingleCRA: public EarlySingleCRA<Domain_Type>, FullMultipCRA<Domain_Type> {

//typedef PID_Integer Integers;
//typedef Integers::Element Integer;
typedef GMPRationalField Rationals;
typedef Rationals::Element Quotient;

	typedef Domain_Type                     Domain;
	typedef typename Domain::Element DomainElement;
	typedef VarPrecEarlySingleCRA<Domain> Self_t;

	public:
		Integer factor_;
		Integer multip_;
	
		VarPrecEarlySingleCRA(const unsigned long EARLY = DEFAULT_EARLY_TERM_THRESHOLD, const double b=0.0, const Integer& f=Integer(1), const Integer& m=Integer(1)): 
			EarlySingleCRA<Domain>(EARLY), FullMultipCRA<Domain>(b), factor_(f), multip_(m) { if (factor_ == 0 ) factor_ = 1;}

		VarPrecEarlySingleCRA(const VarPrecEarlySingleCRA& other) :
                        EarlySingleCRA<Domain>(other.EARLY_TERM_THRESHOLD), FullMultipCRA<Domain>(other.LOGARITHMIC_UPPER_BOUND), factor_(other.factor_), multip_(other.multip_) {factor_ = 1; }

		int getThreshold(int& t) {return t = EarlySingleCRA<Domain>::EARLY_TERM_THRESHOLD;}

		Integer& getModulus(Integer& m) {
			EarlySingleCRA<Domain>::getModulus(m);
			return m;
		}

		/*
		 * Residue is Result * M / F - this is the value that we are reconstructing
		 */ 
		Integer& getResidue(Integer& m) {
			EarlySingleCRA<Domain>::getResidue(m);
                        return m;
                }

		void initialize (const Integer& D, const Integer& e) {
			Integer z;
			inv(z,factor_,D);
			z*=e;
			z*=multip_;
			z%=D;
			
			EarlySingleCRA<Domain>::initialize(D, z);
			std::vector<Integer> v;
			v.push_back(e);
			FullMultipCRA<Domain>::initialize(D, v);
		}

		void initialize (const Domain& D, const DomainElement& e) {
                        DomainElement z;
                        D.init(z,factor_);
                        D.invin(z);
                        D.mulin(z,e);
			DomainElement m; D.init(m, multip_);
			D.mulin(z,m);

                        EarlySingleCRA<Domain>::initialize(D, z);
			std::vector<DomainElement> v;
			v.push_back(e);
			FullMultipCRA<Domain>::initialize(D, v);
		}

	        void progress (const Integer& D, Integer& e) {

		        Integer z;
		        // Could be much faster
		        // - do not compute twice the product of moduli
		        // - reconstruct one element of e until Early Termination,
		        //   then only, try a random linear combination.

			inv(z,factor_,D);
			z*=e;
			z*=multip_;
			z%=D;
			//z = e / factor mod D;
		  	EarlySingleCRA<Domain>::progress(D, z);

                        std::vector<Integer> v;
			v.push_back(e);
			FullMultipCRA<Domain>::progress(D, v);
		}
		
                void progress (const Domain& D, const DomainElement& e) {

			DomainElement z;
			D.init(z,factor_);
			D.invin(z);
			D.mulin(z,e);
			DomainElement m; D.init(m, multip_);
                        D.mulin(z,m);
			//z = (e/ factor mod D) 

                	// Could be much faster
			// - do not compute twice the product of moduli
			// - reconstruct one element of e until Early Termination,
			//   then only, try a random linear combination.
			
			EarlySingleCRA<Domain>::progress(D, z);

                        std::vector<DomainElement> v;
			v.push_back(e);

			FullMultipCRA<Domain>::progress(D, v);
		}

		bool terminated() {
			bool ET = EarlySingleCRA<Domain>::terminated();
			if (FullMultipCRA<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) ET = ET || FullMultipCRA<Domain>::terminated();
	                return ET;
                }


                bool noncoprime(const Integer& i) const {
                        return EarlySingleCRA<Domain>::noncoprime(i);
                }

		Integer& getFactor(Integer& f) {
			return f = factor_; 
		}

		Integer& getMultip(Integer& m) {
			return m=multip_;
		}

		Integer& getPreconditioner(Integer& f, Integer& m) {
			getMultip(m);
			getFactor(f);
			return f;	
		}	

		Quotient& getPreconditioner(Quotient& q) {
			Integer m,f;
			getPreconditioner(f,m);
			return q=Quotient(m,f);
		}

		Integer& result(Integer& r) {
			if ((FullMultipCRA<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( FullMultipCRA<Domain>::terminated() )) {
				std::vector<Integer> v;
				FullMultipCRA<Domain>::result(v);
				return r = v.front();
			} else {
				Integer z;
				EarlySingleCRA<Domain>::result(z);
				return (r=factor_*z/multip_);
			}
		}

		Quotient& result(Quotient& q) {
			if ((FullMultipCRA<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( FullMultipCRA<Domain>::terminated() )) {
	         		std::vector<Integer> v;
	                	FullMultipCRA<Domain>::result(v);
		        	return q = v.front();
			} else {
				Integer z;
				EarlySingleCRA<Domain>::result(z);//residue
				return (q=Quotient(factor_*z,multip_));
			
			}
		}

		Integer& result(Integer& num, Integer& den) {
			if ((FullMultipCRA<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( FullMultipCRA<Domain>::terminated() )) {
				den =1;
				std::vector<Integer> v;
				FullMultipCRA<Domain>::result(v);
				return num = v.front();
			} else {
				Integer z;
				EarlySingleCRA<Domain>::result(z);
				den = multip_;
				num = factor_*z;
				return num;
			}
		}
		/*
		 * changes preconditioners and recomputes the residue
		 */
		bool changePreconditioner(const Integer& f, const Integer& m) {
			if ((factor_ == f) && (multip_==m)) return EarlySingleCRA<Domain>::terminated();

			factor_ = f;
			multip_ = m;
			if (factor_==0) {factor_ = 1;}//no factor if factor==0
			
			Integer e=0;
			
			//clear CRAEarlySingle;
			EarlySingleCRA<Domain>::occurency_ = 0;
			EarlySingleCRA<Domain>::nextM_ = 1UL;
			EarlySingleCRA<Domain>::primeProd_ = 1UL;
			EarlySingleCRA<Domain>::residue_ = 0;
			
			//Computation of residue_
			
		        //std::vector< double >::iterator  _dsz_it = RadixSizes_.begin();//nie wiem
	                std::vector< LazyProduct >::iterator _mod_it = FullMultipCRA<Domain>::RadixPrimeProd_.end();// list of prime products
		        std::vector< std::vector<Integer> >::iterator _tab_it = FullMultipCRA<Domain>::RadixResidues_.end();// list of residues as vectors of size 1
		        std::vector< bool >::iterator    _occ_it = FullMultipCRA<Domain>::RadixOccupancy_.end();//flags of occupied fields
			int n= FullMultipCRA<Domain>::RadixOccupancy_.size();
		        //std::vector<Integer> ri(1); LazyProduct mi; double di;//nie wiem
			// could be much faster if max occupandy is stored
			--_mod_it; --_tab_it; --_occ_it;
			int prev_shelf=0, shelf = 0; Integer prev_residue_ =0;
			for (int i=n; i > 0 ; --i, --_mod_it, --_tab_it, --_occ_it ) {
				++shelf;
				if (*_occ_it) {
					Integer D = _mod_it->operator()();
					Integer e; 

					inv(e,factor_,D);
					e *= (_tab_it->front());
					e *=multip_;
					e %=D ;


					prev_residue_ = EarlySingleCRA<Domain>::residue_;
					EarlySingleCRA<Domain>::progress(D,e);
					
					if (prev_residue_ == EarlySingleCRA<Domain>::residue_ ) 
						EarlySingleCRA<Domain>::occurency_ = EarlySingleCRA<Domain>::occurency_ +  (shelf - prev_shelf);



					if ( EarlySingleCRA<Domain>::terminated() ) {
						return true; 
					}	
					prev_shelf = shelf;
				}
			}

			return false;
		}
							
};

} //namespace LinBox

#endif //__LINBOX_varprec_cra_early_single_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
