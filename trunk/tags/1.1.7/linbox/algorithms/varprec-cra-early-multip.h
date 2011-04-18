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

#ifndef __LINBOX_varprec_cra_multip_single_H
#define __LINBOX_varprec_cra_multip_single_H

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

//typedef PID_Integer Integers;
//typedef Integers::Element Integer;

template<class Domain_Type>
struct VarPrecEarlyMultipCRA: public EarlySingleCRA<Domain_Type>, FullMultipCRA<Domain_Type> {

typedef GMPRationalField Rationals;
typedef Rationals::Element Quotient;

typedef Domain_Type                     Domain;
typedef typename Domain::Element DomainElement;
typedef VarPrecEarlyMultipCRA<Domain> Self_t;

protected:
	std::vector< Integer > vfactor_;
	std::vector< Integer > vmultip_;

	std::vector< unsigned long > randv;
public:
	VarPrecEarlyMultipCRA(const unsigned long EARLY = DEFAULT_EARLY_TERM_THRESHOLD, const double b=0.0, 
			      const std::vector<Integer>& vf = std::vector<Integer>(0), 
			      const std::vector<Integer>& vm = std::vector<Integer>(0)): 
		EarlySingleCRA<Domain>(EARLY), FullMultipCRA<Domain>(b), vfactor_(vf), vmultip_(vm) {
			for (int i=0; i < vfactor_.size(); ++i) {
				if (vfactor_[i]==0) vfactor_[i]=1;
			}
		}

	VarPrecEarlyMultipCRA(VarPrecEarlyMultipCRA& other): EarlySingleCRA<Domain>(other.EARLY_TERM_THRESHOLD), FullMultipCRA<Domain>(other.LOGARITHMIC_UPPER_BOUND), vfactor_(other.vfactor_), vmultip_(other.vmultip_) { 
		for (int i=0; i < vfactor_.size(); ++i) {
			if (vfactor_[i]==0) vfactor_[i]=1;
		}
	} 

	int getThreshold(int& t) {return t = EarlySingleCRA<Domain>::EARLY_TERM_THRESHOLD;}	

	Integer& getModulus(Integer& m) {EarlySingleCRA<Domain>::getModulus(m);return m;}
	Integer& getResidue(Integer& r) {EarlySingleCRA<Domain>::getResidue(r);return r;}
	
	template<class Vect>
	Vect& getResidue(Vect& r) {
		Vect z,vf, vm;
                FullMultipCRA<Domain>::result(z);
		
		typename Vect::const_iterator it,itf,itm;
	        Integer M; getModulus(M);
	        getPreconditioner(vf,vm);

		//r.resize(vf.size());//vector of residues
                inverse(r, vf, M);
		normproductin(r, z, M);
		normproductin(r, vm, M);
	
		//FullMultipCRA<Domain>::getResidue(r);
		return r;
	}

	template<class Vect>
	void initialize (const Integer& D, const Vect& e) {
		srand48(BaseTimer::seed());
		vfactor_.resize( e.size(),1 );
		vmultip_.resize( e.size(),1 );
	        randv. resize ( e.size() );

		for ( std::vector<unsigned long>::iterator int_p = randv. begin(); int_p != randv. end(); ++ int_p)
	                *int_p = ((unsigned long)lrand48()) % 20000 - 10000;
		
		std::vector<Integer> vz(vfactor_.size());
		inverse(vz,vfactor_,D);
		productin(vz,vmultip_,D);
		productin(vz,e,D);
	
		Integer z;
		dot(z,D,vz,randv);
		
		EarlySingleCRA<Domain>::initialize(D, z);
		FullMultipCRA<Domain>::initialize(D, e);
	}

	template<class Vect>
	void initialize (const Domain& D, Vect& e) {
		srand48(BaseTimer::seed());
		vfactor_.resize( e.size(),1 );
		vmultip_.resize( e.size(),1 );
		randv. resize ( e.size() );

		for ( std::vector<unsigned long>::iterator int_p = randv. begin(); int_p != randv. end(); ++ int_p)
                        *int_p = ((unsigned long)lrand48()) % 20000 - 10000;

        	std::vector<DomainElement> vz(vfactor_.size());
		inverse(vz,vfactor_,D);
		productin(vz,vmultip_,D);
		productin(vz,e,D);

		DomainElement z;
		dot(z,D,vz,randv);

                EarlySingleCRA<Domain>::initialize(D, z);
		FullMultipCRA<Domain>::initialize(D, e);
	}

	template<class Vect>
	void progress (const Integer& D, const Vect& e) {

	        // Could be much faster
	        // - do not compute twice the product of moduli
	        // - reconstruct one element of e until Early Termination,
	        //   then only, try a random linear combination.
		
		std::vector<Integer> vz(vfactor_.size());
                inverse(vz,vfactor_,D);
                productin(vz,vmultip_,D);
                productin(vz,e,D);

		Integer z;
                dot(z,D,vz,randv);

	  	EarlySingleCRA<Domain>::progress(D, z);
		FullMultipCRA<Domain>::progress(D, e);
	}
		
	template<class Vect>	
        void progress (const Domain& D, const Vect& e) {
		//z = (e/ factor mod D) 
               	// Could be much faster
		// - do not compute twice the product of moduli
		// - reconstruct one element of e until Early Termination,
		//   then only, try a random linear combination.
		
		std::vector<DomainElement> vz(vfactor_.size());
                inverse(vz,vfactor_,D);
                productin(vz,vmultip_,D);
                productin(vz,e,D);
                DomainElement z;
                dot(z,D,vz,randv);

		EarlySingleCRA<Domain>::progress(D, z);
		FullMultipCRA<Domain>::progress(D, e);
	}

	bool terminated() {
		bool ET = EarlySingleCRA<Domain>::terminated();
		if (FullMultipCRA<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) ET = ET || FullMultipCRA<Domain>::terminated();
                return ET;
        }

        bool noncoprime(const Integer& i) const {
	        return EarlySingleCRA<Domain>::noncoprime(i);
        }

	//Integer& getFactor(Integer& f) {
	//	return f = EarlySingleCRA<Domain>::factor_; 
	//}

	//Integer& getMultip(Integer& m) {
	//	return m=EarlySingleCRA<Domain>::multip_;
	//}

	template<class Vect>
	Vect& getFactor(Vect& vf) const {
		//Integer f = getFactor(f);
		vf.clear();// vf.resize(vfactor_.size());
		typename Vect::const_iterator it = vfactor_.begin();
		for (;it != vfactor_.end(); ++it) {
			vf.push_back((*it));
		}
		return vf;// = vfactor_;
	}


	template<class Vect>
        Vect& getMultip(Vect& vm) const {
		//Integer m = getMultip(m);
		vm.clear(); //vm.resize(vmultip_ .size());
		typename Vect::const_iterator it = vmultip_.begin();
                for (;it != vmultip_.end(); ++it) {
	                vm.push_back((*it));
	        }
		return vm;
	}

	//Integer& getPreconditioner(Integer& f, Integer& m) const {
	//	getMultip(m);
	//	return getFactor(f);	
	//}

	template<class Vect>
        Vect& getPreconditioner(Vect& vf, Vect& vm) const {
		getMultip(vm);
		getFactor(vf);
		return vf;
	}

	//Quotient& getPreconditioner(Quotient& q) {
	//	return q = EarlySingleCRA<Domain>::getPreconditioner(q);
	//}

	template<template<class, class> class Vect, template<class> class Alloc>
	Vect<Quotient, Alloc<Quotient> >& getPreconditioner(Vect<Quotient, Alloc<Quotient> >& vq) const {
		
		Vect<Integer, Alloc<Integer> > vf,vm;
	        getPreconditioner(vf,vm);
		vq.clear();

		typename Vect<Integer, Alloc<Integer> >::const_iterator itf, itm;
		for (itf = vf.begin(), itm=vm.begin() ; itf != vf.end(); ++itf,++itm) {
			vq.push_back(Quotient(*itm,*itf));
		}
		
		return vq;
	}

	template<template<class, class> class Vect, template<class> class Alloc>
	Vect<Integer, Alloc<Integer> >& result(Vect<Integer, Alloc<Integer> >& r) {
		if ((FullMultipCRA<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( FullMultipCRA<Domain>::terminated() )) {
			FullMultipCRA<Domain>::result(r);
			return r ;
		} else {
			//Integer M; getModulus(m);
			Vect<Integer, Alloc<Integer> > z,vf, vm;
			FullMultipCRA<Domain>::result(z);
			
			typename Vect<Integer, Alloc<Integer> >::const_iterator it,itf,itm;

			Integer M; getModulus(M);
			getPreconditioner(vf,vm);

			Vect<Integer, Alloc<Integer> > residue;
			getResidue(residue);
			//vector of residues

			r.clear();

			//getPreconditioner(vf,vm);
			itf = vf.begin(); itm = vm.begin();it = z.begin();

			for (; it!= z.end(); ++it, ++itf,++itm ) {
				r.push_back(*itf * (*it) / (*itm) );
			}
			return r;
		}
	}

	template<template<class,class> class Vect, template<class> class Alloc>
	Vect<Integer, Alloc<Integer> >& result(Vect<Integer, Alloc<Integer> >& num, Integer& den) {
		if ((FullMultipCRA<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( FullMultipCRA<Domain>::terminated() )) {
			FullMultipCRA<Domain>::result(num);
			den = 1;
			return num;
		} else {
			Vect<Integer, Alloc<Integer> > z,vf, vm;
			FullMultipCRA<Domain>::result(z);//vector of non prec results

			typename Vect<Integer, Alloc<Integer> >::const_iterator it,itf,itm;
			typename Vect<Integer, Alloc<Integer> >::iterator itt;

			Integer M; getModulus(M);
			getPreconditioner(vf,vm);

			Vect<Integer, Alloc<Integer> > residue;//vector of residues
			getResidue(residue);

			num.clear();

			//getPreconditioner(vf,vm);
			itf = vf.begin(); itm = vm.begin();it = residue.begin();
			den = 1; Integer old_den =1;
			for (; it!= residue.end(); ++it, ++itf,++itm ) {
				lcm(den,den,*itm);Integer d = den/old_den;

				num.push_back((*itf) * (*it));
				if (den != old_den) {
					for (itt=num.begin(); itt !=num.end()-1 ;++itt) {
						*itt *= d; 
					}
					*itt = *itt * (den/ (*itm));
				}
				old_den = den;
			}
			return num;
		}
	}
 
	template<template<class, class> class Vect, template<class> class Alloc>
	Vect<Quotient, Alloc<Quotient> >& result(Vect<Quotient, Alloc<Quotient> >& q) {
		q.clear();
		if ((FullMultipCRA<Domain>::LOGARITHMIC_UPPER_BOUND> 1.0) && ( FullMultipCRA<Domain>::terminated() )) {
	       		std::vector<Integer> vz;
	               	FullMultipCRA<Domain>::result(vz);

			typename Vect<Integer, Alloc<Integer> >::const_iterator it = vz.begin();
			for (; it!= vz.end(); ++it) {
				q.push_back(Quotient(*it,1UL));	
			}
	        	return q;
		} else {
			Vect<Integer, Alloc<Integer> > z,vf, vm;
                        FullMultipCRA<Domain>::result(z);
			typename Vect<Integer, Alloc<Integer> >::const_iterator it = z.begin(),itf,itm;

			Integer M; getModulus(M);
	                getPreconditioner(vf,vm);

			Vect<Integer, Alloc<Integer> > residue;//vector of residues
                        inverse(residue, vf, M);
                        normproductin(residue, z, M);
                        normproductin(residue, vm, M);

	                itf = vf.begin(); itm = vm.begin();

	                for (; it!= z.end(); ++it, ++itf,++itm ) {
				q.push_back(Quotient(*itf * (*it) , *itm));
			}
			return q;
		}
	}
	

	template<class Vect>
	bool changePreconditioner(const Vect& vf, const Vect& vm) {
		//Warning does not detect unchanged preconditioners !!!
		//if ((factor_ == f) && (multip_==m)) return EarlySingleCRA<Domain>::terminated();
		
		typename Vect::const_iterator itf, itm, itf2, itm2;

		vfactor_ = vf;
		for (int i=0; i < vfactor_.size(); ++i) {
                        if (vfactor_[i]==0) vfactor_[i]=1;	//if factor ==0 set no factor
                }
		vmultip_ = vm;
					
		Vect e(vfactor_.size());
			
		//clear CRAEarlySingle;
		EarlySingleCRA<Domain>::occurency_ = 0;
		EarlySingleCRA<Domain>::nextM_ = 1UL;
		EarlySingleCRA<Domain>::primeProd_ = 1UL;
		EarlySingleCRA<Domain>::residue_ = 0;
			
		//Computation of residue_
			
		//std::vector< double >::iterator  _dsz_it = RadixSizes_.begin();
	        std::vector< LazyProduct >::iterator _mod_it = FullMultipCRA<Domain>::RadixPrimeProd_.end();// list of prime products
		std::vector< std::vector<Integer> >::iterator _tab_it = FullMultipCRA<Domain>::RadixResidues_.end();// list of residues as vectors of size 1
		std::vector< bool >::iterator    _occ_it = FullMultipCRA<Domain>::RadixOccupancy_.end();//flags of occupied fields
		int n = FullMultipCRA<Domain>::RadixOccupancy_.size();
		 //std::vector<Integer> ri(1); LazyProduct mi; double di;
		 //could be much faster if max occupandy is stored
		int prev_shelf=0, shelf = 0; Integer prev_residue_=0;
		--_occ_it; --_mod_it; --_tab_it;//last elemet
		for (int i=n; i > 0; --i, --_mod_it, --_tab_it, --_occ_it ) {
			//--_occ_it; --_mod_it; --_tab_it;
			++shelf;
			if (*_occ_it) {
				Integer D = _mod_it->operator()();
				Vect e(vfactor_.size());
				inverse(e,vfactor_,D);
				productin(e,*_tab_it,D);
				productin(e,vmultip_,D);

				Integer z;
				dot(z,D, e, randv);
			

				prev_residue_ = EarlySingleCRA<Domain>::residue_;
				EarlySingleCRA<Domain>::progress(D,z);

				if (prev_residue_ == EarlySingleCRA<Domain>::residue_ ) 
					EarlySingleCRA<Domain>::occurency_ = EarlySingleCRA<Domain>::occurency_ +  (shelf - prev_shelf);
				if ( EarlySingleCRA<Domain>::terminated() ) {
					return true; 
				}	
				prev_shelf = shelf;
			}
		}
		
		return EarlySingleCRA<Domain>::terminated();

		//forall results in FullMultipCRA do 
		//	precondition
		//	progress to  EarlySingleCRA
		//	check for termination
		//recompute
	}

	bool changeVector() {
		for ( std::vector<unsigned long>::iterator int_p = randv. begin();int_p != randv. end(); ++ int_p)
	                *int_p = ((unsigned long)lrand48()) % 20000;
		return changePreconditioner(vfactor_,vmultip_);
	}

protected:
	template <template<class> class Alloc, template<class, class> class Vect1, class Vect2>
        DomainElement& dot (DomainElement& z, const Domain& D,
        		const Vect1<DomainElement, Alloc<DomainElement> >& v1,
			const Vect2& v2) {
	        D.init(z,0); DomainElement tmp;
	        typename Vect1<DomainElement, Alloc<DomainElement> >::const_iterator v1_p;
		typename Vect2::const_iterator v2_p;
	
		for (v1_p  = v1. begin(), v2_p = v2. begin();v1_p != v1. end();++ v1_p, ++ v2_p)
			D.axpyin(z, (*v1_p), D.init(tmp, (*v2_p)));
		return z;
	}

	template <template<class> class Alloc, template<class, class> class Vect1, class Vect2>
        Integer& dot (Integer& z, const Integer& D, const Vect1<Integer, Alloc<Integer> >& v1, const Vect2& v2) {
	        z = 0;
	        typename Vect1<Integer, Alloc<Integer> >::const_iterator v1_p;
	        typename Vect2::const_iterator v2_p;
	        for (v1_p  = v1. begin(), v2_p = v2. begin(); v1_p != v1. end(); ++ v1_p, ++ v2_p) {
			z = (z + (*v1_p)*(*v2_p))%D;
	        }
		return z;
	}

	template<class Vect1, class Vect2> 
        Vect1& inverse(Vect1& vz,const Vect2 vf, const Domain D) {
		vz.clear();
		typename Vect2::const_iterator it = vf.begin();

		for (; it != vf.end(); ++it) {
			DomainElement z,i;
			D.init(z,1);
			D.init(i,*it);
			if (!D.isZero(i)) D.inv(z,i); 
			vz.push_back(z);
		}
		return vz;

	}

	template<class Vect1, class Vect2>
	Vect1& inverse(Vect1& vz,const Vect2& vf, const Integer& D) {
                vz.clear();
                typename Vect2::const_iterator it = vf.begin();

                for (; it != vf.end(); ++it) {
                	Integer z=1;
			if ((*it) != 0) inv(z,*it,D);
			vz.push_back(z);
		}
		return vz;
	}	
public:
	template<class Vect1, class Vect2>
	Vect1& productin(Vect1& vz, const Vect2 vm, const Domain D) {
		typename Vect1::iterator v1_p;
                typename Vect2::const_iterator v2_p;

		for (v1_p  = vz. begin(), v2_p = vm. begin(); v1_p != vz. end(); ++ v1_p, ++ v2_p) {
			DomainElement i;
			D.init(i,*v2_p);
	                D.mulin(*v1_p , i);
	        }
	        return vz;
	}


	template<class Vect1, class Vect2>
	Vect1& productin(Vect1& vz, const Vect2 vm, const Integer& D) {
		typename Vect1::iterator v1_p;
		typename Vect2::const_iterator v2_p;

		for (v1_p  = vz. begin(), v2_p = vm. begin(); v1_p != vz. end(); ++ v1_p, ++ v2_p) {
                        *v1_p = (*v1_p) * (*v2_p);
			*v1_p = (*v1_p) % D;
		}
		return vz;
	}

	template<class Vect1, class Vect2>
        Vect1& normproductin(Vect1& vz, const Vect2 vm, const Integer& D) {
                typename Vect1::iterator v1_p;
                typename Vect2::const_iterator v2_p;
                for (v1_p  = vz. begin(), v2_p = vm. begin(); v1_p != vz. end(); ++ v1_p, ++ v2_p) {
                        *v1_p = (*v1_p) * (*v2_p);
                        *v1_p = (*v1_p) % D;
			Integer tmp = (*v1_p >0) ? *v1_p - D : *v1_p + D;
			if (absCompare(*v1_p,tmp)> 0) *v1_p = tmp;
                }
	        return vz;
	}

	template<class Vect1, class Vect2>
	        Vect1& productin(Vect1& vz, const Vect2 vm) {
                typename Vect1::iterator v1_p;
                typename Vect2::const_iterator v2_p;

		for (v1_p  = vz. begin(), v2_p = vm. begin(); v1_p != vz. end(); ++ v1_p, ++ v2_p) {
			*v1_p *= (*v2_p);
		}
		
		return vz;
	}

};

} //namespace LinBox

#endif //__LINBOX_varprec_cra_multip_single_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
