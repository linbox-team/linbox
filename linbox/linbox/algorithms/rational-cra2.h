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


#ifndef __LINBOX_rational2_cra_H
#define __LINBOX_rational2_cra_H
#define CRATIMING


#include "linbox/field/PID-integer.h"

#include <linbox/algorithms/rational-reconstruction-base.h>
#include <linbox/algorithms/classic-rational-reconstruction.h>

//#define RCRATIMING

namespace LinBox 
{ 

#if 0
	 template<class T, template <class T> class Container>
	 std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
		 for(typename Container<T>::const_iterator refs =  C.begin();
			 refs != C.end() ;
			 ++refs )
			 o << (*refs) << " " ;
		 return o << std::endl;
	 }
#endif


	/** \brief Chinese remainder of rationals
	 *
	 * Compute the reconstruction of rational numbers
	 * Either by Early Termination see [Dumas, Saunder, Villard, JSC 32 (1/2), pp 71-99, 2001],
	 * Or via a bound on the size of the Integers.
	 */
//typedef PID_Integer Integers;
//typedef Integers::Element Integer;

template<class RatCRABase, class RatRecon = RReconstruction<PID_integer, ClassicMaxQRationalReconstruction<PID_integer> > >
    struct RationalRemainder2 {

	
	typedef typename RatCRABase::Domain		Domain;
        typedef typename RatCRABase::DomainElement	DomainElement;
    protected:
        RatCRABase Builder_;
	RatRecon RR_;        
    public:

	int IterCounter;

        template<class Param>
        RationalRemainder2(const Param& b, const RatRecon& RR = RatRecon()) : Builder_(b), RR_(RR) { 
		IterCounter = 0;
	}

	RationalRemainder2(RatCRABase b, const RatRecon& RR = RatRecon()) : Builder_(b), RR_() {
                IterCounter = 0;
        }

           /** \brief The Rational CRA loop
				
            Given a function to generate residues mod a single prime, this loop produces the residues 
            resulting from the Chinese remainder process on sufficiently many primes to meet the 
            termination condition.
			
            \parameter F - Function object of two arguments, F(r, p), given prime p it outputs residue(s) r.
            This loop may be parallelized.  F must be reentrant, thread safe.
            For example, F may be returning the coefficients of the minimal polynomial of a matrix mod p.
            Warning - we won't detect bad primes.
			
            \parameter genprime - RandIter object for generating primes.
            \result num - the rational numerator
            \result den - the rational denominator
            */
        template<class Function, class RandPrimeIterator>
        Integer & operator() (Integer& num, Integer& den, Function& Iteration, RandPrimeIterator& genprime) {
	
            ++genprime;
            Domain D(*genprime); 
            DomainElement r; D.init(r);
            Builder_.initialize( D, Iteration(r, D) );				
	    ++IterCounter;

            int coprime =0;
	    int maxnoncoprime = 1000;

	    Integer f_in,m_in;
	    Builder_.getPreconditioner(f_in,m_in);

	    while( ! Builder_.terminated() ) {

		++genprime;
		while(Builder_.noncoprime(*genprime) ) {
			++genprime;
			++coprime;
 	                if (coprime > maxnoncoprime) {
		                cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
		                Integer res; return Builder_.result(res);
		        }
		}
		coprime = 0;
                Domain D(*genprime); 
                DomainElement r; D.init(r);
                Builder_.progress( D, Iteration(r, D) );
		if (RR_.scheduled(IterCounter-1)) {
			Integer M ; Builder_.getModulus(M);
			Integer r ; Builder_.getResidue(r);
			if (RR_.reconstructRational(num,den,r,M)) {
				Builder_.changePreconditioner(f_in*num,m_in*den);
				int k ; Builder_.getThreshold(k);
				if (this->operator()(k,num,den,Iteration,genprime)) break; 
				else {
					Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results			
				}
			}
		} 
		++IterCounter;
            }
	    Integer g;
	    Builder_.result(num,den);
	    if (gcd(g,num,den) != 1) { num /=g; den/=g; }
            return num; //Builder_.result(num, den);
        }

/* 
 * progress for k>=0 iterations
 * run until terminated() if k<0   
 * no rational reconstruction!
 */

        template<class Function, class RandPrimeIterator>
        bool operator() (const int k, Integer& num, Integer& den, Function& Iteration, RandPrimeIterator& genprime) {

	    if ((IterCounter==0) && (k != 0)) {
		++IterCounter;
                ++genprime;
	        Domain D(*genprime); 
    	        DomainElement r; D.init(r);
        	Builder_.initialize( D, Iteration(r, D) );				
	    }
            
            int coprime =0;
	    int maxnoncoprime = 1000;

	    Integer f_in,m_in;
            Builder_.getPreconditioner(f_in,m_in);

	    for ( int i=0; ((k< 0) && ! Builder_.terminated()) || (i<k) ; ++i ) {

		++genprime;

		while(Builder_.noncoprime(*genprime) ) {
			++genprime;
			++coprime;
 	                if (coprime > maxnoncoprime) {
		                cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
		                //Integer res; return Builder_.result(res);
				return false;
		        }
		}
		coprime = 0;
                Domain D(*genprime); 
                DomainElement r; D.init(r);
                Builder_.progress( D, Iteration(r, D) );
		//if (RR_.scheduled(IterCounter-1)) {
		++IterCounter;
/*	        
		Integer M ; Builder_.getModulus(M);
	        
		Integer r ; Builder_.getResidue(r);
	        
		if (RR_.reconstructRational(num,den,r,M)) {
		
		if (Builder_.changePreconditioner(f_in*num,m_in*den)) break;
                
		else Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
		
		
		}
*/		//}
	    }

            if (Builder_.terminated() ) {
	            Integer g;
	            //Builder_.getPreconditioner(p_div,p_mul);
	            Builder_.result(num,den);
	            //num = p_div*res;
	            //den = p_mul;
	            if (gcd(g,num,den) != 1) { num /=g; den/=g; }		 
	    	return true;
	    }
	    else return false;
        }

	    template<template <class, class> class Vect, template<class> class Alloc,  class Function, class RandPrimeIterator>
	    Vect<Integer, Alloc<Integer> > & operator() (Vect<Integer, Alloc<Integer> >& num, Integer& den, Function& Iteration, RandPrimeIterator& genprime) {
	    	++IterCounter;
            	++genprime;
            	Domain D(*genprime); 
            	Vect<DomainElement, Alloc<DomainElement>  > r; 
            	Builder_.initialize( D, Iteration(r, D) );

            	int coprime =0;
	    	int maxnoncoprime = 1000;

		Vect<Integer, Alloc<Integer> > f_in,m_in;
            	Builder_.getPreconditioner(f_in,m_in);

            	//while( ! Builder_.terminated() ) {
		while (1) { // in case of terminated() - checks for RR of the whole vector
                	//++IterCounter;
			++genprime; 
                	while(Builder_.noncoprime(*genprime) ) {
		        	++genprime;
                		++coprime;
				if (coprime > maxnoncoprime) {
					cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
					return num;
				}
			}
			coprime = 0;


			Domain D(*genprime); 
                	Vect<DomainElement, Alloc<DomainElement> > r; 
               		Builder_.progress( D, Iteration(r, D) );

			if (RR_.scheduled(IterCounter-1) || Builder_.terminated()) {
				Integer M ; Builder_.getModulus(M);
				if ( Builder_.terminated() ) {//early or full termination occurred, check reconstruction of the whole vector
					//early or full termination
					Vect<Integer, Alloc<Integer> > r ; Builder_.getResidue(r);
					if (RR_.reconstructRational(num,den,r,M) ) {
						Vect<Integer, Alloc<Integer> > vnum(num),vden(m_in.size(),den);
						for (int i=0; i < vnum.size(); ++ i) {
							if (vnum[i]==0) vnum[i] = 1; // no prec
						}
						Builder_.productin(vnum, f_in); Builder_.productin(vden,m_in);
						Builder_.changePreconditioner(vnum,vden) ;
						int k ; Builder_.getThreshold(k);
						if (this->operator()(k,num,den,Iteration,genprime)) {
							break;
						}
						else {	// back to original preconditioners
							Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
							Builder_.changeVector();
						}
					} else {  //back to original preconditioners
						Builder_.changePreconditioner(f_in,m_in);
						Builder_.changeVector();
					}
				} else {
					//heuristics: reconstruction of vector
					Integer r ; Builder_.getResidue(r);
					Integer n,d;
					if (RR_.reconstructRational(n,d,r,M)) {
						Vect<Integer, Alloc<Integer> > vden(m_in.size(),d);
						Builder_.productin(vden,m_in);
						Builder_.changePreconditioner(f_in,vden); 
						int k; Builder_.getThreshold(k);
						if (this->operator()(k,num,den,Iteration,genprime)) { //prob. certify result of RR
							m_in = vden;
						}
						else {	//false result of RR	
							Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
						}
								
					}
				}
			}
			++IterCounter;
		}
		Builder_.result(num,den);
            
            	return num;
        }

/* 
 * progress for k>=0 iterations
 * run until terminated if k <0
 */
	template<template <class, class> class Vect, template<class> class Alloc, class Function, class RandPrimeIterator>
	bool operator() (const int k, Vect<Integer, Alloc<Integer>  >& num, Integer& den, Function& Iteration, RandPrimeIterator& genprime) {
    	    if ((IterCounter==0) && (k != 0)) {
		++IterCounter;
        	++genprime;
		Domain D(*genprime); 
    	        Vect<DomainElement, Alloc<DomainElement>  > r;
    		Builder_.initialize( D, Iteration(r, D) );				
	    }            
            int coprime =0;
	    int maxnoncoprime = 1000;

            Vect<Integer, Alloc<Integer>  > f_in,m_in;
            Builder_.getPreconditioner(f_in,m_in);
            for (int i=0; ((k<0) && Builder_.terminated()) || (i <k); ++i ) {
                //++IterCounter;
		++genprime; 
                while(Builder_.noncoprime(*genprime) ) {
		        ++genprime;
                	++coprime;
			if (coprime > maxnoncoprime) {
				cout << "you are running out of primes. " << maxnoncoprime << " coprime primes found";
			        return false;
				//return Builder_.result(res);
			}
		}
		coprime = 0;
		Domain D(*genprime); 
                Vect<DomainElement, Alloc<DomainElement>  > r; 
                Builder_.progress( D, Iteration(r, D) );
		//if (RR_.scheduled(IterCounter-1)) {
			++IterCounter;
/*
			Integer M ; Builder_.getModulus(M);
			if ( Builder_.terminated() ) {
				Vect<Integer> r ; Builder_.getResidue(r);
				if (RR_.reconstructRational(num,den,r,M) ) {
					Vect<Integer> vnum(num),vden(m_in.size(),den);
					Builder_.productin(vnum, f_in); Builder_.productin(vden,m_in);
					if (Builder_.changePreconditioner(vnum,vden)) break;
					else Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
				}
			} else {
				Integer r ; Builder_.getResidue(r);
			        Integer n,d;
			        if (RR_.reconstructRational(n,d,r,M)) {
				        Vect<Integer > vden(m_in.size(),d);
				        Builder_.productin(vden,m_in);
				        if (Builder_.changePreconditioner(f_in,vden)) break;
					else Builder_.changePreconditioner(f_in,m_in);//not optimal, may restore results
				}
			}
		}
*/  
  	}
            if (Builder_.terminated()) {
		//Vect<Integer> p_mul, p_div,res, g;
                //Builder_.getPreconditioner(p_div,p_mul);
                Builder_.result(num,den);
                //num = Builder_productin(res,p_div);
                //typename Vect<Integer>::iterator itnum,itden,ittmp;
                //den = 1; Integer denold = 1;
                //for (itnum = num.begin(), itden=p_mul.begin(); itnum != num.end(); ++itnum,++itden) {
                //        lcm(den,den,*itden);
                //        if (denold != den) {
                //                Integer h = den/denold;
                //                ittmp = num.begin();
                //                for (;itnum != itnum; ++ittmp)  *ittmp *= h;
		//	  }
                //        denold = den;
		//}

                //den = p_mul;
	    	return true;
	    } else return false;
        }

#ifdef CRATIMING
        std::ostream& reportTimes(std::ostream& os) {
        	//Builder_.reportTimes(os);
		return os <<  "Iterations:" << IterCounter << "\n" ;
	}
#endif
};

#ifdef RCRATIMING

class RCRATimer {
	public: 
		mutable Timer ttInit, tt RRecon, ttIRecon, ttImaging, ttIteration, ttOther;
		void clear() const {
			ttInit.clear();
			ttRRecon.clear();
			ttIRecon.clear();
			ttImaging.clear();
			ttIteration.clear();
			ttother.clear();
		}
		/*
		template<class RR, class LC>
		void update(RR& rr, LC& lc) const {
			ttSetup += lc.ttSetup;
			ttRecon += rr.ttRecon;
			ttGetDigit += lc.ttGetDigit;
			ttGetDigitConvert += lc.ttGetDigitConvert;
			ttRingOther += lc.ttRingOther;
			ttRingApply += lc.ttRingApply;
		}
		*/
	};
#endif

}

#undef RCRATIMING
#undef CRATIMING

#endif // __LINBOX_rational2_cra_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
