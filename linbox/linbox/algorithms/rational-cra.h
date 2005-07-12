/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// ======================================================================= //
// Time-stamp: <12 Jul 05 11:45:46 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_RATIONAL_CRA_H
#define __LINBOX_RATIONAL_CRA_H

#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/cra-domain.h"

namespace LinBox {

template<class T, template <class T> class Container>
std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
          for(typename Container<T>::const_iterator refs =  C.begin();
                                refs != C.end() ;
                                      ++refs )
                          o << (*refs) << " " ;
            return o << std::endl;
}


	/** \brief Chinese remainder of rationals
	 *
	 * Compute the reconstruction of rational numbers
	 * Either by Early Termination see [Dumas, Saunder, Villard, JSC 32 (1/2), pp 71-99, 2001],
	 * Or via a bound on the size of the integers.
	 */
    template<class Domain>
    struct RationalRemainder : public ChineseRemainder<Domain> {
		typedef ChineseRemainder<Domain> Father_t;
		PID_integer _ZZ;
    public:

        RationalRemainder(const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD, const size_t n=1) 
				: Father_t(EARLY, n) {}
        
        RationalRemainder(const double BOUND) 
				: Father_t(BOUND) {}

		
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
        template<class Function, class RandPrime>
        Integer & operator() (Integer& num, Integer& den, const Function& Iteration, RandPrime& genprime) {
            Integer p;
                while( ! this->Full_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Full_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    typename Father_t::DomainElement r; D.init(r);
                    this->Full_progress( D, Iteration(r, D) );
                }
                return this->Full_result(num, den);
        }

        template<template <class T> class Vect, class Function, class RandPrime>
        Vect<Integer> & operator() (Vect<Integer>& num, Integer& den, const Function& Iteration, RandPrime& genprime) {
            Integer p;
                while( ! this->Full_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Full_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    Vect<typename Father_t::DomainElement> r; 
                    this->Full_progress( D, Iteration(r, D) );
                }
                return this->Full_result(num, den);
        }

      protected:

        Integer& Full_result (Integer &num, Integer& den){
            std::vector<Integer> Vd; Vd.push_back(num);
            Full_result(Vd, den);
            return num=Vd.front();
        }
		
        template<template<class T> class Vect>
        Vect<Integer>& Full_result (Vect<Integer> &num, Integer& den){
            num.resize( (Father_t::Table.front()).size() );
            std::vector< LazyProduct >::iterator 			_mod_it = Father_t::Modulo.begin();
            std::vector< std::vector< Integer > >::iterator _tab_it = Father_t::Table.begin();
            std::vector< bool >::iterator    				_occ_it = Father_t::Occupation.begin();
            LazyProduct Product;
            for( ; _occ_it != Father_t::Occupation.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
                if (*_occ_it) {
                    Product = *_mod_it;
                    std::vector<Integer>::iterator t0_it = num.begin();
                    std::vector<Integer>::iterator t_it = _tab_it->begin();
                    if (++_occ_it == Father_t::Occupation.end()) {
						den = 1;
						Integer s, nd; _ZZ.sqrt(s, _mod_it->operator()());
						for( ; t0_it != num.end(); ++t0_it, ++t_it) {
                            ratrecon(*t0_it = *t_it, nd, den, _mod_it->operator()(), s);
							if (nd > 1) {
								std::vector<Integer>::iterator  t02 = num.begin();
								for( ; t02 != t0_it ; ++t02)
									*t02 *= nd;
								den *= nd;
							}
						}
                        return num;
                    } else {
                        for( ; t0_it != num.end(); ++t0_it, ++t_it)
                            *t0_it  = *t_it;
                        ++_mod_it; ++_tab_it; 
                        break;
                    }
                }
            }
            for( ; _occ_it != Father_t::Occupation.end() ; ++_mod_it, ++_tab_it, ++_occ_it) {
                if (*_occ_it) {
                    std::vector<Integer>::iterator t0_it = num.begin();
                    std::vector<Integer>::const_iterator t_it = _tab_it->begin();
                    for( ; t0_it != num.end(); ++t0_it, ++t_it)
                        normalizesmallbigreconstruct(*t0_it, Product(), *t_it, _mod_it->operator()() );
                    Product.mulin(*_mod_it);
                }
            }
			den = 1;
			Integer s, nd; _ZZ.sqrt(s, Product.operator()());
			std::vector<Integer>::iterator t0_it = num.begin();
			for( ; t0_it != num.end(); ++t0_it) {
				ratrecon(*t0_it, nd, den, Product.operator()(), s);
				if (nd > 1) {
					std::vector<Integer>::iterator  t02 = num.begin();
					for( ; t02 != t0_it ; ++t02)
						*t02 *= nd;
					den *= nd;
				}
			}
            return num;
        }
		
        Integer& ratrecon(Integer& u1, Integer& new_den, const Integer& old_den, const Integer& m1, const Integer& s) {
			Integer a;
			PID_integer::reconstructRational(a, new_den, u1*=old_den, m1, s, s);
			return u1=a;
        }
        
		

    };
}

#endif
