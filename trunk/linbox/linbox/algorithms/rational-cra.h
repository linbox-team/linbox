/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// ======================================================================= //
// Time-stamp: <15 Jan 07 14:16:54 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_RATIONAL_CRA_H
#define __LINBOX_RATIONAL_CRA_H

#include "linbox/field/PID-integer.h"
#include "linbox/algorithms/cra-domain.h"

namespace LinBox {

//     template<class T, template <class T> class Container>
//     std::ostream& operator<< (std::ostream& o, const Container<T>& C) {
//         for(typename Container<T>::const_iterator refs =  C.begin();
//             refs != C.end() ;
//             ++refs )
//             o << (*refs) << " " ;
//         return o << std::endl;
//     }


	/** \brief Chinese remainder of rationals
	 *
	 * Compute the reconstruction of rational numbers
	 * Either by Early Termination see [Dumas, Saunder, Villard, JSC 32 (1/2), pp 71-99, 2001],
	 * Or via a bound on the size of the integers.
	 */
    template<class Domain>
    struct RationalRemainder : public ChineseRemainder<Domain> {
        typedef ChineseRemainder<Domain> Father_t;
        typedef typename Father_t::DomainElement DomainElement;
        PID_integer _ZZ;
    public:

        using Father_t::Table; 
        using Father_t::Table0;
        using Father_t::nextm;
        using Father_t::Modulo0;
        using Father_t::occurency;
        using Father_t::EARLY_TERM_THRESHOLD;

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
        Integer & operator() (Integer& num, Integer& den, Function& Iteration, RandPrime& genprime) {
            Integer p;
            while( ! this->Full_terminated() ) {
                genprime.randomPrime(p);
                while(this->Full_noncoprimality(p) )
                    genprime.randomPrime(p);
                Domain D(p); 
                DomainElement r; D.init(r);
                this->Full_progress( D, Iteration(r, D) );
            }
            return this->Full_result(num, den);
        }

        template<template <class T> class Vect, class Function, class RandPrime>
        Vect<Integer> & operator() (Vect<Integer>& num, Integer& den, Function& Iteration, RandPrime& genprime) {
            Integer p;
            if (EARLY_TERM_THRESHOLD) {
                {
                    genprime.randomPrime(p);
                    Domain D(p); 
                    Vect<DomainElement> r; 
                    Father_t::First_Early_progress( D, Iteration(r, D) );				
                }
                while( ! this->Early_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Early_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    Vect<DomainElement> r; 
                    Father_t::Early_progress( D, Iteration(r, D) );
                }
                return this->Early_result(num, den);
            } else {
                while( ! this->Full_terminated() ) {
                    genprime.randomPrime(p);
                    while(this->Full_noncoprimality(p) )
                        genprime.randomPrime(p);
                    Domain D(p); 
                    Vect<DomainElement> r; 
                    this->Full_progress( D, Iteration(r, D) );
                }
                return this->Full_result(num, den);
            }
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
                            iterativeratrecon(*t0_it = *t_it, nd, den, _mod_it->operator()(), s);
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
                        this->normalizesmallbigreconstruct(*t0_it, Product(), *t_it, _mod_it->operator()() );
                    Product.mulin(*_mod_it);
                }
            }
            den = 1;
            Integer s, nd; _ZZ.sqrt(s, Product.operator()());
            std::vector<Integer>::iterator t0_it = num.begin();
            for( ; t0_it != num.end(); ++t0_it) {
                iterativeratrecon(*t0_it, nd, den, Product.operator()(), s);
                if (nd > 1) {
                    std::vector<Integer>::iterator  t02 = num.begin();
                    for( ; t02 != t0_it ; ++t02)
                        *t02 *= nd;
                    den *= nd;
                }
            }
            return num;
        }
		
        Integer& iterativeratrecon(Integer& u1, Integer& new_den, const Integer& old_den, const Integer& m1, const Integer& s) {
            Integer a;
            _ZZ.reconstructRational(a, new_den, u1*=old_den, m1, s);
            return u1=a;
        }

        void Early_progress (const Domain& D, const DomainElement& e) {
            DomainElement u0, m0;

            fieldreconstruct(Table0, D, e, D.init(u0,Table0), D.init(m0,Modulo0), Integer(Table0), Modulo0);
            D.characteristic( nextm );
            Modulo0 *= nextm;
            Integer a, b;
            _ZZ.reconstructRational(a, b, Table0, Modulo0);
            if ((a == Numer0) && (b == Denom0))
                ++occurency;
            else {
                occurency = 0;
                Numer0 = a;
                Denom0 = b;
            } 
        }

        void First_Early_progress (const Domain& D, const DomainElement& e) {
            D.characteristic( Modulo0 );
            D.convert(Table0, e);
            _ZZ.reconstructRational(Numer0, Denom0, Table0, Modulo0);
        }

        template<template<class T> class Vect>
        Vect<Integer>& Early_result(Vect<Integer>& num, Integer& den) {
            return Full_result(num, den);
        }
 
    protected:
        Integer					Numer0;
        Integer					Denom0;
    };
}

#endif
