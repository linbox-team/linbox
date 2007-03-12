/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// ======================================================================= //
// Time-stamp: <12 Mar 07 19:40:17 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_RATIONAL_CRA_H
#define __LINBOX_RATIONAL_CRA_H

#include "linbox/field/PID-integer.h"

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
    template<class RatCRABase>
    struct RationalRemainder {
        typedef typename RatCRABase::Domain		Domain;
        typedef typename RatCRABase::DomainElement	DomainElement;
    protected:
        RatCRABase Builder_;
        
    public:
        template<class Param>
        RationalRemainder(const Param& b) : Builder_(b) { }

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
            while( ! Builder_.terminated() ) {
                ++genprime; while(Builder_.noncoprime(*genprime) ) ++genprime;
                Domain D(*genprime); 
                DomainElement r; D.init(r);
                Builder_.progress( D, Iteration(r, D) );
            }
            return Builder_.result(num, den);
        }

        template<template <class T> class Vect, class Function, class RandPrimeIterator>
        Vect<Integer> & operator() (Vect<Integer>& num, Integer& den, Function& Iteration, RandPrimeIterator& genprime) {
            ++genprime;
            Domain D(*genprime); 
            Vect<DomainElement> r; 
            Builder_.initialize( D, Iteration(r, D) );				
            while( ! Builder_.terminated() ) {
                ++genprime; while(Builder_.noncoprime(*genprime) ) ++genprime;
                Domain D(*genprime); 
                Vect<DomainElement> r; 
                Builder_.progress( D, Iteration(r, D) );
            }
            return Builder_.result(num, den);
        }
    };
}

#endif
