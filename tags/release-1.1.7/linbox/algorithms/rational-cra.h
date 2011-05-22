/* Copyright (C) 2007 LinBox
 * Written by JG Dumas
 *
 *
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


#ifndef __LINBOX_rational_cra_H
#define __LINBOX_rational_cra_H

#include "linbox/field/PID-integer.h"

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

	    template<template <class, class> class Vect, template <class> class Alloc,  class Function, class RandPrimeIterator>
	    Vect<Integer, Alloc<Integer> > & operator() (Vect<Integer, Alloc<Integer> >& num, Integer& den, Function& Iteration, RandPrimeIterator& genprime) {
            ++genprime;
            Domain D(*genprime); 
            Vect<DomainElement, Alloc<DomainElement> > r; 
            Builder_.initialize( D, Iteration(r, D) );				
            while( ! Builder_.terminated() ) {
                ++genprime; while(Builder_.noncoprime(*genprime) ) ++genprime;
                Domain D(*genprime); 
                Vect<DomainElement, Alloc<DomainElement> > r; 
                Builder_.progress( D, Iteration(r, D) );
            }
            return Builder_.result(num, den);
        }
    };
}

#endif //__LINBOX_rational_cra_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
