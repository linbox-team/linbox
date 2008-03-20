// ======================================================================= //
// Time-stamp: <09 Mar 07 18:45:55 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_CRA_H
#define __LINBOX_CRA_H
#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>

namespace LinBox {
    
    template<class CRABase>
    struct ChineseRemainder {
        typedef typename CRABase::Domain	Domain;
        typedef typename CRABase::DomainElement	DomainElement;
    protected:
        CRABase Builder_;
        
    public:
        template<class Param>
        ChineseRemainder(const Param& b) : Builder_(b) {}

            /** \brief The CRA loop
             *
             * Given a function to generate residues mod a single prime, 
             * this loop produces the residues resulting from the Chinese 
             * remainder process on sufficiently many primes to meet the 
             * termination condition.
             *
             * \param Iteration - Function object of two arguments, F(r, p), 
             * given prime p it outputs residue(s) r. This loop may be 
             * parallelized.  F must be reentrant, thread safe. For example, 
             * F may be returning the coefficients of the minimal polynomial 
             * of a matrix mod p.
             *
             * Warning - we won't detect bad primes.
             *
             * \param PrimeIterator - iterator for generating primes.
             * 
             * \result res - an integer
             */
        template<class Int, class Function, class PrimeIterator>
        Int& operator() (Int& res, Function& Iteration, PrimeIterator& primeiter) {
            ++primeiter; 
            Domain D(*primeiter); 
            DomainElement r; D.init(r);
            Builder_.initialize( D, Iteration(r, D) );

            while( ! Builder_.terminated() ) {
                ++primeiter; while(Builder_.noncoprime(*primeiter) ) ++primeiter; 
                Domain D(*primeiter); 
                DomainElement r; D.init(r);
                Builder_.progress( D, Iteration(r, D) );
            }
            return Builder_.result(res);
        }

        
      template<class Int, template <class, class> class Vect, template <class> class Alloc, class Function, class PrimeIterator>
	  Vect<Int,Alloc<Int> > & operator() (Vect<Int,Alloc<Int> >& res, Function& Iteration, PrimeIterator& primeiter) {
            
            ++primeiter; 
            Domain D(*primeiter); 
            Vect<DomainElement, Alloc<DomainElement> > r; 
            Builder_.initialize( D, Iteration(r, D) );

            while( ! Builder_.terminated() ) {
                ++primeiter; while(Builder_.noncoprime(*primeiter) ) ++primeiter; 
                Domain D(*primeiter); 
                Vect<DomainElement, Alloc<DomainElement> > r; 
                Builder_.progress( D, Iteration(r, D) );
            }
            return Builder_.result(res);
        }

    };

}

#endif
