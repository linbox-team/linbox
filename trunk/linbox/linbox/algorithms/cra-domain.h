// ======================================================================= //
// Time-stamp: <07 Jul 09 13:35:17 Jean-Guillaume.Dumas@imag.fr> 
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
    
    template<class Function, class Element> struct CRATemporaryVectorTrait {
        typedef std::vector<Element> Type_t;
    };        


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
        template<class Function, class PrimeIterator>
        Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter) {
            commentator.start ("Modular iteration", "mmcrait");
            ++primeiter; 
            Domain D(*primeiter); 
            commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
            DomainElement r; D.init(r);
            Builder_.initialize( D, Iteration(r, D) );
            
            while( ! Builder_.terminated() ) {
                ++primeiter; while(Builder_.noncoprime(*primeiter) ) ++primeiter; 
                Domain D(*primeiter); 
                commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
                DomainElement r; D.init(r);
                Builder_.progress( D, Iteration(r, D) );
            }
            commentator.stop ("done", NULL, "mmcrait");
            return Builder_.result(res);
        }
        
        
      template<class Iterator, class Function, class PrimeIterator>
	  Iterator& operator() (Iterator& res, Function& Iteration, PrimeIterator& primeiter) {
          commentator.start ("Modular vectorized iteration", "mmcravit");
          
          ++primeiter; 
          Domain D(*primeiter); 
          commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
          typename CRATemporaryVectorTrait<Function, DomainElement>::Type_t r; 
          Builder_.initialize( D, Iteration(r, D) );
          
          while( ! Builder_.terminated() ) {
              ++primeiter; while(Builder_.noncoprime(*primeiter) ) ++primeiter; 
              Domain D(*primeiter); 
              commentator.report(Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION) << "With prime " << *primeiter << std::endl;
              typename CRATemporaryVectorTrait<Function, DomainElement>::Type_t r; 
              Builder_.progress( D, Iteration(r, D) );
          }
          commentator.stop ("done", NULL, "mmcravit");
          return Builder_.result(res);
      }
        
    };
    
}

#endif
