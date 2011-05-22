// ======================================================================= //
// Time-stamp: <15 Mar 07 17:04:48 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#ifndef __LINBOX_CRA_EARLY_SINGLE_H
#define __LINBOX_CRA_EARLY_SINGLE_H

#include "linbox/util/timer.h"
#include <stdlib.h>
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include <vector>
#include <utility>

namespace LinBox {
    
    template<class Domain_Type>
    struct EarlySingleCRA {
        typedef Domain_Type			Domain;
        typedef typename Domain::Element DomainElement;
        typedef EarlySingleCRA<Domain> Self_t;

    protected:
            // PrimeProd*nextM_ is the modulus
        Integer 	primeProd_;
        Integer		nextM_;
        Integer 	residue_; 	// remainder to be reconstructed
        unsigned int    occurency_;	// number of equalities
        
        const unsigned int    EARLY_TERM_THRESHOLD;
      

    public:

        EarlySingleCRA(const unsigned long EARLY=DEFAULT_EARLY_TERM_THRESHOLD)
                : primeProd_(1UL), nextM_(1UL), occurency_(0), EARLY_TERM_THRESHOLD(EARLY-1) {
        }

//         EarlySingleCRA(const Self_t& c) 
//                 : primeProd_(c.primeProd_),
//                   nextM_(c.nextM_),
//                   residue_(c.residue_),
//                   occurency_(c.occurency_), 
//                   EARLY_TERM_THRESHOLD(c.EARLY_TERM_THRESHOLD) {}
        

        virtual void initialize (const Domain& D, const DomainElement& e) {
            D.characteristic( primeProd_ );
            nextM_ = 1UL;
            D.convert( residue_, e);
            occurency_ = 1;
        }

        virtual Integer& result(Integer& d) {
            return d=residue_;
        }
        
        virtual bool terminated() {            
            return occurency_>EARLY_TERM_THRESHOLD;
        }

        virtual void progress (const Domain& D, const DomainElement& e) {
                // Precondition : initialize has been called once before
            primeProd_ *= nextM_; 
            D.characteristic( nextM_ );
            
            DomainElement u0;
            if (D.areEqual( D.init(u0, residue_), e)) {
                ++occurency_;
            } else {
                occurency_ = 1;
                
                D.negin(u0);       	// u0  <-- -u0
                D.addin(u0, e);    	// u0  <-- e-u0
                
                DomainElement m0; 
                D.init(m0, primeProd_);
                D.invin(m0);       	// m0  <-- m0^{-1} mod nextM_
                D.mulin(u0, m0);   	// u0  <-- (e-u0)( m0^{-1} mod nextM_ )
                
                Integer res;
                D.convert(res, u0);	// res <-- (e-u0)( m0^{-1} mod nextM_ )		
                    			// and res < nextM_ 
                
                Integer tmp(res);
                tmp -= nextM_;
                if (absCompare(res,tmp)>0) res = tmp; // Normalize
                
                res *= primeProd_;	// res <-- (e-u0)( m0^{-1} mod nextM_ ) m0	
                    			// and res <= (m0.nextM_-m0)
                
                residue_ += res;	// <-- u0 + (e-u0)( m0^{-1} mod nextM_ ) m0
					// and res <  m0.nextM_
            }
        }

        virtual bool noncoprime(const Integer& i) const {
            Integer g;
            return ( (gcd(g, i, nextM_) != 1) || (gcd(g, i, primeProd_) != 1) );
        }


        virtual ~EarlySingleCRA() {}

    };

    
}
#endif

