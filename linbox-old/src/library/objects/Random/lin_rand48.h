// ===================================================================
// (C) The Linbox Group 2000
// Random Generator
// Time-stamp: <20 Mar 00 17:53:12 Jean-Guillaume.Dumas@imag.fr> 
// ===================================================================
#ifndef _LIN_RANDOM_48_H_
#define _LIN_RANDOM_48_H_


#include <lin_rand.h>

class Random48 {
public:

    Random48(const unsigned long s=0) {
        if (s)  srand48( s ); else srand48(Random::seed()); 
    }

    unsigned long setseed(const unsigned long s = 0) {
        if (s)  srand48( s ); else srand48(Random::seed()); 
    }

    unsigned short precision() {
        return 32;
    }

    unsigned long operator() () {
        return ((unsigned long) (INT32_MAX) - lrand48());
    }    

    unsigned long& operator() (unsigned long& r) {
        return r = this->operator() ();
    }    
    
    long& operator() (long& r) {
        return r = mrand48();
    }    
    
    double& operator() (double& d) {
        return d = drand48();
    }    

        /// Default implementation for any type.
    template<class XXX>
    XXX& operator() (XXX& x) {
        return x = this->operator() ();
    }                        
    
};

#endif
