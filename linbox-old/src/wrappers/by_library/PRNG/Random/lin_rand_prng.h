// ===================================================================
// (C) The Linbox Group 2000
// Random Generator wrapper for PNRG library
// http://random.mat.sbg.ac.at/ftp/pub/software/gen
// Time-stamp: <17 Mar 00 12:43:19 Jean-Guillaume.Dumas@imag.fr> 
// ===================================================================

#ifndef _RAN_PRNG_WRAP_
#define _RAN_PRNG_WRAP_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "lin_rand.h"

extern "C" {
#line __LINE__ "PRNG library is needed in lin_rand_prng.h"
#include "prng.h"
}

class PRNG_Random {
    struct prng *_g;
public:

    PRNG_Random(char * method = 0) {
        if (method == 0) {
            method = new char[30+12];
            sprintf(method,"EICG(2147483647,7,0,%ld)",Random::seed());
        }
        _g = prng_new(method); 
        prng_get_next(_g); 
        if (prng_is_congruential(_g)) prng_get_next_int(_g);
    }

    ~PRNG_Random() {
       prng_free(_g);
    }
            
    unsigned short precision() { return 32; }

    unsigned long setseed(const unsigned long s = 0) {
        unsigned long ns = ((s == 0)?Random::seed():s);
        if (prng_can_seed(_g)) prng_seed(_g,ns);
        return ns;
    }    

    unsigned long operator() () {
        long r; this->operator()(r);
        return (r<0)?(unsigned long)(r-INT32_MIN):(unsigned long)((unsigned long)(r)-INT32_MIN);
    }
    
    unsigned long& operator() (unsigned long& r) {
        return r = this->operator()();
    }
    
    long& operator() (long& l) {
        if (prng_is_congruential(_g))
            return l = prng_get_next_int(_g);
        else {
            unsigned long r, s; double d;
            this->operator()(d);
            r = (unsigned long)(this->operator()(d) * UINT32_MAX) ;
            return l = ( (r>(unsigned long)INT32_MAX) ? (long)(r+INT32_MIN): ((long)(r)+INT32_MIN) );
        }
    }
    
    
    double& operator() (double& d) {
        return d = prng_get_next(_g);
    }    

        /// Default implementation
    template< class XXX>
    XXX& operator() (XXX& x) {
        double d;
        return x = this->operator()(d);
    }        
    
};

#endif
