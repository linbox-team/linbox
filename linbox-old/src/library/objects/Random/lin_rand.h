// ===================================================================
// (C) The Linbox Group 2000
// Random Generator
// Time-stamp: <10 Apr 00 17:45:06 Jean-Guillaume.Dumas@imag.fr> 
// ===================================================================
#ifndef _LIN_RANDOM_H_
#define _LIN_RANDOM_H_


#include <math.h>
#include <limits.h>
extern "C" {
# include <sys/time.h>
}

// -----------------------------------------------------
// Fishman, G.S. "Multiplicative congruential random
// number generators ..." Math. Comp. 54:331-344 (1990)
#define _RAN_MULTIPLYER_ 950706376
#define _RAN_MODULO_     2147483647
// -----------------------------------------------------

class Random {
    unsigned long _seed;
    const unsigned long _mul, _mod;
public:

    Random(const unsigned long s=0, const unsigned long m=_RAN_MULTIPLYER_, const unsigned long p=_RAN_MODULO_) 
            : _seed(s), _mul(m), _mod(p) { if (!s) _seed = seed(); }

    Random( const Random& R) : _seed(R._seed), _mul(R._mul), _mod(R._mod) {}
    
    Random& operator= (const Random& R) { _seed = R._seed; }

        // Returns a value to initialize random generator 
    static unsigned long seed() {
        struct timeval tp;
        gettimeofday(&tp, 0) ;
        return (tp.tv_usec);
    }
    
    unsigned long setseed(const unsigned long s = 0) {
        return _seed = (s == 0)?seed():s;
    }

    unsigned short precision() {
        return (short)log10(_mod);
    }

    unsigned long operator() () {
        return _seed = (unsigned long) ( (unsigned long long)_mul * (unsigned long long)(_seed) % (unsigned long long)_mod );
    }    

    unsigned long& operator() (unsigned long& r) {
        return r = this->operator() ();
    }    
    
    long& operator() (long& r) {
        unsigned long u; this->operator()(u);
        return r = ((u>INT_MAX)? (long)(u - INT_MAX ) : (long)(u) - INT_MAX);
    }    
    
    double& operator() (double& d) {
        unsigned long r; this->operator()(r);
        return d = double(r)/double(_mod);
    }    

        /// Default implementation for any type.
    template<class XXX>
    XXX& operator() (XXX& x) {
        return x = this->operator() ();
    }                        
    
};

#endif
