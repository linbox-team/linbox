#ifndef _LINBOX_SQRT_H_
#define _LINBOX_SQRT_H_
// ==========================================================================
// (c) 2000 The Linbox Group
// Time-stamp: <05 Apr 00 13:44:41 Jean-Guillaume.Dumas@imag.fr>
// Description: square roots
// ==========================================================================


#include <math.h>

// long sqrt(long a) {
//     return long( sqrt(double(a)) );
// }

// long long sqrt(long long a) {
//     return long long( sqrt(double(a)) );
// }

#ifndef LINABS
#define LINABS(a) ((a)<0?-(a):(a))
#endif

long sqrt(long a) {
    unsigned long a0 = LINABS(a), a1 = a0 >> 1;
    while( a1 != a0 ) { cerr << "a0: " << a0 << ", a1: " << a1 << endl;  a0 = a1; a1 = (a1 + a/a1) >> 1; }
    return a0;
}

long long sqrt(long long a) {
    unsigned long long a0 = LINABS(a), a1 = a0 >> 1;
    while( a1 != a0 ) { cerr << "a0: " << a0 << ", a1: " << a1 << endl; a0 = a1; a1 = (a1 + a/a1) >> 1; }
    return a0;
}

#endif
    
