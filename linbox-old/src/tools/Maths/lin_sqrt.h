#ifndef _LINBOX_SQRT_H_
#define _LINBOX_SQRT_H_
// ==========================================================================
// (c) 2000 The Linbox Group
// Time-stamp: <05 Apr 00 13:44:41 Jean-Guillaume.Dumas@imag.fr>
// Description: square roots
// ==========================================================================


#include <math.h>

unsigned long sqrt(unsigned long a) {
    unsigned long a0 = a, a1 = a0 >> 1;
    while( a1 != a0 ) { cerr << "a0: " << a0 << ", a1: " << a1 << endl;  a0 = a1; a1 = (a1 + a/a1) >> 1; }
    return a0;
}

unsigned long long sqrt(unsigned long long a) {
    unsigned long long a0 = a, a1 = a0 >> 1;
    while( a1 != a0 ) { cerr << "a0: " << a0 << ", a1: " << a1 << endl; a0 = a1; a1 = (a1 + a/a1) >> 1; }
    return a0;
}

#ifndef LINABS
#define LINABS(a) ((a)<0?-(a):(a))
#endif

long sqrt(long a) { return sqrt( (unsigned long)LINABS(a) ); }
long long sqrt(long long a) { return sqrt( (unsigned long long)LINABS(a) ); }

int sqrt(int a) { return sqrt( (unsigned long)LINABS(a) ); }
unsigned int sqrt(unsigned int a) { return sqrt( (unsigned long)a ); }

#endif
    
