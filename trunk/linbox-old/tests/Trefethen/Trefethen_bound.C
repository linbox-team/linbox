// ========================================================== //
// Trefethen problem
// Time-stamp: <27 Jan 02 19:55:57 Jean-Guillaume.Dumas@imag.fr>// ========================================================== //
 
#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------
#include <givtimer.h>
#include <givinit.h>
#include "givintprime.h"

#define ABS(a) ((a)>0?(a):-(a))

int main(int argc, char ** argv) {
    
    IntPrimeDom IPD; Integer prime = 2, prevprime = 2;
    Integer bound = 1;
    
    unsigned radius = 0, prevradius = 0;
    unsigned long size = 20000;
    unsigned long start = 1;
    if (argc > 1) size = atoi( argv[1] );
    if (argc > 2) start = atoi( argv[2] );
    
    
    unsigned long expo = 0, totalexpo = 0;
    Timer tim; tim.clear(); tim.start();
    
    for (unsigned long i = 1; i<=size; ++i) {
        prevprime = prime;
        prevradius = radius;
        ++expo;
        IPD.nextprimein(prime);
        
        radius = 0;
        for( unsigned long j = 1; j<size; j*= 2) {
            if ( (i-j) > 0 ) ++radius;
            if ( (i+j) < size ) ++radius;
        }

        if ((prime-radius) > (prevprime+prevradius) ) {
            bound *= pow(prevprime+prevradius,expo);
            totalexpo += expo;
            cerr << prevprime << "^" << expo << " * ";
            expo = 0;
        }
    }

    bound *= pow(prevprime,expo);
    totalexpo += expo;
    cerr << prevprime << "^" << expo << ";" << endl;
    tim.stop();
    cout << bound << "; evalf(log[10](%);" << endl;
    cerr << logp(bound,10) << endl;
    
    cerr << "total expo: " << totalexpo << " [" << tim << "]" << endl;
    return 0;
    
}

	
