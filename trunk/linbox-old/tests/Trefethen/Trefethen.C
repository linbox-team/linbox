// ========================================================== //
// Trefethen problem
// Time-stamp: <27 Jan 02 19:35:32 Jean-Guillaume.Dumas@imag.fr>// ========================================================== //
 
#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------
#include <givtimer.h>
#include <givinit.h>
#include "givintprime.h"

#define ABS(a) ((a)>0?(a):-(a))

int main(int argc, char ** argv) {
    
    IntPrimeDom IPD; Integer prime = 2;
    
    unsigned long size = 20000;
    unsigned long start = 1;
    if (argc > 1) size = atoi( argv[1] );
    if (argc > 2) start = atoi( argv[2] );
    
    
    Timer tim; tim.clear(); tim.start();
    
    printf("%ld %ld M\n",size-start+1,size-start+1);
    for (unsigned long i = 1; i<=size; ++i) {
        for (unsigned long j = 1; j<=size; ++j) {
            long diff = i-j;
            unsigned long imj = ABS(diff);
            switch (imj) {
                case 0:
                    if ( (i>=start) && (j>=start)) printf("%ld %ld %ld\n",i-start+1,j-start+1,(long)prime);
                    break;

                case 1:  case 2: case 4:   case 8:   case 16:   case 32:   case 64:   case 128:   case 256:   case 512:   case 1024:   case 2048:   case 4096:   case 8192:   case 16384:
                    if ( (i>=start) && (j>=start)) printf("%ld %ld 1\n",i-start+1,j-start+1);
                    break;
            }
        }
        IPD.nextprimein(prime);
    }


    printf("0 0 0\n");

    tim.stop();
    cerr << tim << endl;
    return 0;
    
}

	
