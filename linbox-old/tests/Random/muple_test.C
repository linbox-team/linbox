// ===================================================================
// (C) The Linbox Group 2000
// M-uple test for random generators
// Uses PRNG library and the associated wrapper
// Time-stamp: <20 Mar 00 17:50:55 Jean-Guillaume.Dumas@imag.fr> 
// ===================================================================
#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include "lin_rand.h"
#include "lin_rand48.h"
#include "lin_rand_prng.h"
#include <set.h>
#include <vector.h>

#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif
long nb = 10000, vect=1; 
long size = 16, stop = size;
double Space;

template< class RG > void test_m_uples(RG generator) {
     for(unsigned long div=1,pow=2; div <= stop; ++div, pow*=2)
        if ((size/div)*div == size) {
            set< vector<unsigned long> > S; 
            for(unsigned long j=nb;j--;) {
                generator.setseed();
                vector<unsigned long> u(vect);
                for(unsigned long k =vect; k--;) {
                    u[k] = generator() % (pow);
                    for(unsigned long i=size/div; --i;) 
                        u[k] = (pow)*u[k]+(generator() % (pow) );
                }
                S.insert(u);
            }
            cerr << "[ " << div << " : " 
                 << (S.size()*100.0)/Space << " % ] ";
        }
}
    
// ---------------------------------------------
// argv[1] : Number of different seeds tested
// argv[2] : Number of uples
// argv[3] : Number of bits per number
// argv[4] : bits are generated 1 by 1, 2 by 2, .., argv[4] by argv[4]
//
// Size of the number space : (2^argv[3])^argv[2]
// Shows percentages as the number of distinct generated values
// among all possibilities.
// rand48 is very bad with this test for small bits.  

int main(int argc, char* argv[]) {
    if (argc > 1)
        nb = atoi( argv[1] );
    if (argc > 2)
        vect = atoi( argv[2] );
    if (argc > 3)
        size = atoi( argv[3] );
    if (argc > 4)
        stop = atoi( argv[4] );

    unsigned long p = 1; 
    for(unsigned long j=0;(j<vect) && (p<nb); ++j)
        for(unsigned long i=0; i<size; ++i) p*=2;
    Space = (double)(GIVMIN(nb,p));
    
    { cerr << "r48 : ";
    test_m_uples(Random48());
    cerr << endl; }
    
    { cerr << "lin : ";
    test_m_uples(Random());
    cerr << endl; }
    
    { char * method = new char[30+12];
    sprintf(method,"LCG(2147483647,1078318381,0,%ld)",Random::seed());
    cerr << "lcg : ";
    test_m_uples(PRNG_Random(method));
    cerr << endl; }
    
    { cerr << "tt8 : ";
    test_m_uples(PRNG_Random("tt800"));
    cerr << endl; }
    
    { cerr << "prn : ";
    test_m_uples(PRNG_Random());
    cerr << endl; }
    
    { char * method = new char[30+12];
    sprintf(method,"ICG(2147483647,1,1,%ld)",Random::seed());
    cerr << "icg : ";
    test_m_uples(PRNG_Random(method));
    cerr << endl; }
    
    { char * method = new char[30+12];
    sprintf(method,"DICG(30,17928205,6,%ld)",Random::seed());
    cerr << "dic : ";
    test_m_uples(PRNG_Random(method));
    cerr << endl; }
    
    
    

    return 0;
}
