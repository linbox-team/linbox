#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
#include "lin_rand.h"
#include "lin_rand_prng.h"
#include <set.h>
#include <vector.h>
#include <algo.h>

#ifndef GIVMIN
#define GIVMIN(a,b) ((a)<(b)?(a):(b))
#endif
unsigned long n = 16777216, m=512, nb=500; 
double tm;

template< class RG > void test_birthday(RG generator) {
    double mean = 0.0;
    for(unsigned long i=0;i<nb;) {
        set< unsigned long > S; 
        for(unsigned long k =m; k--;) S.insert( generator() % (n) );
        vector< unsigned long > V(S.begin(), S.end());
        sort(V.begin(),V.end());
        set< unsigned long > Gaps; set< unsigned long > G2;
        vector< unsigned long >::const_iterator si=V.begin(), sp=V.begin();
        unsigned long s1, s2 = 0;
        for(++si;si!=V.end();++si, ++sp) {
            s1 = s2;
            Gaps.insert( *si - *sp );
            if ( (s2 = Gaps.size()) == s1)
                G2.insert( *si - *sp );
        }
        mean = (mean * double(i) + G2.size() )/double(++i);
    }
    cerr << "mean = " << mean << " --> " << (mean * 100.0)/tm << " %";
}


    
// ---------------------------------------------

int main(int argc, char* argv[]) {
    if (argc > 1)
        nb = atoi( argv[1] );
    if (argc > 2)
        n = atoi( argv[2] );
    if (argc > 3)
        m = atoi( argv[3] );
    srand48(Random::seed());

    tm = (((unsigned long long)m*(unsigned long long)m*(unsigned long long)m)/((unsigned long long)4*(unsigned long long)n));

    cerr << "Theory mean = " << tm << endl;

    cerr << "r48  : ";
    double mean = 0.0;
    for(unsigned long i=0;i<nb;) {
        set< unsigned long > S; 
        for(unsigned long k =m; k--;) S.insert( lrand48() % (n) );
        vector< unsigned long > V(S.begin(), S.end());
        sort(V.begin(),V.end());
        set< unsigned long > Gaps; set< unsigned long > G2;
        vector< unsigned long >::const_iterator si=V.begin(), sp=V.begin();
        unsigned long s1, s2 = 0;
        for(++si;si!=V.end();++si, ++sp) {
            s1 = s2;
            Gaps.insert( *si - *sp );
            if ( (s2 = Gaps.size()) == s1)
                G2.insert( *si - *sp );
        }
        mean = (mean * double(i) + G2.size() )/double(++i);
    }
    cerr << "mean = " << mean << " --> " << (mean * 100.0)/tm << " %";
    cerr << endl;
            
  

    { cerr << "lin  : ";
    test_birthday(Random());
    cerr << endl; }
            
    { char * method = new char[30+12];
    sprintf(method,"LCG(2147483647,1078318381,0,%ld)",Random::seed());
    cerr << "lcg  : ";
    test_birthday(PRNG_Random(method));
    cerr << endl; }
    
    { cerr << "prn  : ";
    test_birthday(PRNG_Random());
    cerr << endl; }
    
    { char * method = new char[30+12];
    sprintf(method,"ICG(2147483647,1,1,%ld)",Random::seed());
    cerr << "icg  : ";
    test_birthday(PRNG_Random(method));
    cerr << endl; }
    
    { cerr << "tt80 : ";
    test_birthday(PRNG_Random("tt800"));
    cerr << endl; }
            
    { cerr << "cmrg : ";
    test_birthday(PRNG_Random("cmrg"));
    cerr << endl; }
            
    { cerr << "mrg  : ";
    test_birthday(PRNG_Random("mrg"));
    cerr << endl; }
            
    { cerr << "ctg  : ";
    test_birthday(PRNG_Random("ctg"));
    cerr << endl; }
            
    { char * method = new char[30+12];
    sprintf(method,"DICG(30,17928205,6,%ld)",Random::seed());
    cerr << "dic  : ";
    test_birthday(PRNG_Random(method));
    cerr << endl; }
    
   

    return 0;
}
