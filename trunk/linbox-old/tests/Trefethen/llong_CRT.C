// ======================================================================= //
// Wiedemann algorithm using Massey
// With diagonal Scaling and Transpose Computation
// Time-stamp: <30 Jan 02 11:02:01 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------
#include <givtimer.h>
#include "givinteger.h"
#include <givintfactor.h>
#include "givintrns.h"
#include <givinit.h>


// ---------------------------------------------
// MAIN
int main(int argc, char* argv[]) {
    // Givaro::Init(&argc, &argv);

    unsigned long long rk = atoi(argv[1]);
    
    Timer init; init.clear();init.start();

    
    vector<unsigned long long> P;
    vector<long long> V;
    unsigned long long p, r; long long v;

    while(! cin.eof()) {
	cin >> p >> v >> r;
	if (r == rk) {
	   bool isbr = 0;
	   if (P.size())
	        for (vector<unsigned long long>::const_iterator pi = P.begin(); pi != P.end(); ++pi)
	             if ( p == *pi ) { isbr = 1; break; }
	   if (isbr) {
           	cerr << "isbr excluded: " << p << " " << v << " " << r << endl;
		continue;
	   }
	   P.push_back( p );
	   V.push_back( v );
	} else {
           cerr << "r=rk excluded: " << p << " " << v << " " << r << endl;
	}
    }
    init.stop();
    cerr << "P size [" << init << "]: " << P.size() << endl ;
	
    unsigned long s = P.size();
    

    vector< IntRNSsystem<vector>::element > CheckP;
    vector< IntRNSsystem<vector>::element > CheckV;
    unsigned long checks = s/1000; checks = (checks>0? checks : 10);
    for(unsigned long i = 0; i<checks ; ++i) {
        CheckP.push_back( P.back() ); P.pop_back();
        CheckV.push_back( V.back() ); V.pop_back();
    }
	    
    Timer tim; tim.clear();tim.start();

    IntRNSsystem<vector>::element Valence;
    IntRNSsystem<vector> RNs(P);
    RNs.RnsToRing(Valence, V);

    CheckP.push_back( RNs.product() );
    CheckV.push_back( Valence );
    
     if (rk & 0x1) {
        if (Valence > 0)
            Valence -= RNs.product();
    } else {
        if (Valence < 0)
            Valence += RNs.product();
    }
    
    tim.stop();
    cerr << tim << endl;

    Timer check; check.clear(); check.start();
    IntRNSsystem<vector>::element CheckValence;
    IntRNSsystem<vector> CheckRNs(CheckP);
    CheckRNs.RnsToRing(CheckValence, CheckV);
    check.stop();

    if (rk & 0x1) {
        if (CheckValence > 0)
            CheckValence -= CheckRNs.product();
    } else {
        if (CheckValence < 0)
            CheckValence += CheckRNs.product();
    }
    
    cerr << (check+tim+init) << endl;
   
    if (CheckValence == Valence) 
        cout << "Valence checked(" << checks << ") : " << Valence << endl;
    else
        cout << "Current valence : " << CheckValence << endl;

//    Givaro::End();
    return 0;
};

