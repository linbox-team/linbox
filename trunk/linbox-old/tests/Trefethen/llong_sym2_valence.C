// ======================================================================= //
// Wiedemann algorithm using Massey
// With diagonal Scaling and Transpose Computation
// Time-stamp: <28 Jan 02 18:07:41 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------
#include <givtimer.h>
#include "givinteger.h"
#include <givintfactor.h>
#include "givintrns.h"
#include <givinit.h>
#include "givzpz64std.C"
typedef ZpzDom<Std64>               GFqDomain;
typedef GFqDomain::Residu_t        Residu;
   
#include <givpoly1.h>
typedef Poly1Dom< GFqDomain, Dense > Polys;

#include "LinBox/lin_cassini.h"                   // Cassini bound

#include "LinBox/lin_rand.h"                      // Random Iterator
#include "lin_spv_modbb.h"                    // BB Wrapper for sparse vectors
#include "LinBox/lin_symmetric_bbit.h"           // BB iterator
#include "lin_llmassey.C"                // massey reccuring sequ536870912ence solver

typedef SparseBlackBoxModularDom< GFqDomain > SPBB ;
typedef BB_Symmetric_Container< SPBB > SzCBB;

#define LOWERPRIME 1073741824ULL

struct New_Prime { 
    void operator() (Random& generator, vector<unsigned long long>& vp) {
    srand48(BaseTimer::seed());
    srandom(BaseTimer::seed());
        IntPrimeDom I;
        IntPrimeDom::element P;
        vector<unsigned long long>::const_iterator vi;
        if (vp.size())
            do {
//	                I.random(generator, P, IntPrimeDom::element( LOWERPRIME ) );toto.getvalue(),
		unsigned long long ss; generator(ss);
                I.nextprime(P, IntPrimeDom::element(LOWERPRIME+(ss%LOWERPRIME)));
                for (vi = vp.begin(); vi != vp.end(); ++vi)
                    if ( I.areEqual(*vi, P) ) break;
            } while ( vi != vp.end() ) ;
        else {
	    unsigned long long ss; generator(ss);
            I.nextprime(P, IntPrimeDom::element(LOWERPRIME+(ss%LOWERPRIME)));
            //I.random(generator, P, IntPrimeDom::element( LOWERPRIME ) );
            //I.nextprime(P, I.addin(P, IntPrimeDom::element(LOWERPRIME)));
        }
        vp.push_back( I.Integer2long(P) );
    }
};


struct One_Wiedemann {
    void operator()(IntPrimeDom::element& V, unsigned long long& rk,  char * Matrix_File_Name, Random& generator, long long p, unsigned long long Early_Term) {
        GFqDomain F(p);
        SPBB MF( F , Matrix_File_Name );  
        SzCBB TF( &MF, generator); 
//         MasseyDom< SzCBB >  WD(Commentator(PRINT_EVERYTHING,PRINT_EVERYTHING,2,4,EXPO_ESTIMATE), &TF, Early_Term );
        MasseyDom< SzCBB >  WD( &TF, Early_Term );
        GFqDomain::Rep v;
        WD.valence(v,rk);
        F.write(V, v);
    }
};




// ---------------------------------------------
// MAIN
int main(int argc, char* argv[]) {
    // Givaro::Init(&argc, &argv);
        // Must have a matrix name    
    Timer init; init.clear();init.start();
    char* Matrix_File_Name = argv[1];

    unsigned long long valBOUND = 10;
    if (argc > 2) valBOUND = atoi( argv[2] );
    

    Cassini<double>::Type_t OCB;
    {
        Cassini<double>(Matrix_File_Name).operator() (OCB);
        init.stop();
        cerr << "Cassini bound : " << OCB << endl << init << endl;
    }

    Timer tim; tim.clear();tim.start();
        // Default seed is timeofday.
    Random generator;

    
   unsigned long long Early_Term = 20;
    if (argc > 2)
        Early_Term = atoi( argv[2] );


    vector<unsigned long long> p;
    New_Prime() (generator, p);

    vector<IntPrimeDom::element> V(1);
    vector<unsigned long long> rk(1);

    One_Wiedemann() (V[0], rk[0], Matrix_File_Name, generator,p[0], Early_Term);
    if (rk[0] & 0x1) {
        if (V[0] > 0)
            V[0] -= Integer(p[0]);
    } else
        if (V[0] < 0)
            V[0] +=  Integer(p[0]);

    

    unsigned long long NBLAUNCH = (unsigned long long)ceil( (double)rk[0] * log(2.0*OCB) / log(LOWERPRIME) );

    cerr << "nb loops : " << NBLAUNCH << endl;
    cerr << "P[0]: " << p[0] 
         << ", V[0]: " << V[0] 
         << ", rk[0]: " << rk[0] << endl;

    IntRNSsystem<vector>::element Valence;
    unsigned long long diff_bound = 0;

    for(unsigned long long i=1;i<NBLAUNCH;++i) {
        V.resize(V.size()+1);
        rk.resize(rk.size()+1);
        New_Prime() (generator, p);
        One_Wiedemann() (V[i], rk[i], Matrix_File_Name, generator, p[i], Early_Term);
        cerr << "P[" << i << "]: " << p[i] 
             << ", V[" << i << "]: " << V[i] 
             << ", rk[" << i << "]: " << rk[i] << endl;
   
        IntRNSsystem<vector> RNs(p);
        IntRNSsystem<vector>::element V0;
        RNs.RnsToRing(V0, V);
        if (rk[i] & 0x1) {
            if (V0 > 0)
                V0 -= RNs.product();
        } else
            if (V0 < 0)
                V0 += RNs.product();
        if ((V0 == Valence) && (valBOUND > 0) ) if (++diff_bound > valBOUND) { cerr << "Early Termination" << endl; break; }
        Valence = V0;
    }
    
    tim.stop();
    cerr << "Valence: " << Valence << endl << init+tim << endl;

//    Givaro::End();
    return 0;
};

