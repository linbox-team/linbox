// ======================================================================= //
// Wiedemann algorithm using Massey
// With diagonal Scaling and Transpose Computation
// Time-stamp: <06 Apr 00 16:46:28 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------
#include <givtimer.h>
#include "givinteger.h"
#include <givintfactor.h>
#include "givintrns.h"
#include <givinit.h>
#include "givgfq.h"
typedef GFqDom<long>               GFqDomain;
typedef GFqDomain::Residu_t        Residu;
   
#include <givpoly1.h>
typedef Poly1Dom< GFqDomain, Dense > Polys;

#include "LinBox/lin_cassini.h"                   // Cassini bound

#include "LinBox/lin_rand.h"                      // Random Iterator
#include "LinBox/lin_spv_bb.h"                    // BB Wrapper for sparse vectors
#include "LinBox/lin_symmetrize_bbit.h"           // BB iterator
#include "LinBox/lin_massey.C"                // massey reccuring sequence solver

typedef SparseBlackBoxDom< GFqDomain > SPBB ;
typedef BB_Symmetrize_Container< SPBB > SzCBB;

#define LOWERPRIME 32768    

struct New_Prime { 
    void operator() (Random& generator, vector<unsigned long>& vp) {
        IntPrimeDom I;
        IntPrimeDom::element P;
        vector<unsigned long>::const_iterator vi;
        if (vp.size())
            do {
                I.random(generator, P, IntPrimeDom::element( LOWERPRIME ) );
                I.nextprime(P, I.addin(P, IntPrimeDom::element(LOWERPRIME)));
                for (vi = vp.begin(); vi != vp.end(); ++vi)
                    if ( I.areEqual(*vi, P) ) break;
            } while ( vi != vp.end() ) ;
        else {
            I.random(generator, P, IntPrimeDom::element( LOWERPRIME ) );
            I.nextprime(P, I.addin(P, IntPrimeDom::element(LOWERPRIME)));
        }
        vp.push_back( I.Integer2long(P) );
    }
};


struct One_Wiedemann {
    void operator()(IntPrimeDom::element& V, unsigned long& rk,  char * Matrix_File_Name, Random& generator, long p, unsigned long Early_Term) {
        GFqDomain F(p,1);
        SPBB MF( F, Matrix_File_Name );  
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

    Cassini<double>::Type_t OCB;
    {
        Cassini<double>(Matrix_File_Name).operator() (OCB);
        init.stop();
        cerr << "Cassini bound : " << OCB << endl << init << endl;
    }

    Timer tim; tim.clear();tim.start();
        // Default seed is timeofday.
    Random generator;

    
   unsigned long Early_Term = 20;
    if (argc > 2)
        Early_Term = atoi( argv[2] );


    vector<unsigned long> p;
    New_Prime() (generator, p);

    vector<IntPrimeDom::element> V(1);
    vector<unsigned long> rk(1);

    One_Wiedemann() (V[0], rk[0], Matrix_File_Name, generator,p[0], Early_Term);

    unsigned long NBLAUNCH = (unsigned long)ceil( (double)rk[0] * log(2.0*OCB) / log(LOWERPRIME) );

    cerr << "nb loops : " << NBLAUNCH << endl;
    cerr << "P[0]: " << p[0] 
         << ", V[0]: " << V[0] 
         << ", rk[0]: " << rk[0] << endl;
   
    V.resize(NBLAUNCH);
    rk.resize(NBLAUNCH);

    for(unsigned long i=1;i<NBLAUNCH;++i) {
        New_Prime() (generator, p);
        One_Wiedemann() (V[i], rk[i], Matrix_File_Name, generator, p[i], Early_Term);
        cerr << "P[" << i << "]: " << p[i] 
             << ", V[" << i << "]: " << V[i] 
             << ", rk[" << i << "]: " << rk[i] << endl;
    }
   
    IntRNSsystem<vector> RNs(p);
    IntRNSsystem<vector>::element Valence;
    RNs.RnsToRing(Valence, V);
    if (rk[0] & 0x1) {
        if (Valence > 0)
            Valence = Valence - RNs.product();
    } else
        if (Valence < 0)
            Valence = Valence + RNs.product();
    
    tim.stop();
    cerr << "Valence: " << Valence << endl << init+tim << endl;

//    Givaro::End();
    return 0;
};

