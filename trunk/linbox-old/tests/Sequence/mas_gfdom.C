// ======================================================================= //
// Wiedemann algorithm using Massey
// With diagonal Scaling and Transpose Computation
// Time-stamp: <05 Apr 00 10:06:13 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------
#include <givtimer.h>
#include <givinit.h>
#include "givgfq.h"
#include <givpoly1.h>

#include "lin_rand.h"                      // Random Iterator
#include "lin_spv_bb.h"                    // BB Wrapper for sparse vectors
#include "lin_symmetrize_bbit.h"           // BB iterator
#include "lin_massey.C"                // massey reccuring sequence solver

// ---------------------------------------------
// MAIN
int main(int argc, char* argv[]) {
    Givaro::Init(&argc, &argv);
//     try{
        
    Timer init; init.clear();init.start();
    
    typedef GFqDom<long>               GFqDomain;
    typedef GFqDomain::Residu_t        Residu;
   
        // Must have a matrix name    
    char* Matrix_File_Name = argv[1];
    Residu MOD=65521,expo=1;
    unsigned long Early_Term = 20;
    if (argc > 2)
        MOD = atoi( argv[2] );
    if (argc > 3)
        expo = atoi( argv[3] );
    if (argc > 4)
        Early_Term = atoi( argv[4] );
        // Field


    typedef SparseBlackBoxDom< GFqDomain > SPBB ;
    unsigned long rk1, rk2;
    GFqDomain::Rep valence;

    typedef Poly1Dom< GFqDomain, Dense > Polys;

    Indeter X("X");
    GFqDomain F(MOD,expo);
    Polys PolD( F, X);
//     Polys PolD( GFqDomain(MOD,expo) , X);
    Polys::element Q,G,P1,P2;
    Degree dq;

    init.stop();
    cerr << "Init : " << init.usertime() << " (" << init.realtime() << ")" << endl;

    SPBB MF( F );
//     SPBB MF( PolD.getdomain() );
    {

            // Default seed is timeofday.
        Random generator;
        
        MF.read(Matrix_File_Name);
        MF.precondition(generator);

        typedef BB_Symmetrize_Container< SPBB > SzCBB;


    {
        SzCBB TF( &MF, generator); 
        // MasseyDom< SzCBB >  WD(Commentator(PRINT_NOTHING,PRINT_NOTHING,2,4,EXPO_ESTIMATE), &TF );
        MasseyDom< SzCBB >  WD(Commentator(PRINT_EVERYTHING,PRINT_EVERYTHING,2,4,EXPO_ESTIMATE), &TF, Early_Term );
        WD.pseudo_minpoly(P1, rk1);
    }
    {
            // A new initial vector is used via the constructor
        SzCBB TF( &MF, generator);
        MasseyDom< SzCBB >  WD(Commentator(PRINT_EVERYTHING,PRINT_EVERYTHING,2,4,EXPO_ESTIMATE), &TF, Early_Term );
        WD.pseudo_minpoly(P2, rk2);
    }

    }
    
    
    PolD.gcd(G,P1,P2);
    PolD.div(Q,P1,G);
    PolD.mulin(Q,P2);
    PolD.degree(dq,Q);
    
//     PolD.write(cerr << "P1:=", P1) << ";" << endl;
//     PolD.write(cerr << "P2:=", P2) << ";" << endl;
//     PolD.write(cerr << "Gcd(P1,P2) mod " << MOD << " =", G) << ";" << endl;
//     PolD.write(cerr << "Lcm(P1,P2) mod " << MOD << " =", Q) << ";" << endl;

    cerr << "Degree ppcm : " << dq << endl;
    if ((dq > MF.n_col()) || (dq > MF.n_row())) cerr << "****** Failure of early terminaison parameter : " << Early_Term << endl;
    
//     } catch (GivError e) {
//         cerr << e << endl;
//     }

    Givaro::End();
    return 0;
};

