// ======================================================================= //
// Wiedemann algorithm using Massey
// With diagonal Scaling and Transpose Computation
// Time-stamp: <04 Aug 01 13:41:47 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#include <stdio.h>
#include <stdlib.h>

// ---------------------------------------------
#include <givtimer.h>
#include <givinit.h>
#include "givgfq.h"
#include <givpoly1.h>
typedef GFqDom<long>               GFqDomain;
typedef GFqDomain::Residu_t        Residu;
   

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
        // Compute the rank of B^T B
        // with B == A, and then B == D1 A^T D2^2 A D1
    typedef BB_Symmetrize_Container< SPBB > SzCBB;

    unsigned long rk1, rk2;
    GFqDomain::element valence;

    typedef Poly1Dom< GFqDomain, Dense > Polys;

    Indeter X("X");
    GFqDomain F(MOD,expo);
    Polys PolD( F, X);
    Polys::element Q,G,P1,P2;
    Degree dq;

    init.stop();
    cerr << "Init : " << init.usertime() << " (" << init.realtime() << ")" << endl;

    SPBB MF( F );
    {

            // Default seed is timeofday.
        Random generator;
        
        MF.read(Matrix_File_Name);

        cerr << "------------------------------------" << endl;
        cerr << "No preconditionning" << endl<< endl;
        

        GFqDomain::element tr;
        F.write( cerr << "trace_ata : ", MF.trace_ata(tr) ) << endl;




        {
            SzCBB TF( &MF, generator); 
            MasseyDom< SzCBB >  WD(Commentator(PRINT_NOTHING,PRINT_NOTHING,2,4,EXPO_ESTIMATE), &TF );

#ifdef __GIVARO_COUNT__
            F.info(); F.clear();
#endif // __GIVARO_COUNT__

            init.start();
            WD.pseudo_minpoly(P1, rk1);
            init.stop();

#ifdef __GIVARO_COUNT__
            F.info(); 
#endif // __GIVARO_COUNT__
            GFqDomain::element tr1; F.neg(tr1, P1[P1.size()-2]);
            F.write( cerr << "check trace : ", tr1 ) << endl;
            cerr << "rank : " << rk1 << endl;
            cerr << "time : " << init << endl;
        }

        cerr << "------------------------------------" << endl << endl;
        cerr << "------------------------------------" << endl;
        cerr << "Double diagonal preconditionning" << endl;
        cerr << "Two wiedemann calls + ppcm" << endl << endl;

        init.start();
        MF.precondition(generator);

        F.write( cerr << "trace_ata : ", MF.trace_ata(tr) ) << endl;


        {
            SzCBB TF( &MF, generator); 
            MasseyDom< SzCBB >  WD(Commentator(PRINT_NOTHING,PRINT_NOTHING,2,4,EXPO_ESTIMATE), &TF );
#ifdef __GIVARO_COUNT__
            F.info(); F.clear();
#endif // __GIVARO_COUNT__

            WD.pseudo_minpoly(P1, rk1);


#ifdef __GIVARO_COUNT__
            F.info(); 
#endif // __GIVARO_COUNT__
            GFqDomain::element tr1; F.neg(tr1, P1[P1.size()-2]);
            F.write( cerr << "check trace 1 : ", tr1 ) << endl;
            cerr << "rank 1 : " << rk1 << endl;
        }
        {
                // A new initial vector is used via the constructor
            SzCBB TF( &MF, generator);
            MasseyDom< SzCBB >  WD(Commentator(PRINT_EVERYTHING,PRINT_EVERYTHING,2,4,EXPO_ESTIMATE), &TF, Early_Term );
            WD.pseudo_minpoly(P2, rk2);
            GFqDomain::element tr2; F.neg(tr2, P2[P2.size()-2]);
            F.write( cerr << "check trace 2 : ", tr2 ) << endl;
        }

    }
    
    
    PolD.gcd(G,P1,P2);
    PolD.div(Q,P1,G);
    PolD.mulin(Q,P2);
    PolD.degree(dq,Q);

    GFqDomain::element tr3; F.neg(tr3, Q[Q.size()-2]);
    F.write( cerr << "check trace cm : ", tr3 ) << endl;
/*    
    PolD.write(cerr << "P1:=", P1) << ";" << endl;
    PolD.write(cerr << "P2:=", P2) << ";" << endl;
    PolD.write(cerr << "Gcd(P1,P2) mod " << MOD << " =", G) << ";" << endl;
    PolD.write(cerr << "Lcm(P1,P2) mod " << MOD << " =", Q) << ";" << endl;
*/
    cerr << "Degree lcm : " << dq << endl;
    if ((dq > MF.n_col()) || (dq > MF.n_row())) cerr << "****** Failure of early terminaison parameter : " << Early_Term << endl;
    init.stop();
    cerr << "global time: " << init << endl;
    cerr << "------------------------------------" << endl << endl;
   
 
    Givaro::End();
    return 0;
};

