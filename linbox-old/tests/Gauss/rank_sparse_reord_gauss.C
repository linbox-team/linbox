// ======================================================================= //
// Givaro / Athapascan-1
// Algorithme de Coppersmith : n=m=1 --> Wiedemann ameliore, Creux
// With diagonal Scaling and Transpose Computation
// Multiply computations are stopped when the polynomials remains the same
// for more than DIFF_BOUND
// Time-stamp: <08 Mar 00 16:21:40 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
// ---------------------------------------------
#include <givinit.h>
#include "givgfq.h"

#include "lin_spv_bb.h"                // BB Wrapper for sparse vectors
#include "lin_density.h"
#include "lin_dom_gauss.C"   
// ---------------------------------------------
// MAIN

int main(int argc, char* argv[]) {
    Givaro::Init(&argc, &argv);
    Timer init; init.clear();init.start();
    
    typedef GFqDom<long>               GFqDomain;
    typedef GFqDomain::Residu_t        Residu;
   
        // Must have a matrix name    
    char* Matrix_File_Name = argv[1];
    Residu MOD=65521,expo=1,cardinal;
    short transp = -1;
    if (argc > 2)
        MOD = atoi( argv[2] );
    if (argc > 3)
        expo = atoi( argv[3] );
    if (argc > 4)
        transp = atoi( argv[4] );
        // Première étape : le corps de base
    GFqDomain F(MOD,expo);
    cardinal = F.size();
    GaussDom< GFqDomain > GF( F, Commentator(PRINT_EVERYTHING,PRINT_EVERYTHING) );    

    typedef SparseBlackBoxDom< GFqDomain > SPBB ;
    SPBB TheMat(F);
    unsigned long ni, nj; TheMat.read(ni,nj,Matrix_File_Name);
    if (transp < 0)
    	if (ni > nj)
            TheMat.read_transpose(Matrix_File_Name);
    	else
            TheMat.read(Matrix_File_Name);
    else 
        if (transp == 1)
            TheMat.read_transpose(Matrix_File_Name);
    	else
            TheMat.read(Matrix_File_Name);

    unsigned long rank;
    init.stop();
    cerr << "Init, " << TheMat.n_elem() << " elts : " << init.usertime() 
         << " (" << init.realtime() << ")" << endl;
    
    {
        
    GF.gauss_rankin(rank,TheMat,Density());
  
   }
    
    Givaro::End();
    return 0;
};
