// ======================================================================= //
// Givaro / Athapascan-1
// Algorithme de Coppersmith : n=m=1 --> Wiedemann ameliore, Creux
// With diagonal Scaling and Transpose Computation
// Multiply computations are stopped when the polynomials remains the same
// for more than DIFF_BOUND
// Time-stamp: <08 Mar 00 16:21:40 Jean-Guillaume.Dumas@imag.fr> 
// ======================================================================= //
#define _REENTRANT

#include <stdio.h>
#include <stdlib.h>
#include <iostream.h>
// ---------------------------------------------

#include "LinBox/integer.h" 
using namespace LinBox;
#include "LinBox/field_archetype.h"

#include "LinBox/gmp-rational-field.C"
#include "LinBox/gmp-rational-random.h"



#include "LinBox/lin_spv_bb.h"                // BB Wrapper for sparse vectors
#include "LinBox/lin_density.h"
#include "LinBox/lin_dom_gauss.C"   
// ---------------------------------------------
// MAIN

int main(int argc, char* argv[]) {
    Timer init; init.clear();init.start();
    
   
        // Must have a matrix name    
    char* Matrix_File_Name = argv[1];
    short transp = -1;
    if (argc > 4)
        transp = atoi( argv[4] );
        // Première étape : le corps de base


    typedef GMP_Rational_Field MyField;
    MyField F;

	
#ifdef __ARCH__
    Field_archetype Af( & F );

//  works also ...
//    Field_archetype Af( new Field_envelope<MyField>(F) );

    GaussDom< Field_archetype > elimination_routines_with_field_stored( Af, Commentator(PRINT_EVERYTHING,PRINT_EVERYTHING) );    
    SparseBlackBoxDom< Field_archetype > TheMat(  Af );
#else
    GaussDom< MyField > elimination_routines_with_field_stored( F, Commentator(PRINT_EVERYTHING,PRINT_EVERYTHING) );    
    SparseBlackBoxDom< MyField > TheMat(F);
#endif



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
        
    elimination_routines_with_field_stored.gauss_rankin(rank,TheMat,Density());
  
   }
    
    return 0;
};
