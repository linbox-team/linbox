//      Gen_SuperLU package (version 1.0)
// gensuperlu.C was made generic by Ahmet Duran and David Saunders 
// so that it works for LinBox::Modular Field , LinBox::UnparametricField 
// and LinBox::GivaroZpz Field.
// Computer & Information Science Department
// University of Delaware
// August 2002

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 */
#include <iostream.h>

#include "../src/util.h"

#include "linbox/field/unparametric.h"
#include "linbox/integer.h"
#include "linbox/field/modular.h"
#include "linbox/field/givaro-zpz.h"

#include "sp_defs.h"
#include "gssv.h"
#include "fmemory.h"
#include "futil.h"
#include "readtriple.h"
#include "determinant.h"


typedef LinBox::Modular<LinBox::uint32> Field;
Field F(101);

// typedef LinBox::GivaroZpz<Std32> Field;
// Field F(101);

// typedef LinBox::UnparametricField<double> Field;
// Field F;



main(int argc, char *argv[])
{
/*
 * Purpose
 * =======
 * 
 * This is the small 5x5 example used in the Sections 1 and 2 of the 
 * User's Guide to illustrate how to call a SuperLU routine, and the
 * matrix data structures used by SuperLU.
 *
 */
    SuperMatrix<Field> A, L, U, B;
    Field::Element *a, *rhs;
    Field::Element  s, u, p, e, r, l;
    int      *asub, *xa;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      nrhs, info, i, m, n, nnz, permc_spec;
    
    /* Initialize matrix A. */
    
    m = n = 5;
    nnz = 12;
    if ( !(a = FieldMalloc(nnz, F)) ) ABORT("Malloc fails for a[].");
    if ( !(asub = intMalloc(nnz)) ) ABORT("Malloc fails for asub[].");
    if ( !(xa = intMalloc(n+1)) ) ABORT("Malloc fails for xa[].");
    F.init(s, 19); F.init(u, 21); F.init(p, 16); F.init(e, 5); F.init(r, 18); 
    F.init(l, 12);
    F.assign(a[0], s); F.assign(a[1], l); F.assign(a[2], l); F.assign(a[3], u); 
    F.assign(a[4], l); F.assign(a[5], l); F.assign(a[6], u); F.assign(a[7], p); 
    F.assign(a[8], u); F.assign(a[9], e); F.assign(a[10], u); F.assign(a[11], r);
    asub[0] = 0; asub[1] = 1; asub[2] = 4; asub[3] = 1;
    asub[4] = 2; asub[5] = 4; asub[6] = 0; asub[7] = 2;
    asub[8] = 0; asub[9] = 3; asub[10]= 3; asub[11]= 4;
    xa[0] = 0; xa[1] = 3; xa[2] = 6; xa[3] = 8; xa[4] = 10; xa[5] = 12;
        
    /* Create matrix A in the format expected by SuperLU. */
    FCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, NC, D_, GE, F);
    
    /* Create right-hand side matrix B. */
    nrhs = 1;
    if ( !(rhs = FieldMalloc(m * nrhs, F)) ) ABORT("Malloc fails for rhs[].");
    F.init(rhs[0], 1);
    for (i = 1; i < m; ++i) F.init(rhs[i], 1);
    FCreate_Dense_Matrix(&B, m, nrhs, rhs, m, DN, D_, GE, F);

    if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
    if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");

    /*
     * Get column permutation vector perm_c[], according to permc_spec:
     *   permc_spec = 0: use the natural ordering 
     *   permc_spec = 1: use minimum degree ordering on structure of A'*A
     *   permc_spec = 2: use minimum degree ordering on structure of A'+A
     *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
     */    	
    permc_spec = 2;
    get_perm_c(permc_spec, &A, perm_c);

    gssv(&A, perm_c, perm_r, &L, &U, &B, &info, F);
    
    /* dPrint_CompCol_Matrix("A", &A); */
    Print_Dense_Matrix("B", &B, F); /* A.Duran */
    Print_CompCol_Matrix("U", &U, F);
    Print_SuperNode_Matrix("L", &L, F);
    print_int_vec("\nperm_r", m, perm_r);
    print_int_vec("\nperm_c", m, perm_c);
    determinant("A", &L, m, perm_r, perm_c, F); // A. Duran

    /* De-allocate storage */
    SUPERLU_FREE (rhs);
    SUPERLU_FREE (perm_r);
    SUPERLU_FREE (perm_c);
    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
}
