//      Gen_SuperLU package (version 1.0)
// genlinsol.C was made generic by Ahmet Duran and David Saunders 
// so that it works for LinBox::Modular Field , LinBox::UnparametricField 
// and LinBox::GivaroZpz Field.
// Computer & Information Science Department
// University of Delaware
// August 2002

/* Evolved from 
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 */
#include <iostream.h>
#include <iomanip.h>
#include <stdlib.h>
#include <time.h>
#include "linbox/util/timer.h"

#include "util.h"

// #include "linbox/field/unparametric.h"
// #include "linbox/integer.h"
// #include "linbox/field/modular.h"
// #include "linbox/field/givaro-zpz.h"

#include "linbox-config.h"
#include "linbox/field/ntl-zz_p.h"
#include "NTL/lzz_p.h"

#include "sp_defs.h"
// #include "gssv.h" // include if both factorization and solution are needed.
#include "gssv_rank.h" // include if rank and LU factorization are needed only
#include "fmemory.h"
#include "futil.h"
#include "readtriple1.h"
#include "determinant.h"
// #include "rank.h"

// "linbox/field/modular.h" should be included
// typedef LinBox::Modular<LinBox::uint32> Field;
// Field F(101);
// Field F(32717);

// "linbox/field/givaro-zpz.h" should be included
// typedef LinBox::GivaroZpz<Std32> Field;
// Field F(101);

// "linbox/field/ntl-zz_p.h" should be included
using namespace NTL;
typedef LinBox::UnparametricField<NTL::zz_p> Field;
Field F;
#define NTLZZP // provides rep

// typedef LinBox::UnparametricField<double> Field;
// Field F;
// #define UNPARAMETRICFIELDFD // float or double

using namespace LinBox;
main(int argc, char *argv[])
{
  long q = 101;  // For LinBox::UnparametricField<NTL::zz_p>
  zz_p::init(q); // For LinBox::UnparametricField<NTL::zz_p>
  
  SuperMatrix<Field> A, L, U, B;
  Field::Element *a, *rhs;
  int      *asub, *xa;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      nrhs, info, i, m, n, nnz, permc_spec;
  int rank; // A. Duran    
  int det;  // A. Duran

  UserTimer T;
  time_t t;

  T.start();

  /* Initialize matrix A. */
  readtriple(&m, &n, &nnz, &a, &asub, &xa, F);
        
  /* Create matrix A in the format expected by SuperLU. */
  FCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, NC, D_, GE, F);
  
  /* Create right-hand side matrix B. */
  // nrhs = 1;
  // if ( !(rhs = FieldMalloc(m * nrhs, F)) ) ABORT("Malloc fails for rhs[].");
  // F.init(rhs[0], 1);
  // for (i = 1; i < m; ++i) F.init(rhs[i], 1);
  // FCreate_Dense_Matrix(&B, m, nrhs, rhs, m, DN, D_, GE, F);
  
  if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
  if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
  
  /*
   * Get column permutation vector perm_c[], according to permc_spec:
   *   permc_spec = 0: use the natural ordering 
   *   permc_spec = 1: use minimum degree ordering on structure of A'*A
   *   permc_spec = 2: use minimum degree ordering on structure of A'+A
   *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
   */    	
  permc_spec = 1;
  get_perm_c(permc_spec, &A, perm_c);
  
  // gssv(&A, perm_c, perm_r, &L, &U, &B, &info, &rank, F);
  // Use this if only rank is needed, becaude we don't need B.
  gssv_rank(&A, perm_c, perm_r, &L, &U, &info, &rank, F); 
  T.stop();
  
  cout << "rank of A :"<< rank <<"\n";     
  cout << T << endl;
  /* dPrint_CompCol_Matrix("A", &A); */
  // Print_Dense_Matrix("B", &B, F); /* A.Duran */
  // print_int_vec("\nperm_r", m, perm_r);
  // print_int_vec("\nperm_c", m, perm_c);
  // Print_CompCol_Matrix("U", &U, F);
  // Print_SuperNode_Matrix("L", &L, F);
  // print_int_vec("\nperm_r", m, perm_r);
  // print_int_vec("\nperm_c", m, perm_c);
  
  // determinant("A", &L, m, perm_r, perm_c, rank, F); // A. Duran
  
  // rank("A", &L, &U, F); // A. Duran
  
  /* De-allocate storage */
  SUPERLU_FREE (rhs);
  SUPERLU_FREE (perm_r);
  SUPERLU_FREE (perm_c);
  Destroy_CompCol_Matrix(&A);
  Destroy_SuperMatrix_Store(&B);
  Destroy_SuperNode_Matrix(&L);
  Destroy_CompCol_Matrix(&U);
  
}
