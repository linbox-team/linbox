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
#include <typeinfo>
#include "linbox/gen_superlu/util.h"
#include "linbox/util/commentator.h"
#include "test-common.h"
#include "test-common.C"

#include "linbox/field/ntl-zz_p.h"

#include "linbox-config.h"
#include "linbox/gen_superlu/sp_defs.h"
// #include "linbox/gen_superlu/gssv.h" // include if both factorization and solution are needed.
#include "linbox/gen_superlu/gssv_rank.h" // include if rank and LU factorization are needed only
#include "linbox/gen_superlu/fmemory.h"
#include "linbox/gen_superlu/futil.h"
#include "linbox/gen_superlu/readtriple.h"
// #include "linbox/gen_superlu/determinant.h"
// #include "linbox/gen_superlu/rank.h"

using namespace LinBox;
using namespace NTL;
template<class Field>
bool  test_gssv_rank(Field& F);

int main(int argc, char *argv[])
{
	typedef LinBox::UnparametricField<NTL::zz_p> Field;
 	Field F;
		
    	int modulus[]  = {65521,1048573,16777259,268435459,1073741789};
	int trial = 0;

	long q = modulus[trial];  // For LinBox::UnparametricField<NTL::zz_p>
	zz_p::init(q); // For LinBox::UnparametricField<NTL::zz_p>

	bool pass = true;

	static Argument args[] = {};

	parseArguments (argc, argv, args);
	  
	pass = test_gssv_rank<Field>(F);

	return pass;
}


template <class Field>
bool test_gssv_rank(Field& F) {	
	
    UserTimer T;
    time_t t;
  
    int      *asub, *xa;
    int      *perm_r; /* row permutations from partial pivoting */
    int      *perm_c; /* column permutation vector */
    int      info, i, m, n, nnz, permc_spec;
    int rank; // A. Duran    
    int det;  // A. Duran
  
    FILE *fp;
  
    bool pass = true;

    ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

    report << endl << "gen-superlu test suite" << endl;

    for (permc_spec = 0; permc_spec < 4; ++permc_spec)
	{
   	  // min deg ordering on A'+A doesn't work
	  if (permc_spec == 2) continue; 

          // construct field
	  SuperMatrix<Field> A, L, U;
	  typename Field::Element *a;
	  
	  /* Initialize matrix A. */
	  fp = fopen("data/gssv_rank_data", "r");
	  readtriple(&m, &n, &nnz, &a, &asub, &xa, fp, F);
	  fclose(fp);
	  
	  T.start();
	  
	  /* Create matrix A in the format expected by SuperLU. */
	  FCreate_CompCol_Matrix(&A, m, n, nnz, a, asub, xa, NC, D_, GE, F);
	  
	  if ( !(perm_r = intMalloc(m)) ) ABORT("Malloc fails for perm_r[].");
	  if ( !(perm_c = intMalloc(n)) ) ABORT("Malloc fails for perm_c[].");
	  
	  /*
	   * Get column permutation vector perm_c[], according to permc_spec:
	   *   permc_spec = 0: use the natural ordering 
	   *   permc_spec = 1: use minimum degree ordering on structure of A'*A
	   *   permc_spec = 2: use minimum degree ordering on structure of A'+A
	   *   permc_spec = 3: approximate minimum degree for unsymmetric matrices
	   */    	
	  // permc_spec = 1;
	  get_perm_c(permc_spec, &A, perm_c);
	  
	  // gssv(&A, perm_c, perm_r, &L, &U, &B, &info, &rank, F);
	  // Use this if only rank is needed, becaude we don't need B.
	  gssv_rank(&A, perm_c, perm_r, &L, &U, &info, &rank, F); 
	  pass &= (rank == n); // assuming full rank example.

	  T.stop();
	  
	  F.write(report);

	  report << T << " sec\n\n";
	  /* dPrint_CompCol_Matrix("A", &A); */
	  // Print_Dense_Matrix("B", &B, F); /* A.Duran */
	  // print_int_vec("\nperm_r", m, perm_r);
	  // print_int_vec("\nperm_c", m, perm_c);
	  // Print_CompCol_Matrix("U", &U, F);
	  // Print_SuperNode_Matrix("L", &L, F);
	  // print_int_vec("\nperm_r", m, perm_r);
	  // print_int_vec("\nperm_c", m, perm_c);
	  
	  /* De-allocate storage */
	  //  SUPERLU_FREE (rhs);
	  SUPERLU_FREE (perm_r);
	  SUPERLU_FREE (perm_c);
	  Destroy_CompCol_Matrix(&A);
	  // Destroy_SuperMatrix_Store(&B);
	  Destroy_SuperNode_Matrix(&L);
	  Destroy_CompCol_Matrix(&U);
	}

    return pass ? 0 : -1;
}
