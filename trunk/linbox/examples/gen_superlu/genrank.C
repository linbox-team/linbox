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

// #define MODULARFIELD
#define UNPARAMETRICZZP
// #define GIVAROFIELD
// #define UNPARAMETRICFIELDFLOAT
// #define UNPARAMETRICFIELDDOUBLE

#ifdef UNPARAMETRICFIELDFLOAT
  #include "linbox/field/unparametric.h"
#endif
#ifdef UNPARAMETRICFIELDDOUBLE
 #include "linbox/field/unparametric.h"
#endif

// #include "linbox/integer.h"

#ifdef MODULARFIELD
  #include "linbox/field/modular.h"
#endif

#ifdef GIVAROFIELD
  #include "linbox/field/givaro-zpz.h"
#endif

#include "linbox-config.h"

#ifdef UNPARAMETRICZZP
  #include "linbox/field/ntl-zz_p.h"
  #include "NTL/lzz_p.h"
#endif

#include "linbox/gen_superlu/sp_defs.h"
// #include "linbox/gen_superlu/gssv.h" // include if both factorization and solution are needed.
#include "linbox/gen_superlu/gssv_rank.h" // include if rank and LU factorization are needed only
#include "linbox/gen_superlu/fmemory.h"
#include "linbox/gen_superlu/futil.h"
#include "linbox/gen_superlu/readtriple.h"
// #include "linbox/gen_superlu/determinant.h"
// #include "linbox/gen_superlu/rank.h"

using namespace LinBox;
main(int argc, char *argv[])
{
  UserTimer T;
  time_t t;
  int modulus[]  = {65521,1048573,16777259,268435459,1073741789};
  // int modulus[]  = {1048573, 1073741789}; 
  
  int      *asub, *xa;
  int      *perm_r; /* row permutations from partial pivoting */
  int      *perm_c; /* column permutation vector */
  int      info, i, m, n, nnz, permc_spec;
  int rank; // A. Duran    
  int det;  // A. Duran
  
  FILE *fp;
  

  int trial;
  
  for (trial = 0; trial < 2; ++trial)
    {
      for (permc_spec = 0; permc_spec < 4; ++permc_spec)
	{
#ifdef MODULARFIELD  
	  // "linbox/field/modular.h" should be included
	  typedef LinBox::Modular<LinBox::uint32> Field;
	  Field F(modulus[trial]);
	  // Field F(32717);
#endif
	  
#ifdef GIVAROFIELD
	  // "linbox/field/givaro-zpz.h" should be included
	  typedef LinBox::GivaroZpz<Std32> Field;
	  Field F(modulus[trial]);	  
#endif
	  
#ifdef UNPARAMETRICZZP
	  // "linbox/field/ntl-zz_p.h" should be included
	  using namespace NTL;
	  typedef LinBox::UnparametricField<NTL::zz_p> Field;
	  Field F;

	  long q = modulus[trial];  // For LinBox::UnparametricField<NTL::zz_p>
	  zz_p::init(q); // For LinBox::UnparametricField<NTL::zz_p>
#define NTLZZP // provides rep
#endif
	  
#ifdef UNPARAMETRICFIELDFLOAT
	  typedef LinBox::UnparametricField<float> Field;
	  Field F;
#define UNPARAMETRICFIELDFD // float or double
#endif
	  
#ifdef UNPARAMETRICFIELDDOUBLE
	  typedef LinBox::UnparametricField<double> Field;
	  Field F;
#define UNPARAMETRICFIELDFD // float or double
#endif  
	  
	  SuperMatrix<Field> A, L, U;
	  Field::Element *a;
	  
	  /* Initialize matrix A. */
	  fp = fopen(argv[1], "r");
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

	  T.stop();
	  
	  cout << typeid(F).name()<<"\n";
	  cout << "File name : " << argv[1] <<"\n";
	  cout << "Modulus : " << modulus[trial] << "\nRank of A : "<< rank <<"\n";     
	  cout << T << " sec\n\n";
	  if (permc_spec == 3) 
	    cout << "**********************\n";
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
	  //  SUPERLU_FREE (rhs);
	  SUPERLU_FREE (perm_r);
	  SUPERLU_FREE (perm_c);
	  Destroy_CompCol_Matrix(&A);
	  // Destroy_SuperMatrix_Store(&B);
	  Destroy_SuperNode_Matrix(&L);
	  Destroy_CompCol_Matrix(&U);
	}
      if (trial == 2) 
	cout << "$$$$$$$$$$$$$$$$$$$$$$$$$\n";
    }
}
