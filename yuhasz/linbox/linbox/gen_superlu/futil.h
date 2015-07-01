

/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
/*
  Copyright (c) 1994 by Xerox Corporation.  All rights reserved.
 
  THIS MATERIAL IS PROVIDED AS IS, WITH ABSOLUTELY NO WARRANTY
  EXPRESSED OR IMPLIED.  ANY USE IS AT YOUR OWN RISK.
 
  Permission is hereby granted to use or copy this program for any
  purpose, provided the above notices are retained on all copies.
  Permission to modify the code and to distribute modified code is
  granted, provided the above notices are retained, and a notice that
  the code was modified is included with the above copyright notice.
*/

#include <math.h>
#include "linbox/gen_superlu/sp_defs.h"
#include "linbox/gen_superlu/util.h"


// A.Duran 8/14/2002
template <class Field>
void
FCreate_CompCol_Matrix(SuperMatrix<Field> *A, int m, int n, int nnz, 
		       typename Field::Element *nzval, int *rowind, int *colptr,
		       Stype_t stype, Dtype_t dtype, Mtype_t mtype, Field& F)
{
    NCformat<Field> *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    // A->Store = (void *) SUPERLU_MALLOC( sizeof(NCformat) );
    A->Store = (typename Field::Element *) SUPERLU_MALLOC( sizeof(NCformat<Field>) ); // A.Duran
    if ( !(A->Store) ) ABORT("SUPERLU_MALLOC fails for A->Store");
    Astore = (NCformat<Field> *)A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->rowind = rowind;
    Astore->colptr = colptr;
}


// A.Duran 8/14/2002
template <class Field>
void
FCreate_CompRow_Matrix(SuperMatrix<Field> *A, int m, int n, int nnz, 
		       typename Field::Element *nzval, int *colind, int *rowptr,
		       Stype_t stype, Dtype_t dtype, Mtype_t mtype, Field& F)
{
    NRformat<Field> *Astore;

    A->Stype = stype;
    A->Dtype = dtype;
    A->Mtype = mtype;
    A->nrow = m;
    A->ncol = n;
    // A->Store = (void *) SUPERLU_MALLOC( sizeof(NRformat) );
    A->Store = (typename Field::Element *) SUPERLU_MALLOC( sizeof(NRformat<Field>) ); // A.Duran
    if ( !(A->Store) ) ABORT("SUPERLU_MALLOC fails for A->Store");
    Astore = (NRformat<Field> *)A->Store;
    Astore->nnz = nnz;
    Astore->nzval = nzval;
    Astore->colind = colind;
    Astore->rowptr = rowptr;
}


/* Copy matrix A into matrix B. */
template <class Field>
void
FCopy_CompCol_Matrix(SuperMatrix<Field> *A, SuperMatrix<Field> *B, Field& F)
{
    NCformat<Field> *Astore, *Bstore;
    int      ncol, nnz, i;

    B->Stype = A->Stype;
    B->Dtype = A->Dtype;
    B->Mtype = A->Mtype;
    B->nrow  = A->nrow;;
    B->ncol  = ncol = A->ncol;
    Astore   = (NCformat<Field> *) A->Store;
    Bstore   = (NCformat<Field> *) B->Store;
    Bstore->nnz = nnz = Astore->nnz;
    for (i = 0; i < nnz; ++i)
	((typename Field::Element *)Bstore->nzval)[i] = ((typename Field::Element *)Astore->nzval)[i];
    for (i = 0; i < nnz; ++i) Bstore->rowind[i] = Astore->rowind[i];
    for (i = 0; i <= ncol; ++i) Bstore->colptr[i] = Astore->colptr[i];
}


// A.Duran 8/14/2002
template <class Field>
void
FCreate_Dense_Matrix(SuperMatrix<Field> *X, int m, int n, typename Field::Element *x, int ldx,
		    Stype_t stype, Dtype_t dtype, Mtype_t mtype, Field& F)
{
    DNformat<Field>    *Xstore;
    
    X->Stype = stype;
    X->Dtype = dtype;
    X->Mtype = mtype;
    X->nrow = m;
    X->ncol = n;
    // X->Store = (void *) SUPERLU_MALLOC( sizeof(DNformat) );
    X->Store = (typename Field::Element *) SUPERLU_MALLOC( sizeof(DNformat<Field>) );
    if ( !(X->Store) ) ABORT("SUPERLU_MALLOC fails for X->Store");
    Xstore = (DNformat<Field> *) X->Store;
    Xstore->lda = ldx;
    Xstore->nzval = (typename Field::Element *) x;
}

template <class Field>
void
FCopy_Dense_Matrix(int M, int N, typename Field::Element *X, int ldx,
			typename Field::Element *Y, int ldy, Field& F)
{
/*
 *
 *  Purpose
 *  =======
 *
 *  Copies a two-dimensional matrix X to another matrix Y.
 */
    int    i, j;
    
    for (j = 0; j < N; ++j)
        for (i = 0; i < M; ++i)
            Y[i + j*ldy] = X[i + j*ldx];
}


template <class Field>
void
FCreate_SuperNode_Matrix(SuperMatrix<Field> *L, int m, int n, int nnz, 
			typename Field::Element *nzval, int *nzval_colptr, int *rowind,
			int *rowind_colptr, int *col_to_sup, int *sup_to_col,
			Stype_t stype, Dtype_t dtype, Mtype_t mtype, Field& F)
{
    SCformat<Field> *Lstore;

    L->Stype = stype;
    L->Dtype = dtype;
    L->Mtype = mtype;
    L->nrow = m;
    L->ncol = n;
    // L->Store = (void *) SUPERLU_MALLOC( sizeof(SCformat) );
    L->Store = (typename Field::Element *) SUPERLU_MALLOC( sizeof(SCformat<Field>) );
    if ( !(L->Store) ) ABORT("SUPERLU_MALLOC fails for L->Store");
    Lstore = (SCformat<Field> *)L->Store;
    Lstore->nnz = nnz;
    Lstore->nsuper = col_to_sup[n];
    Lstore->nzval = nzval;
    Lstore->nzval_colptr = nzval_colptr;
    Lstore->rowind = rowind;
    Lstore->rowind_colptr = rowind_colptr;
    Lstore->col_to_sup = col_to_sup;
    Lstore->sup_to_col = sup_to_col;

}

/*
 * Convert a row compressed storage into a column compressed storage.
 */
template <class Field>
void
CompRow_to_CompCol(int m, int n, int nnz, 
		    typename Field::Element *a, int *colind, int *rowptr,
		    typename Field::Element **at, int **rowind, int **colptr, Field& F)
{
    register int i, j, col, relpos;
    int *marker;

    /* Allocate storage for another copy of the matrix. */
    *at = (typename Field::Element *) FieldMalloc(nnz, F);
    *rowind = (int *) intMalloc(nnz);
    *colptr = (int *) intMalloc(n+1);
    marker = (int *) intCalloc(n);
    
    /* Get counts of each column of A, and set up column pointers */
    for (i = 0; i < m; ++i)
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) ++marker[colind[j]];
    (*colptr)[0] = 0;
    for (j = 0; j < n; ++j) {
	(*colptr)[j+1] = (*colptr)[j] + marker[j];
	marker[j] = (*colptr)[j];
    }

    /* Transfer the matrix into the compressed column storage. */
    for (i = 0; i < m; ++i) {
	for (j = rowptr[i]; j < rowptr[i+1]; ++j) {
	    col = colind[j];
	    relpos = marker[col];
	    (*rowind)[relpos] = i;
	    (*at)[relpos] = a[j];
	    ++marker[col];
	}
    }

    SUPERLU_FREE(marker);
}

template <class Field>
void
Print_CompCol_Matrix(char *what, SuperMatrix<Field> *A, Field& F)
{
    NCformat<Field>     *Astore;
    register int i,n;
    typename Field::Element *dp;
    
    // printf("\nCompCol matrix %s:\n", what);
    cout << "\nCompCol matrix "<<what<<endl;
    // printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    cout<<"Stype "<<A->Stype <<" Dtype "<< A->Dtype<<" Mtype "<<A->Mtype<<endl;
    n = A->ncol;
    Astore = (NCformat<Field> *) A->Store;
    dp = (typename Field::Element *) Astore->nzval;
    // printf("nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
    cout<<"nrow "<<A->nrow<<" ncol "<<A->ncol<<" nnz "<<Astore->nnz<<endl;
    // printf("nzval: ");
    cout<<"nzval: ";
    for (i = 0; i < Astore->colptr[n]; ++i) // printf("%f  ", dp[i]);
      cout <<dp[i]<<"  ";
    // printf("\nrowind: ");
    cout<<"\nrowind: ";
    for (i = 0; i < Astore->colptr[n]; ++i) // printf("%d  ", Astore->rowind[i]);
      cout<<Astore->rowind[i]<<"  ";
    // printf("\ncolptr: ");
    cout<<"\ncolptr: ";
    for (i = 0; i <= n; ++i) // printf("%d  ", Astore->colptr[i]);
      cout<<Astore->colptr[i]<<" ";
    // printf("\n");
    cout<<endl;
    fflush(stdout);
}


// A. Duran added Print_CompCol_NCP_Matrix function
// Actually, colend should be printed also.
template <class Field>
void
Print_CompCol_NCP_Matrix(char *what, SuperMatrix<Field> *A, Field& F)
{
    NCPformat<Field> *Astore;
    register int i,n;
    typename Field::Element *dp;
     fflush(stdout);
    // printf("\nCompCol matrix %s:\n", what);
    cout << "\nCompCol matrix "<<what<<endl;
    // printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    cout<<"Stype "<<A->Stype <<" Dtype "<< A->Dtype<<" Mtype "<<A->Mtype<<endl;
    n = A->ncol;
    Astore = (NCPformat<Field> *) A->Store;
    dp = (typename Field::Element *) Astore->nzval;
    // printf("nrow %d, ncol %d, nnz %d\n", A->nrow,A->ncol,Astore->nnz);
    cout<<"nrow "<<A->nrow<<" ncol "<<A->ncol<<" nnz "<<Astore->nnz<<endl;
    // printf("nzval: ");
    cout<<"nzval: ";
    for (i = 0; i < Astore->nnz; ++i) // printf("%f  ", dp[i]);
      cout <<dp[i]<<"  ";
    // printf("\nrowind: ");
    cout<<"\nrowind: ";
    for (i = 0; i < Astore->nnz; ++i) // printf("%d  ", Astore->rowind[i]);
      cout<<Astore->rowind[i]<<"  ";
    // printf("\ncolptr: ");
    cout<<"\n\ncolbeg: ";
    for (i = 0; i <= n; ++i) // printf("%d  ", Astore->colptr[i]);
      cout<<Astore->colbeg[i]<<" ";
    // printf("\n");
    cout<<"\ncolend: ";
    for (i = 0; i <= n; ++i) // printf("%d  ", Astore->colptr[i]);
      cout<<Astore->colend[i]<<" ";
    // printf("\n");
    cout<<endl;

    fflush(stdout);
}




template <class Field>
void
Print_SuperNode_Matrix(char *what, SuperMatrix<Field> *A, Field& F)
{
    SCformat<Field>     *Astore;
    register int i, j, k, c, d, n, nsup;
    typename Field::Element *dp;

    int *col_to_sup, *sup_to_col, *rowind, *rowind_colptr;
    
    // printf("\nSuperNode matrix %s:\n", what);
    cout<<"\nSuperNode matrix "<<what<<endl;
    // printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    cout<<"Stype "<<A->Stype<<" Dtype "<<A->Dtype<<" Mtype "<<A->Mtype<<endl;
    n = A->ncol;
    Astore = (SCformat<Field> *) A->Store;
    dp = (typename Field::Element *) Astore->nzval;
    col_to_sup = Astore->col_to_sup;
    sup_to_col = Astore->sup_to_col;
    rowind_colptr = Astore->rowind_colptr;
    rowind = Astore->rowind;
    // printf("nrow %d, ncol %d, nnz %d, nsuper %d\n", 
    //	   A->nrow,A->ncol,Astore->nnz,Astore->nsuper);
    cout<<"nrow "<<A->nrow<<" ncol "<<A->ncol<<" nnz "<<Astore->nnz<<" nsuper "<<Astore->nsuper<<endl;
    // printf("nzval:\n");
    cout<<"nzval:"<<endl;
    for (k = 0; k <= Astore->nsuper+1; ++k) {
      c = sup_to_col[k];
      nsup = sup_to_col[k+1] - c;
      // A. Duran added  condition (j < n) . Otherwise, there was overflow!
      for (j = c; (j < c + nsup) && (j < n); ++j) { 
	d = Astore->nzval_colptr[j];
	// A. Duran added  condition (rowind[i] > -1 )
	for (i = rowind_colptr[c]; (i < rowind_colptr[c+1]) && (rowind[i] > -1 ); ++i) {
	  // printf("%d\t%d\t%e\n", rowind[i], j, dp[d++]);
	  cout<<rowind[i]<<"\t"<< j<<"\t"<< dp[d++]<<endl;
	}
      }
    }

#if 0
    for (i = 0; i < Astore->nzval_colptr[n]; ++i) // printf("%f  ", dp[i]);
      cout<<dp[i]<<"  ";
#endif
    // printf("\nnzval_colptr: ");
    cout<<"\nnzval_colptr: ";
    for (i = 0; i <= n; ++i) // printf("%d  ", Astore->nzval_colptr[i]);
      cout<<Astore->nzval_colptr[i]<<"  ";
    // printf("\nrowind: ");
    cout<<"\nrowind: ";
    for (i = 0; i < Astore->rowind_colptr[n]; ++i) 
      // printf("%d  ", Astore->rowind[i]);
      cout<<Astore->rowind[i]<<"  ";
    // printf("\nrowind_colptr: ");
    cout<<"\nrowind_colptr: ";
    for (i = 0; i <= n; ++i) // printf("%d  ", Astore->rowind_colptr[i]);
      cout<<Astore->rowind_colptr[i]<<"  ";
    // printf("\ncol_to_sup: ");
    cout<<"\ncol_to_sup: ";
    for (i = 0; i < n; ++i) // printf("%d  ", col_to_sup[i]);
      cout<<col_to_sup[i]<<"  ";
    // printf("\nsup_to_col: ");
    cout<<"\nsup_to_col: ";
    for (i = 0; i <= Astore->nsuper+1; ++i) 
      // printf("%d  ", sup_to_col[i]);
      cout<<sup_to_col[i]<<"  ";
    // printf("\n");
    cout<<endl;
    fflush(stdout);
}

template <class Field>
void
Print_Dense_Matrix(char *what, SuperMatrix<Field> *A, Field& F)
{
    DNformat<Field> *Astore;
    register int i;
    typename Field::Element *dp;
    
    // printf("\nDense matrix %s:\n", what);
    cout<<"\nDense matrix "<<  what<<":\n";
    // printf("Stype %d, Dtype %d, Mtype %d\n", A->Stype,A->Dtype,A->Mtype);
    cout<<"Stype "<<A->Stype<<" Dtype "<<A->Dtype<<" Mtype "<<A->Mtype<<"\n";
    Astore = (DNformat<Field> *) A->Store;
    dp = (typename Field::Element *) Astore->nzval;
    // printf("nrow %d, ncol %d, lda %d\n", A->nrow,A->ncol,Astore->lda);
    cout<<"nrow "<<A->nrow<<" ncol "<<A->ncol<<" lda "<<Astore->lda<<"\n";
    // printf("\nnzval: ");
    cout<<"\nnzval: ";

    // Only the first entry of the solution for Trefethen problem:    
    // cout<<dp[0];

    // Only the entries to check the correctness of the solution for Trefethen     // problem of Matrix (20000x20000)
    // cout<<dp[0]<<"  ";
    // cout<<dp[1]<<"  ";
    // cout<<dp[2]<<"  ";
    // cout<<dp[4]<<"  ";
    // cout<<dp[8]<<"  ";
    // cout<<dp[16]<<"  ";
    // cout<<dp[32]<<"  ";
    // cout<<dp[64]<<"  ";
    // cout<<dp[128]<<"  ";
    // cout<<dp[256]<<"  ";
    // cout<<dp[512]<<"  ";
    // cout<<dp[1024]<<"  ";
    // cout<<dp[2048]<<"  ";
    // cout<<dp[4096]<<"  ";
    // cout<<dp[8192]<<"  ";
    // cout<<dp[16384]<<"  ";
    
    // These are for full solution vector:
    for (i = 0; i < A->nrow; ++i) // printf("%f  ", dp[i]);
      cout<<dp[i]<<"  ";
    
    // printf("\n");
    cout<<"\n";
    fflush(stdout);
}

/*
 * Diagnostic print of column "jcol" in the U/L factor.
 */

template <class Field>
void
print_lu_col(char *msg, int jcol, int pivrow, int *xprune, GlobalLU_t<Field> *Glu, Field& F)
{
    int     i, k, fsupc;
    int     *xsup, *supno;
    int     *xlsub, *lsub;
    typename Field::Element *lusup;
    int     *xlusup;
    typename Field::Element *ucol;
    int     *usub, *xusub;

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    lusup   = Glu->lusup;
    xlusup  = Glu->xlusup;
    ucol    = Glu->ucol;
    usub    = Glu->usub;
    xusub   = Glu->xusub;
    
    cout << msg;
    cout<<"col" << jcol<<" pivrow "<<pivrow<<" supno "<<supno[jcol] <<" xprune "<<xprune[jcol]<<"\n";
    
    cout<<"\tU-col:\n";
    // cout<<"xusub[jcol] "<<xusub[jcol]<<"\n";       For debugging
    // cout<<"xusub[jcol+1] "<<xusub[jcol+1]<<"\n";   For debugging
    // cout<<"\t"<< usub[0] <<"\t"<< ucol[0]<< "\n";  For debugging
    for (i = xusub[jcol]; i < xusub[jcol+1]; i++)
      printf("\t%d%10.4f\n", usub[i], ucol[i]);
      // cout<<"\t"<< usub[i] <<"\t"<< ucol[i]<< "\n";
    // printf("\tL-col in rectangular snode:\n");
    cout<<"\tL-col in rectangular snode:\n";
    fsupc = xsup[supno[jcol]];	/* first col of the snode */
    i = xlsub[fsupc];
    k = xlusup[jcol];
    while ( i < xlsub[fsupc+1] && k < xlusup[jcol+1] ) {
      printf("\t%d\t%10.4f\n", lsub[i], lusup[k]);
      // cout<<"\t"<<lsub[i]<<"\t"<<lusup[k]<<"\n";
      i++; k++;
    }
    fflush(stdout);
}


/*
 * Check whether tempv[] == 0. This should be true before and after 
 * calling any numeric routines, i.e., "panel_bmod" and "column_bmod". 
 */
template <class Field>
void Fcheck_tempv(int n, typename Field::Element *tempv, Field& F)
{
    int i;
	
    for (i = 0; i < n; i++) {
      typename Field::Element zero;
      F.init(zero, 0);
      if (! F.areEqual(tempv[i], zero)) // A.Duran
	{
	  // fprintf(stderr,"tempv[%d] = %f\n", i,tempv[i]);
	  // ABORT("check_tempv");
	  cout << "stderr: tempv[" << i << "] ="<< tempv[i] << endl;
          cout << "Fcheck_tempv\n";
	}
    }
}

template <class Field>
void
GenXtrue(int n, int nrhs, typename Field::Element *x, int ldx, Field& F)
{
    int  i, j;
    for (j = 0; j < nrhs; ++j)
	for (i = 0; i < n; ++i) {
	    x[i + j*ldx] = 1;/* + (double)(i+1.)/n;*/
	}
}

template <class Field>
void
FGenXtrue(int n, int nrhs, typename Field::Element *x, int ldx,Field& F )
{
    int  i, j;
    for (j = 0; j < nrhs; ++j)
	for (i = 0; i < n; ++i) {
	  x[i + j*ldx] = 1;/* + (double)(i+1.)/n;*/ // A.Duran
	}
}


/*
 * Let rhs[i] = sum of i-th row of A, so the solution vector is all 1's
 */
template <class Field>
void
FFillRHS(char *trans, int nrhs, typename Field::Element *x, int ldx,
		SuperMatrix<Field> *A, SuperMatrix<Field> *B, Field& F)
{
    NCformat<Field> *Astore;
    typename Field::Element   *Aval;
    DNformat<Field> *Bstore;
    typename Field::Element   *rhs;
    typename Field::Element one = 1; // A.Duran
    typename Field::Element zero = 0;//A.Duran
    int      ldc;

    Astore = (NCformat<Field> *)A->Store;
    Aval   = (typename Field::Element *) Astore->nzval;
    Bstore = (DNformat<Field> *)B->Store;
    rhs    = ( typename Field::Element  *)Bstore->nzval;
    ldc    = Bstore->lda;
    
    sp_gemm(trans, "N", A->nrow, nrhs, A->ncol, one, A,
	     x, ldx, zero, rhs, ldc, F);

}

/* 
 * Fills a double precision array with a given value.
 */
template <class Field>
void 
Ffill(typename Field::Element *a, int alen, typename Field::Element dval, Field& F)
{
    register int i;
    for (i = 0; i < alen; i++) a[i] = dval;
}

/* 
 * Check the inf-norm of the error vector 
 */
template <class Field>
void dinf_norm_error(int nrhs, SuperMatrix<Field> *X, double *xtrue)
{
    DNformat<Field> *Xstore;
    double err, xnorm;
    double *Xmat, *soln_work;
    int i, j;

    Xstore = (DNformat<Field> *)X->Store;
    Xmat = (double *)Xstore->nzval;

    for (j = 0; j < nrhs; j++) {
      soln_work = &Xmat[j*Xstore->lda];
      err = xnorm = 0.0;
      for (i = 0; i < X->nrow; i++) {
	err = SUPERLU_MAX(err, fabs(soln_work[i] - xtrue[i]));
	xnorm = SUPERLU_MAX(xnorm, fabs(soln_work[i]));
      }
      err = err / xnorm;
      printf("||X - Xtrue||/||X|| = %e\n", err);
    }
}



/* Print performance of the code. */
template <class Field>
void
dPrintPerf(SuperMatrix<Field> *L, SuperMatrix<Field> *U, mem_usage_t *mem_usage,
	       double rpg, double rcond, double *ferr,
	       double *berr, char *equed)
{
    SCformat<Field> *Lstore;
    NCformat<Field> *Ustore;
    extern SuperLUStat_t SuperLUStat;
    double   *utime;
    flops_t  *ops;
    
    utime = SuperLUStat.utime;
    ops   = SuperLUStat.ops;
    
    if ( utime[FACT] != 0. )
	printf("Factor flops = %e\tMflops = %8.2f\n", ops[FACT],
	       ops[FACT]*1e-6/utime[FACT]);
    printf("Identify relaxed snodes	= %8.2f\n", utime[RELAX]);
    if ( utime[SOLVE] != 0. )
	printf("Solve flops = %.0f, Mflops = %8.2f\n", ops[SOLVE],
	       ops[SOLVE]*1e-6/utime[SOLVE]);
    
    Lstore = (SCformat<Field> *) L->Store;
    Ustore = (NCformat<Field> *) U->Store;
    printf("\tNo of nonzeros in factor L = %d\n", Lstore->nnz);
    printf("\tNo of nonzeros in factor U = %d\n", Ustore->nnz);
    printf("\tNo of nonzeros in L+U = %d\n", Lstore->nnz + Ustore->nnz);
	
    printf("L\\U MB %.3f\ttotal MB needed %.3f\texpansions %d\n",
	   mem_usage->for_lu/1e6, mem_usage->total_needed/1e6,
	   mem_usage->expansions);
	
    printf("\tFactor\tMflops\tSolve\tMflops\tEtree\tEquil\tRcond\tRefine\n");
    printf("PERF:%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f\n",
	   utime[FACT], ops[FACT]*1e-6/utime[FACT],
	   utime[SOLVE], ops[SOLVE]*1e-6/utime[SOLVE],
	   utime[ETREE], utime[EQUIL], utime[RCOND], utime[REFINE]);
    
    printf("\tRpg\t\tRcond\t\tFerr\t\tBerr\t\tEquil?\n");
    printf("NUM:\t%e\t%e\t%e\t%e\t%s\n",
	   rpg, rcond, ferr[0], berr[0], equed);
    
}




int print_double_vec(char *what, int n, double *vec)
{
    int i;
    printf("%s: n %d\n", what, n);
    for (i = 0; i < n; ++i) printf("%d\t%f\n", i, vec[i]);
    return 0;
}

