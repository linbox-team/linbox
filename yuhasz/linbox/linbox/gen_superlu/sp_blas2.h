
/*
 * -- SuperLU routine (version 2.0) --
 * Univ. of California Berkeley, Xerox Palo Alto Research Center,
 * and Lawrence Berkeley National Lab.
 * November 15, 1997
 *
 */
/*
 * File name:		sp_blas2.c
 * Purpose:		Sparse BLAS 2, using some dense BLAS 2 operations.
 */

#include "linbox/gen_superlu/sp_defs.h"
#include "linbox/gen_superlu/util.h"
#include "linbox/gen_superlu/trsv.h"

// int dtrsv_(char *uplo, char *trans, char *diag, int *n, 
//   double *a, int *lda, double *x, int *incx);

/* 
 * Function prototypes 
 */
template <class Field>
int trsv(char *, char *, char *, int *, typename Field::Element*, int *,  typename Field::Element*, int *, Field& F);

template <class Field>
void usolve(int, int, typename Field::Element*, typename Field::Element*, Field& F);
template <class Field>
void lsolve(int, int, typename Field::Element*, typename Field::Element*, Field& F);
template <class Field>
void matvec(int, int, int, typename Field::Element*, typename Field::Element*, typename Field::Element*, Field& F);


template <class Field>
int
sp_trsv(char *uplo, char *trans, char *diag, SuperMatrix<Field> *L, 
	 SuperMatrix<Field> *U, typename Field::Element *x, int *info, Field& F)
{
/*
 *   Purpose
 *   =======
 *
 *   sp_trsv() solves one of the systems of equations   
 *       A*x = b,   or   A'*x = b,
 *   where b and x are n element vectors and A is a sparse unit , or   
 *   non-unit, upper or lower triangular matrix.   
 *   No test for singularity or near-singularity is included in this   
 *   routine. Such tests must be performed before calling this routine.   
 *
 *   Parameters   
 *   ==========   
 *
 *   uplo   - (input) char*
 *            On entry, uplo specifies whether the matrix is an upper or   
 *             lower triangular matrix as follows:   
 *                uplo = 'U' or 'u'   A is an upper triangular matrix.   
 *                uplo = 'L' or 'l'   A is a lower triangular matrix.   
 *
 *   trans  - (input) char*
 *             On entry, trans specifies the equations to be solved as   
 *             follows:   
 *                trans = 'N' or 'n'   A*x = b.   
 *                trans = 'T' or 't'   A'*x = b.   
 *                trans = 'C' or 'c'   A'*x = b.   
 *
 *   diag   - (input) char*
 *             On entry, diag specifies whether or not A is unit   
 *             triangular as follows:   
 *                diag = 'U' or 'u'   A is assumed to be unit triangular.   
 *                diag = 'N' or 'n'   A is not assumed to be unit   
 *                                    triangular.   
 *	     
 *   L       - (input) SuperMatrix*
 *	       The factor L from the factorization Pr*A*Pc=L*U. Use
 *             compressed row subscripts storage for supernodes,
 *             i.e., L has types: Stype = SC, Dtype = D_, Mtype = TRLU.
 *
 *   U       - (input) SuperMatrix*
 *	        The factor U from the factorization Pr*A*Pc=L*U.
 *	        U has types: Stype = NC, Dtype = D_, Mtype = TRU.
 *    
 *   x       - (input/output) double*
 *             Before entry, the incremented array X must contain the n   
 *             element right-hand side vector b. On exit, X is overwritten 
 *             with the solution vector x.
 *
 *   info    - (output) int*
 *             If *info = -i, the i-th argument had an illegal value.
 *
 */
#ifdef _CRAY
  _fcd ftcs1 = _cptofcd("L", strlen("L")),
    ftcs2 = _cptofcd("N", strlen("N")),
    ftcs3 = _cptofcd("U", strlen("U"));
#endif
  SCformat<Field> *Lstore;
  NCformat<Field> *Ustore;
  typename Field::Element   *Lval, *Uval;
  int incx = 1, incy = 1;
  typename Field::Element alpha, beta; // A.Duran
  F.init(alpha, 1);
  F.init(beta, 1);
  int nrow;
  int fsupc, nsupr, nsupc, luptr, istart, irow;
  int i, k, iptr, jcol;
  typename Field::Element *work;
  typename Field::Element temp_mul; // A.Duran 8/13/2002
  flops_t solve_ops;
  extern SuperLUStat_t SuperLUStat;
  
  /* Test the input parameters */
  *info = 0;
  if ( !lsame_(uplo,"L") && !lsame_(uplo, "U") ) *info = -1;
  else if ( !lsame_(trans, "N") && !lsame_(trans, "T") ) *info = -2;
  else if ( !lsame_(diag, "U") && !lsame_(diag, "N") ) *info = -3;
  else if ( L->nrow != L->ncol || L->nrow < 0 ) *info = -4;
  else if ( U->nrow != U->ncol || U->nrow < 0 ) *info = -5;
  if ( *info ) {
    i = -(*info);
    printf("** On entry to %6s, parameter number %2d had an illegal value\n", "sp_trsv", &i); // Ahmet 8/12/2002
    // xerbla_("sp_trsv", &i);
    return 0;
  }
  
  Lstore = (SCformat<Field> *)L->Store;
  Lval = (typename Field::Element *)Lstore->nzval;
  Ustore = (NCformat<Field> *)U->Store;
  Uval = (typename Field::Element *)Ustore->nzval;
  solve_ops = 0;
  
  if ( !(work = FieldCalloc(L->nrow, F)) )
    ABORT("Malloc fails for work in sp_trsv().");
  
  if ( lsame_(trans, "N") ) {	/* Form x := inv(A)*x. */
    
    if ( lsame_(uplo, "L") ) {
      /* Form x := inv(L)*x */
      if ( L->nrow == 0 ) return 0; /* Quick return */
      
      for (k = 0; k <= Lstore->nsuper; k++) {
	fsupc = L_FST_SUPC(k);
	istart = L_SUB_START(fsupc);
	nsupr = L_SUB_START(fsupc+1) - istart;
	nsupc = L_FST_SUPC(k+1) - fsupc;
	luptr = L_NZ_START(fsupc);
	nrow = nsupr - nsupc;
	
	solve_ops += nsupc * (nsupc - 1);
	solve_ops += 2 * nrow * nsupc;
	
	if ( nsupc == 1 ) {
	  for (iptr=istart+1; iptr < L_SUB_START(fsupc+1); ++iptr) {
	    irow = L_SUB(iptr);
	    ++luptr;
	    // x[irow] -= x[fsupc] * Lval[luptr];
	    F.mul(temp_mul, x[fsupc], Lval[luptr]); // A.Duran 8/30/2002
	    F.subin(x[irow], temp_mul);

	  }
	} else {
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
	  STRSV(ftcs1, ftcs2, ftcs3, &nsupc, &Lval[luptr], &nsupr,
		&x[fsupc], &incx);
	  
	  SGEMV(ftcs2, &nrow, &nsupc, &alpha, &Lval[luptr+nsupc], 
		&nsupr, &x[fsupc], &incx, &beta, &work[0], &incy);
#else
	  trsv("L", "N", "U", &nsupc, &Lval[luptr], &nsupr, &x[fsupc], &incx, F); // A.Duran 8/14/2002
	  
	  dgemv_("N", &nrow, &nsupc, &alpha, &Lval[luptr+nsupc], 
		 &nsupr, &x[fsupc], &incx, &beta, &work[0], &incy);
#endif
#else
	  lsolve ( nsupr, nsupc, &Lval[luptr], &x[fsupc], F);
	  
	  matvec ( nsupr, nsupr-nsupc, nsupc, &Lval[luptr+nsupc],
		    &x[fsupc], &work[0], F);
#endif		
	  
	  iptr = istart + nsupc;
	  for (i = 0; i < nrow; ++i, ++iptr) {
	    irow = L_SUB(iptr);
	    // x[irow] -= work[i];	/* Scatter */
	    F.subin(x[irow], work[i]); // A.Duran 8/13/2002
	    work[i] = 0; // A.Duran
	    
	  }
	}
      } /* for k ... */
      
    } else {
      /* Form x := inv(U)*x */
      
      if ( U->nrow == 0 ) return 0; /* Quick return */
      
      for (k = Lstore->nsuper; k >= 0; k--) {
	fsupc = L_FST_SUPC(k);
	nsupr = L_SUB_START(fsupc+1) - L_SUB_START(fsupc);
	nsupc = L_FST_SUPC(k+1) - fsupc;
	luptr = L_NZ_START(fsupc);
	
	solve_ops += nsupc * (nsupc + 1);
	
	if ( nsupc == 1 ) {
	  // x[fsupc] /= Lval[luptr]; // A.Duran 8/13/2002
	  F.divin(x[fsupc], Lval[luptr]); 
	  for (i = U_NZ_START(fsupc); i < U_NZ_START(fsupc+1); ++i) {
	    irow = U_SUB(i);
	    // x[irow] -= x[fsupc] * Uval[i];
	    F.mul(temp_mul, x[fsupc], Uval[i]); // A.Duran 8/13/2002
	    F.subin(x[irow], temp_mul); // A.Duran 8/13/2002
	  }
	} else {
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
	  STRSV(ftcs3, ftcs2, ftcs2, &nsupc, &Lval[luptr], &nsupr,
		&x[fsupc], &incx);
#else
	  trsv("U", "N", "N", &nsupc, &Lval[luptr], &nsupr, &x[fsupc], &incx, F); // A.Duran 8/14/2002
#endif
#else		
	  usolve ( nsupr, nsupc, &Lval[luptr], &x[fsupc], F);
#endif		
	  
	  for (jcol = fsupc; jcol < L_FST_SUPC(k+1); jcol++) {
	    solve_ops += 2*(U_NZ_START(jcol+1) - U_NZ_START(jcol));
	    for (i = U_NZ_START(jcol); i < U_NZ_START(jcol+1); 
		 i++) {
	      irow = U_SUB(i);
	      // x[irow] -= x[jcol] * Uval[i];
	      F.mul(temp_mul, x[jcol], Uval[i]); // A.Duran 8/13/2002
	      F.subin(x[irow], temp_mul); // A.Duran 8/13/2002
	    }
	  }
	}
      } /* for k ... */
      
    }
  } else { /* Form x := inv(A')*x */
    
    if ( lsame_(uplo, "L") ) {
      /* Form x := inv(L')*x */
      if ( L->nrow == 0 ) return 0; /* Quick return */
      
      for (k = Lstore->nsuper; k >= 0; --k) {
	fsupc = L_FST_SUPC(k);
	istart = L_SUB_START(fsupc);
	nsupr = L_SUB_START(fsupc+1) - istart;
	nsupc = L_FST_SUPC(k+1) - fsupc;
	luptr = L_NZ_START(fsupc);
	
	solve_ops += 2 * (nsupr - nsupc) * nsupc;
	
	for (jcol = fsupc; jcol < L_FST_SUPC(k+1); jcol++) {
	  iptr = istart + nsupc;
	  for (i = L_NZ_START(jcol) + nsupc; 
	       i < L_NZ_START(jcol+1); i++) {
	    irow = L_SUB(iptr);
	    // x[jcol] -= x[irow] * Lval[i];
	    F.mul(temp_mul, x[irow], Lval[i]); // A.Duran 8/13/2002
	    F.subin(x[jcol], temp_mul); // A.Duran 8/13/2002
	    iptr++;
	  }
	}
	
	if ( nsupc > 1 ) {
	  solve_ops += nsupc * (nsupc - 1);
#ifdef _CRAY
	  ftcs1 = _cptofcd("L", strlen("L"));
	  ftcs2 = _cptofcd("T", strlen("T"));
	  ftcs3 = _cptofcd("U", strlen("U"));
	  STRSV(ftcs1, ftcs2, ftcs3, &nsupc, &Lval[luptr], &nsupr,
		&x[fsupc], &incx);
#else
	  trsv("L", "T", "U", &nsupc, &Lval[luptr], &nsupr, &x[fsupc], &incx, F); // A.Duran
#endif
	}
      }
    } else {
      /* Form x := inv(U')*x */
      if ( U->nrow == 0 ) return 0; /* Quick return */
      
      for (k = 0; k <= Lstore->nsuper; k++) {
	fsupc = L_FST_SUPC(k);
	nsupr = L_SUB_START(fsupc+1) - L_SUB_START(fsupc);
	nsupc = L_FST_SUPC(k+1) - fsupc;
	luptr = L_NZ_START(fsupc);
	
	for (jcol = fsupc; jcol < L_FST_SUPC(k+1); jcol++) {
	  solve_ops += 2*(U_NZ_START(jcol+1) - U_NZ_START(jcol));
	  for (i = U_NZ_START(jcol); i < U_NZ_START(jcol+1); i++) {
	    irow = U_SUB(i);
	    // x[jcol] -= x[irow] * Uval[i];
	    F.mul(temp_mul, x[irow], Uval[i]); // A.Duran 8/13/2002
	    F.subin(x[jcol], temp_mul); // A.Duran 8/13/2002
	  }
	}
	
	solve_ops += nsupc * (nsupc + 1);
	
	if ( nsupc == 1 ) {
	  // x[fsupc] /= Lval[luptr]; 
	  F.divin(x[fsupc], Lval[luptr]); // A.Duran 8/13/2002
	} else {
#ifdef _CRAY
	  ftcs1 = _cptofcd("U", strlen("U"));
	  ftcs2 = _cptofcd("T", strlen("T"));
	  ftcs3 = _cptofcd("N", strlen("N"));
	  STRSV( ftcs1, ftcs2, ftcs3, &nsupc, &Lval[luptr], &nsupr,
		 &x[fsupc], &incx);
#else
	  trsv("U", "T", "N", &nsupc, &Lval[luptr], &nsupr, &x[fsupc], &incx, F); // A.Duran
#endif
	}
      } /* for k ... */
    }
  }
  
  SuperLUStat.ops[SOLVE] += solve_ops;
  SUPERLU_FREE(work);
  return 0;
}



// ??
// int
// sp_dgemv(char *trans, double alpha, SuperMatrix *A, double *x, 
//	 int incx, double beta, double *y, int incy)
// {
/*  Purpose   
    =======   

    sp_dgemv()  performs one of the matrix-vector operations   
       y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,   
    where alpha and beta are scalars, x and y are vectors and A is a
    sparse A->nrow by A->ncol matrix.   

    Parameters   
    ==========   

    TRANS  - (input) char*
             On entry, TRANS specifies the operation to be performed as   
             follows:   
                TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.   
                TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.   
                TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.   

    ALPHA  - (input) double
             On entry, ALPHA specifies the scalar alpha.   

    A      - (input) SuperMatrix*
             Matrix A with a sparse format, of dimension (A->nrow, A->ncol).
             Currently, the type of A can be:
                 Stype = NC or NCP; Dtype = D_; Mtype = GE. 
             In the future, more general A can be handled.

    X      - (input) double*, array of DIMENSION at least   
             ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.   
             Before entry, the incremented array X must contain the   
             vector x.   

    INCX   - (input) int
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   

    BETA   - (input) double
             On entry, BETA specifies the scalar beta. When BETA is   
             supplied as zero then Y need not be set on input.   

    Y      - (output) double*,  array of DIMENSION at least   
             ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'   
             and at least   
             ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.   
             Before entry with BETA non-zero, the incremented array Y   
             must contain the vector y. On exit, Y is overwritten by the 
             updated vector y.
	     
    INCY   - (input) int
             On entry, INCY specifies the increment for the elements of   
             Y. INCY must not be zero.   

    ==== Sparse Level 2 Blas routine.   
*/
