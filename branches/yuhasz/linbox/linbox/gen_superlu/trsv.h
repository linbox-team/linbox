
/*  -- translated by f2c (version 19940927).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "linbox/gen_superlu/f2c.h"
#include "linbox/gen_superlu/xerbla.h"
// #include "linbox/gen_superlu/lsame.c"


// #define DEBUG

/* Subroutine */ 
template <class Field>
int trsv(char *uplo, char *trans, char *diag, int *n, typename Field::Element *a, int *lda, typename Field::Element *x, int *incx, Field& F)
{


    /* System generated locals */
    int a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    static int info;
    static typename Field::Element temp;
    static typename Field::Element temp_mul; // A.Duran
    static int i, j;
    // extern int lsame_(char *, char *); // A.Duran
    static int ix, jx, kx;
    // extern /* Subroutine */ int xerbla_(char *, integer *); // A.Duran
    static int nounit;


/*  Purpose   
    =======   

    DTRSV  solves one of the systems of equations   

       A*x = b,   or   A'*x = b,   

    where b and x are n element vectors and A is an n by n unit, or   
    non-unit, upper or lower triangular matrix.   

    No test for singularity or near-singularity is included in this   
    routine. Such tests must be performed before calling this routine.   

    Parameters   
    ==========   

    UPLO   - CHARACTER*1.   
             On entry, UPLO specifies whether the matrix is an upper or   
             lower triangular matrix as follows:   

                UPLO = 'U' or 'u'   A is an upper triangular matrix.   

                UPLO = 'L' or 'l'   A is a lower triangular matrix.   

             Unchanged on exit.   

    TRANS  - CHARACTER*1.   
             On entry, TRANS specifies the equations to be solved as   
             follows:   

                TRANS = 'N' or 'n'   A*x = b.   

                TRANS = 'T' or 't'   A'*x = b.   

                TRANS = 'C' or 'c'   A'*x = b.   

             Unchanged on exit.   

    DIAG   - CHARACTER*1.   
             On entry, DIAG specifies whether or not A is unit   
             triangular as follows:   

                DIAG = 'U' or 'u'   A is assumed to be unit triangular.   

                DIAG = 'N' or 'n'   A is not assumed to be unit   
                                    triangular.   

             Unchanged on exit.   

    N      - INTEGER.   
             On entry, N specifies the order of the matrix A.   
             N must be at least zero.   
             Unchanged on exit.   

    A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).   
             Before entry with  UPLO = 'U' or 'u', the leading n by n   
             upper triangular part of the array A must contain the upper 
  
             triangular matrix and the strictly lower triangular part of 
  
             A is not referenced.   
             Before entry with UPLO = 'L' or 'l', the leading n by n   
             lower triangular part of the array A must contain the lower 
  
             triangular matrix and the strictly upper triangular part of 
  
             A is not referenced.   
             Note that when  DIAG = 'U' or 'u', the diagonal elements of 
  
             A are not referenced either, but are assumed to be unity.   
             Unchanged on exit.   

    LDA    - INTEGER.   
             On entry, LDA specifies the first dimension of A as declared 
  
             in the calling (sub) program. LDA must be at least   
             max( 1, n ).   
             Unchanged on exit.   

    X      - DOUBLE PRECISION array of dimension at least   
             ( 1 + ( n - 1 )*abs( INCX ) ).   
             Before entry, the incremented array X must contain the n   
             element right-hand side vector b. On exit, X is overwritten 
  
             with the solution vector x.   

    INCX   - INTEGER.   
             On entry, INCX specifies the increment for the elements of   
             X. INCX must not be zero.   
             Unchanged on exit.   


    Level 2 Blas routine.   

    -- Written on 22-October-1986.   
       Jack Dongarra, Argonne National Lab.   
       Jeremy Du Croz, Nag Central Office.   
       Sven Hammarling, Nag Central Office.   
       Richard Hanson, Sandia National Labs.   



       Test the input parameters.   

    
   Parameter adjustments   
       Function Body */
#define X(I) x[(I)-1]

#define A(I,J) a[(I)-1 + ((J)-1)* ( *lda)]

#ifdef DEBUG
    printf("trsv is visited\n");
#endif DEBUG

    info = 0;
    if (! lsame_(uplo, "U") && ! lsame_(uplo, "L")) {
	info = 1;
    } else if (! lsame_(trans, "N") && ! lsame_(trans, "T") &&
	     ! lsame_(trans, "C")) {
	info = 2;
    } else if (! lsame_(diag, "U") && ! lsame_(diag, "N")) {
	info = 3;
    } else if (*n < 0) {
	info = 4;
    } else if (*lda < max(1,*n)) {
	info = 6;
    } else if (*incx == 0) {
	info = 8;
    }
    if (info != 0) {
	xerbla_("DTRSV ", &info);
	return 0;
    }

/*     Quick return if possible. */

    if (*n == 0) {
	return 0;
    }

    nounit = lsame_(diag, "N");

/*     Set up the start point in X if the increment is not unity. This   
       will be  ( N - 1 )*INCX  too small for descending loops. */

    if (*incx <= 0) {
	kx = 1 - (*n - 1) * *incx;
    } else if (*incx != 1) {
	kx = 1;
    }

/*     Start the operations. In this version the elements of A are   
       accessed sequentially with one pass through A. */

    if (lsame_(trans, "N")) {

/*        Form  x := inv( A )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		  // if (X(j) != 0.) {
		  if (!F.isZero(X(j))) { // A. Duran
			if (nounit) {
			  // X(j) /= A(j,j);
			  F.divin( X(j), A(j,j) ); // A. Duran
			}
			// temp = X(j);
			F.assign(temp, X(j) ); // A. Duran
			for (i = j - 1; i >= 1; --i) {
			  // X(i) -= temp * A(i,j);
			  F.subin (X(i), F.mul(temp_mul, temp, A(i,j) ) ); // A. Duran
/* L10: */
			}
		    }
/* L20: */
		}
	    } else {
		jx = kx + (*n - 1) * *incx;
		for (j = *n; j >= 1; --j) {
		  // if (X(jx) != 0.) {
		  if (!F.isZero(X(jx))) { // A. Duran
			if (nounit) {
			  // X(jx) /= A(j,j);
			  F.divin(X(jx), A(j,j)); // A. Duran
			}
			// temp = X(jx);
			F.assign(temp, X(jx)); // A. Duran
			ix = jx;
			for (i = j - 1; i >= 1; --i) {
			    ix -= *incx;
			    // X(ix) -= temp * A(i,j);
			    F.subin(X(ix), F.mul(temp_mul, temp, A(i,j))); // A. Duran
/* L30: */
			}
		    }
		    jx -= *incx;
/* L40: */
		}
	    }
	} else {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		  // if (X(j) != 0.) {
		  if (!F.isZero(X(j))) { // A. Duran
			if (nounit) {
			  // X(j) /= A(j,j);
			  F.divin(X(j), A(j,j)); // A. Duran
			}
			F.assign(temp, X(j));
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			  // X(i) -= temp * A(i,j);
			  F.subin(X(i), F.mul(temp_mul, temp, A(i,j))); // A. Duran
/* L50: */
			}
		    }
/* L60: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		  // if (X(jx) != 0.) {
		  if (!F.isZero(X(jx))) { // A. Duran
			if (nounit) {
			  // X(jx) /= A(j,j);
			  F.divin(X(jx), A(j,j)); // A. Duran
			}
			F.assign(temp, X(jx));
			ix = jx;
			i__2 = *n;
			for (i = j + 1; i <= *n; ++i) {
			    ix += *incx;
			    // X(ix) -= temp * A(i,j);
			    F.subin(X(ix), F.mul(temp_mul, temp, A(i,j))); // A. Duran
/* L70: */
			}
		    }
		    jx += *incx;
/* L80: */
		}
	    }
	}
    } else {

/*        Form  x := inv( A' )*x. */

	if (lsame_(uplo, "U")) {
	    if (*incx == 1) {
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		  // temp = X(j);
		  F.assign(temp, X(j)); // A. Duran
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
		      // temp -= A(i,j) * X(i);
		      F.subin(temp, F.mul(temp_mul, A(i,j), X(i))); // A. Duran
/* L90: */
		    }
		    if (nounit) {
		      // temp /= A(j,j);
		      F.divin(temp, A(j,j)); // A. Duran
		    }
		    // X(j) = temp;
		    F.assign(X(j), temp); // A. Duran
/* L100: */
		}
	    } else {
		jx = kx;
		i__1 = *n;
		for (j = 1; j <= *n; ++j) {
		  // temp = X(jx);
		  F.assign(temp, X(jx)); // A. Duran
		    ix = kx;
		    i__2 = j - 1;
		    for (i = 1; i <= j-1; ++i) {
		      // temp -= A(i,j) * X(ix);
		      F.subin(temp, F.mul(temp_mul, A(i,j), X(ix))); // A. Duran
			ix += *incx;
/* L110: */
		    }
		    if (nounit) {
		      // temp /= A(j,j);
		      F.divin(temp, A(j,j)); // A. Duran
		    }
		    // X(jx) = temp;
		    F.assign(X(jx), temp); // A. Duran
		    jx += *incx;
/* L120: */
		}
	    }
	} else {
	    if (*incx == 1) {
		for (j = *n; j >= 1; --j) {
		  // temp = X(j);
		  F.assign(temp, X(j)); // A. Duran
		    i__1 = j + 1;
		    for (i = *n; i >= j+1; --i) {
			F.subin(temp, F.mul(temp_mul, A(i,j), X(i))); // A. Duran
/* L130: */
		    }
		    if (nounit) {
		      // temp /= A(j,j);
		      F.divin(temp, A(j,j)); // A. Duran
		    }
		    // X(j) = temp;
		    F.assign(X(j), temp); // A. Duran
/* L140: */
		}
	    } else {
		kx += (*n - 1) * *incx;
		jx = kx;
		for (j = *n; j >= 1; --j) {
		  // temp = X(jx);
		  F.assign(temp, X(jx));
		    ix = kx;
		    i__1 = j + 1;
		    for (i = *n; i >= j+1; --i) {
		      // temp -= A(i,j) * X(ix);
		      F.subin(temp, F.mul(temp_mul, A(i,j), X(ix))); // A. Duran
			ix -= *incx;
/* L150: */
		    }
		    if (nounit) {
		      // temp /= A(j,j);
		      F.divin(temp, A(j,j));
		    }
		    // X(jx) = temp;
		    F.assign(X(jx), temp);
		    jx -= *incx;
/* L160: */
		}
	    }
	}
    }

    return 0;

/*     End of TRSV . */

} /* trsv_ */

