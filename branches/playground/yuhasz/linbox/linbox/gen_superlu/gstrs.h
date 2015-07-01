

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

#include "linbox/gen_superlu/sp_defs.h"
#include "linbox/gen_superlu/util.h"
#include "linbox/gen_superlu/sp_blas2.h"
#include "linbox/gen_superlu/myblas2.h"

// #define DEBUG
/* 
 * Function prototypes 
 */

template <class Field>
void usolve(int, int, typename Field::Element*, typename Field::Element*, Field& F);
template <class Field>
void lsolve(int, int, typename Field::Element*, typename Field::Element*, Field& F);
template <class Field>
void matvec(int, int, int, typename Field::Element*, typename Field::Element*, typename Field::Element*, Field& F);

template <class Field>
void
gstrs (char *trans, SuperMatrix<Field> *L, SuperMatrix<Field> *U,
	int *perm_r, int *perm_c, SuperMatrix<Field> *B, int *info, Field& F)
{
/*
 * Purpose
 * =======
 *
 * GSTRS solves a system of linear equations A*X=B or A'*X=B
 * with A sparse and B dense, using the LU factorization computed by
 * GSTRF.
 *
 * See supermatrix.h for the definition of 'SuperMatrix' structure.
 *
 * Arguments
 * =========
 *
 * trans   (input) char*
 *          Specifies the form of the system of equations:
 *          = 'N':  A * X = B  (No transpose)
 *          = 'T':  A'* X = B  (Transpose)
 *          = 'C':  A**H * X = B  (Conjugate transpose)
 *
 * L       (input) SuperMatrix*
 *         The factor L from the factorization Pr*A*Pc=L*U as computed by
 *         dgstrf(). Use compressed row subscripts storage for supernodes,
 *         i.e., L has types: Stype = SC, Dtype = D_, Mtype = TRLU.
 *
 * U       (input) SuperMatrix*
 *         The factor U from the factorization Pr*A*Pc=L*U as computed by
 *         dgstrf(). Use column-wise storage scheme, i.e., U has types:
 *         Stype = NC, Dtype = D_, Mtype = TRU.
 *
 * perm_r  (input) int*, dimension (L->nrow)
 *         Row permutation vector, which defines the permutation matrix Pr; 
 *         perm_r[i] = j means row i of A is in position j in Pr*A.
 *
 * perm_c  (input) int*, dimension (L->ncol)
 *	   Column permutation vector, which defines the 
 *         permutation matrix Pc; perm_c[i] = j means column i of A is 
 *         in position j in A*Pc.
 *
 * B       (input/output) SuperMatrix*
 *         B has types: Stype = DN, Dtype = D_, Mtype = GE.
 *         On entry, the right hand side matrix.
 *         On exit, the solution matrix if info = 0;
 *
 * info    (output) int*
 * 	   = 0: successful exit
 *	   < 0: if info = -i, the i-th argument had an illegal value
 *
 */
#ifdef _CRAY
    _fcd ftcs1, ftcs2, ftcs3, ftcs4;
#endif
    int      incx = 1, incy = 1;
    typename Field::Element   alpha, beta; // A.Duran
    F.init(alpha, 1);
    F.init(beta, 1);
    DNformat<Field> *Bstore;
    typename Field::Element *Bmat;
    SCformat<Field> *Lstore;
    NCformat<Field> *Ustore;
    typename Field::Element   *Lval, *Uval;
    int      nrow, notran;
    int      fsupc, nsupr, nsupc, luptr, istart, irow;
    int      i, j, k, iptr, jcol, n, ldb, nrhs;
    typename Field::Element   *work, *work_col, *rhs_work, *soln;
    typename Field::Element   temp_mul; // A.Duran 8/8/2002
    flops_t  solve_ops;
    extern SuperLUStat_t SuperLUStat;
    // template <class Field>
    void dprint_soln();

    /* Test input parameters ... */
    *info = 0;
    Bstore = (DNformat<Field> *)B->Store;
    ldb = Bstore->lda;
    nrhs = B->ncol;
    notran = lsame_(trans, "N");
    if ( !notran && !lsame_(trans, "T") && !lsame_(trans, "C") ) *info = -1;
    else if ( L->nrow != L->ncol || L->nrow < 0 ||
	      L->Stype != SC || L->Dtype != D_ || L->Mtype != TRLU )
	*info = -2;
    else if ( U->nrow != U->ncol || U->nrow < 0 ||
	      U->Stype != NC || U->Dtype != D_ || U->Mtype != TRU )
	*info = -3;
    else if ( ldb < SUPERLU_MAX(0, L->nrow) ||
	      B->Stype != DN || B->Dtype != D_ || B->Mtype != GE )
	*info = -6;
    if ( *info ) {
	i = -(*info);
printf("** On entry to %6s, parameter number %2d had an illegal value\n",
		"dgstrs", &i);
// xerbla_("dgstrs", &i);
	return;
    }

    n = L->nrow;
    work = FieldCalloc(n * nrhs, F);
    if ( !work ) ABORT("Malloc fails for local work[].");
    soln = FieldMalloc(n, F);
    if ( !soln ) ABORT("Malloc fails for local soln[].");

    Bmat = (typename Field::Element *)Bstore->nzval;
    Lstore = (SCformat<Field> *)L->Store;
    Lval = (typename Field::Element *)Lstore->nzval;
    Ustore = (NCformat<Field> *)U->Store;
    Uval = (typename Field::Element *)Ustore->nzval;
    solve_ops = 0;
    
    if ( notran ) {
	/* Permute right hand sides to form Pr*B */
	for (i = 0; i < nrhs; i++) {
	    rhs_work = &Bmat[i*ldb];
	    for (k = 0; k < n; k++) soln[perm_r[k]] = rhs_work[k];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}
	
	/* Forward solve PLy=Pb. */
	for (k = 0; k <= Lstore->nsuper; k++) {
	    fsupc = L_FST_SUPC(k);
	    istart = L_SUB_START(fsupc);
	    nsupr = L_SUB_START(fsupc+1) - istart;
	    nsupc = L_FST_SUPC(k+1) - fsupc;
	    nrow = nsupr - nsupc;

	    solve_ops += nsupc * (nsupc - 1) * nrhs;
	    solve_ops += 2 * nrow * nsupc * nrhs;
	    
	    if ( nsupc == 1 ) {
		for (j = 0; j < nrhs; j++) {
		    rhs_work = &Bmat[j*ldb];
	    	    luptr = L_NZ_START(fsupc);
		    for (iptr=istart+1; iptr < L_SUB_START(fsupc+1); iptr++){
			irow = L_SUB(iptr);
			++luptr;
			/* rhs_work[irow] -= rhs_work[fsupc] * Lval[luptr]; */
			// A.Duran 8/8/2002
			F.mul(temp_mul, rhs_work[fsupc], Lval[luptr]);
			F.subin(rhs_work[irow],temp_mul); 

		    }
		}
	    } else {
	    	luptr = L_NZ_START(fsupc);
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
		ftcs1 = _cptofcd("L", strlen("L"));
		ftcs2 = _cptofcd("N", strlen("N"));
		ftcs3 = _cptofcd("U", strlen("U"));
		STRSM( ftcs1, ftcs1, ftcs2, ftcs3, &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
		
		SGEMM( ftcs2, ftcs2, &nrow, &nrhs, &nsupc, &alpha, 
			&Lval[luptr+nsupc], &nsupr, &Bmat[fsupc], &ldb, 
			&beta, &work[0], &n );
#else
		dtrsm_("L", "L", "N", "U", &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
		
		dgemm_( "N", "N", &nrow, &nrhs, &nsupc, &alpha, 
			&Lval[luptr+nsupc], &nsupr, &Bmat[fsupc], &ldb, 
			&beta, &work[0], &n );
#endif
		for (j = 0; j < nrhs; j++) {
		    rhs_work = &Bmat[j*ldb];
		    work_col = &work[j*n];
		    iptr = istart + nsupc;
		    for (i = 0; i < nrow; i++) {
			irow = L_SUB(iptr);
			rhs_work[irow] -= work_col[i]; /* Scatter */
			work_col[i] = 0.0;
			iptr++;
		    }
		}
#else		
		for (j = 0; j < nrhs; j++) {
		    rhs_work = &Bmat[j*ldb];
		    lsolve (nsupr, nsupc, &Lval[luptr], &rhs_work[fsupc], F);
		    matvec (nsupr, nrow, nsupc, &Lval[luptr+nsupc],
			    &rhs_work[fsupc], &work[0], F );

		    iptr = istart + nsupc;
		    for (i = 0; i < nrow; i++) {
			irow = L_SUB(iptr);
			// rhs_work[irow] -= work[i];
			F.subin(rhs_work[irow], work[i]); // A.Duran 8/13/2002
			work[i] = 0; // A.Duran
			iptr++;
		    }
		}
#endif		    
	    } /* else ... */
	} /* for L-solve */

#ifdef DEBUG
  	cout<<"After L-solve: y=\n";
	print_soln(n, nrhs, Bmat, F);
#endif

	/*
	 * Back solve Ux=y.
	 */
	for (k = Lstore->nsuper; k >= 0; k--) {
	    fsupc = L_FST_SUPC(k);
	    istart = L_SUB_START(fsupc);
	    nsupr = L_SUB_START(fsupc+1) - istart;
	    nsupc = L_FST_SUPC(k+1) - fsupc;
	    luptr = L_NZ_START(fsupc);

	    solve_ops += nsupc * (nsupc + 1) * nrhs;

	    if ( nsupc == 1 ) {
		rhs_work = &Bmat[0];
		for (j = 0; j < nrhs; j++) {
		  cout << "At gstrs.h before F.divin\n";
		  // rhs_work[fsupc] /= Lval[luptr];
		  F.divin(rhs_work[fsupc], Lval[luptr]); // A.Duran 8/13/2002
		  cout << "At gstrs.h after F.divin\n";
		  rhs_work += ldb;
		}
	    } else {
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
		ftcs1 = _cptofcd("L", strlen("L"));
		ftcs2 = _cptofcd("U", strlen("U"));
		ftcs3 = _cptofcd("N", strlen("N"));
		STRSM( ftcs1, ftcs2, ftcs3, ftcs3, &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
#else
		dtrsm_("L", "U", "N", "N", &nsupc, &nrhs, &alpha,
		       &Lval[luptr], &nsupr, &Bmat[fsupc], &ldb);
#endif
#else		
		for (j = 0; j < nrhs; j++)
		    usolve ( nsupr, nsupc, &Lval[luptr], &Bmat[fsupc+j*ldb], F );
#endif		
	    }

	    for (j = 0; j < nrhs; ++j) {
		rhs_work = &Bmat[j*ldb];
		for (jcol = fsupc; jcol < fsupc + nsupc; jcol++) {
		    solve_ops += 2*(U_NZ_START(jcol+1) - U_NZ_START(jcol));
		    for (i = U_NZ_START(jcol); i < U_NZ_START(jcol+1); i++ ){
			irow = U_SUB(i);
			// rhs_work[irow] -= rhs_work[jcol] * Uval[i];
			F.mul(temp_mul, rhs_work[jcol], Uval[i]); // A.Duran 8/13/2002
			F.subin(rhs_work[irow], temp_mul); // A.Duran 8/13/2002
		    }
		}
	    }
	    
	} /* for U-solve */

#ifdef DEBUG
  	printf("After U-solve: x=\n");
	print_soln(n, nrhs, Bmat, F);
#endif

	/* Compute the final solution X := Pc*X. */
	for (i = 0; i < nrhs; i++) {
	    rhs_work = &Bmat[i*ldb];
	    for (k = 0; k < n; k++) soln[k] = rhs_work[perm_c[k]];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}
	
        SuperLUStat.ops[SOLVE] = solve_ops;

    } else { /* Solve A'*X=B */
	/* Permute right hand sides to form Pc'*B. */
	for (i = 0; i < nrhs; i++) {
	    rhs_work = &Bmat[i*ldb];
	    for (k = 0; k < n; k++) soln[perm_c[k]] = rhs_work[k];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}

	SuperLUStat.ops[SOLVE] = 0;
	
	for (k = 0; k < nrhs; ++k) {
	    
	    /* Multiply by inv(U'). */
	    sp_trsv("U", "T", "N", L, U, &Bmat[k*ldb], info, F);
	    
	    /* Multiply by inv(L'). */
	    sp_trsv("L", "T", "U", L, U, &Bmat[k*ldb], info, F);
	    
	}
	
	/* Compute the final solution X := Pr'*X (=inv(Pr)*X) */
	for (i = 0; i < nrhs; i++) {
	    rhs_work = &Bmat[i*ldb];
	    for (k = 0; k < n; k++) soln[k] = rhs_work[perm_r[k]];
	    for (k = 0; k < n; k++) rhs_work[k] = soln[k];
	}

    }

    SUPERLU_FREE(work);
    SUPERLU_FREE(soln);
}

/*
 * Diagnostic print of the solution vector 
 */
template <class Field>
void
print_soln(int n, int nrhs, typename Field::Element *soln, Field& F)
{
    int i;

    for (i = 0; i < n; i++) 
  	cout<<"\t"<<i<<": "<<soln[i]<<"\n";
}
