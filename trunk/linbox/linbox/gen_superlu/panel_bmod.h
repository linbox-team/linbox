

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

#include <stdio.h>
#include <stdlib.h>
#include "linbox/gen_superlu/sp_defs.h"
#include "linbox/gen_superlu/sp_ienv.h"
#include "linbox/gen_superlu/util.h"

/* 
 * Function prototypes 
 */
template <class Field>
void lsolve(int, int, typename Field::Element *, typename Field::Element *, Field& F);
template <class Field>
void matvec(int, int, int, typename Field::Element *, typename Field::Element *, typename Field::Element *, Field& F);
extern void Fcheck_tempv(); // ?

template <class Field>
void
panel_bmod (
	    const int  m,          /* in - number of rows in the matrix */
	    const int  w,          /* in */
	    const int  jcol,       /* in */
	    const int  nseg,       /* in */
	    typename Field::Element *dense, /* out, of size n by w */
	    typename Field::Element *tempv, /* working array */
	    int        *segrep,    /* in */
	    int        *repfnz,    /* in, of size n by w */
	    GlobalLU_t<Field> *Glu,        /* modified */
	    Field& F
	    )
{
  /* 
   * Purpose
   * =======
   *
   *    Performs numeric block updates (sup-panel) in topological order.
   *    It features: col-col, 2cols-col, 3cols-col, and sup-col updates.
   *    Special processing on the supernodal portion of L\U[*,j]
   *
   *    Before entering this routine, the original nonzeros in the panel 
   *    were already copied into the spa[m,w].
   *
   *    Updated/Output parameters-
   *	dense[0:m-1,w]: L[*,j:j+w-1] and U[*,j:j+w-1] are returned 
   *      collectively in the m-by-w vector dense[*]. 
   *
   */
  
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
  _fcd ftcs1 = _cptofcd("L", strlen("L")),
    ftcs2 = _cptofcd("N", strlen("N")),
    ftcs3 = _cptofcd("U", strlen("U"));
#endif
  int            incx = 1, incy = 1;
  Field::Element alpha, beta;
#endif
  
  register int            k, ksub;
  int                     fsupc, nsupc, nsupr, nrow;
  int                     krep, krep_ind;
  typename Field::Element ukj, ukj1, ukj2;
  int                     luptr, luptr1, luptr2;
  int                     segsze;
  int                     block_nrow;  /* no of rows in a block row */
  register int            lptr;	 /* Points to the row subscripts of a supernode */
  int                     kfnz, irow, no_zeros; 
  register int            isub, isub1, i;
  register int            jj;	         /* Index through each column in the panel */
  int                     *xsup, *supno;
  int                     *lsub, *xlsub;
  typename Field::Element *lusup;
  int                     *xlusup;
  int                     *repfnz_col; /* repfnz[] for a column in the panel */
  typename Field::Element *dense_col;  /* dense[] for a column in the panel */
  typename Field::Element *tempv1;             /* Used in 1-D update */
  typename Field::Element *TriTmp, *MatvecTmp; /* used in 2-D update */
  typename Field::Element temp_mul;
  typename Field::Element temp_mul1;
  typename Field::Element temp_add;
  typename Field::Element temp_sub;
  typename Field::Element zero;
  F.init(zero, 0);
  typename Field::Element one;
  F.init(one, 1);
  register int            ldaTmp;
  register int            r_ind, r_hi;
  static   int first = 1, maxsuper, rowblk, colblk;
  extern SuperLUStat_t SuperLUStat;
  flops_t  *ops = SuperLUStat.ops;
  
  xsup    = Glu->xsup;
  supno   = Glu->supno;
  lsub    = Glu->lsub;
  xlsub   = Glu->xlsub;
  lusup   = Glu->lusup; // ?
  xlusup  = Glu->xlusup;
  
  if ( first ) {
    maxsuper = sp_ienv(3);
    rowblk   = sp_ienv(4);
    colblk   = sp_ienv(5);
    first = 0;
  }
  ldaTmp = maxsuper + rowblk;
  
  /* 
   * For each nonz supernode segment of U[*,j] in topological order 
   */
  k = nseg - 1;
  for (ksub = 0; ksub < nseg; ksub++) { /* for each updating supernode */
    
    /* krep = representative of current k-th supernode
     * fsupc = first supernodal column
     * nsupc = no of columns in a supernode
     * nsupr = no of rows in a supernode
     */
    krep = segrep[k--];
    fsupc = xsup[supno[krep]];
    nsupc = krep - fsupc + 1;
    nsupr = xlsub[fsupc+1] - xlsub[fsupc];
    nrow = nsupr - nsupc;
    lptr = xlsub[fsupc];
    krep_ind = lptr + nsupc - 1;
    
    repfnz_col = repfnz;
    dense_col = dense; // ?
    
    if ( nsupc >= colblk && nrow > rowblk ) { /* 2-D block update */
      
      TriTmp = tempv;
      
      /* Sequence through each column in panel -- triangular solves */
      for (jj = jcol; jj < jcol + w; jj++,
	     repfnz_col += m, dense_col += m, TriTmp += ldaTmp ) {
	
	kfnz = repfnz_col[krep];
	if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */
	
	segsze = krep - kfnz + 1;
	luptr = xlusup[fsupc];
	
	ops[TRSV] += segsze * (segsze - 1);
	ops[GEMV] += 2 * nrow * segsze;
	
	/* Case 1: Update U-segment of size 1 -- col-col update */
	if ( segsze == 1 ) {
	  F.assign(ukj, dense_col[lsub[krep_ind]]);
	  luptr += nsupr*(nsupc-1) + nsupc;
	  
	  for (i = lptr + nsupc; i < xlsub[fsupc+1]; i++) {
	    irow = lsub[i];
	    // dense_col[irow] -= ukj * lusup[luptr];
	    F.mul(temp_mul, ukj, lusup[luptr]); // A.Duran
	    F.subin(dense_col[irow], temp_mul);
	    ++luptr;
	  }
	  
	} else if ( segsze <= 3 ) {
	  F.assign(ukj, dense_col[lsub[krep_ind]]);
	  F.assign(ukj1, dense_col[lsub[krep_ind - 1]]);
	  luptr += nsupr*(nsupc-1) + nsupc-1;
	  luptr1 = luptr - nsupr;
	  
	  if ( segsze == 2 ) {
	    // ukj -= ukj1 * lusup[luptr1];
	    F.mul(temp_mul, ukj1, lusup[luptr1]); // A.Duran
	    F.subin(ukj, temp_mul);
	    F.assign(dense_col[lsub[krep_ind]], ukj);
	    for (i = lptr + nsupc; i < xlsub[fsupc+1]; ++i) {
	      irow = lsub[i];
	      luptr++; luptr1++;
	      // dense_col[irow] -= (ukj*lusup[luptr]
	      //	  + ukj1*lusup[luptr1]);
	      F.mul(temp_mul, ukj, lusup[luptr]); // A.Duran
	      F.mul(temp_mul1,ukj1, lusup[luptr1]);
	      F.add(temp_add, temp_mul, temp_mul1);
	      F.subin(dense_col[irow], temp_add);
	    }
	  } else {
	    F.assign(ukj2, dense_col[lsub[krep_ind - 2]]);
	    luptr2 = luptr1 - nsupr;
	    // ukj1 -= ukj2 * lusup[luptr2-1];
	    F.mul(temp_mul, ukj2, lusup[luptr2-1]); // A.Duran
	    F.subin(ukj1, temp_mul);
	    // ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
	    F.mul(temp_mul, ukj1, lusup[luptr1]); // A.Duran
	    F.sub(temp_sub, ukj, temp_mul);
	    F.mul(temp_mul, ukj2,lusup[luptr2]);
	    F.sub(ukj, temp_sub, temp_mul);
	    F.assign(dense_col[lsub[krep_ind]], ukj);
	    F.assign(dense_col[lsub[krep_ind-1]], ukj1);
	    for (i = lptr + nsupc; i < xlsub[fsupc+1]; ++i) {
	      irow = lsub[i];
	      luptr++; luptr1++; luptr2++;
	      // dense_col[irow] -= ( ukj*lusup[luptr]
	      //	   + ukj1*lusup[luptr1] + ukj2*lusup[luptr2] );
	      F.mul(temp_mul, ukj, lusup[luptr]); // A.Duran
	      F.mul(temp_mul1, ukj1, lusup[luptr1]);
	      F.add(temp_add, temp_mul, temp_mul1);
	      F.mul(temp_mul, ukj2, lusup[luptr2]);
	      F.addin(temp_add, temp_mul);
	      F.subin(dense_col[irow], temp_add);
	    }
	  }
	  
	} else  {	/* segsze >= 4 */
	  
	  /* Copy U[*,j] segment from dense[*] to TriTmp[*], which
	     holds the result of triangular solves.    */
	  no_zeros = kfnz - fsupc;
	  isub = lptr + no_zeros;
	  for (i = 0; i < segsze; ++i) {
	    irow = lsub[isub];
	    F.assign(TriTmp[i], dense_col[irow]); /* Gather */
	    ++isub;
	  }
	  
	  /* start effective triangle */
	  luptr += nsupr * no_zeros + no_zeros;
	  
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
	  STRSV( ftcs1, ftcs2, ftcs3, &segsze, &lusup[luptr], 
		 &nsupr, TriTmp, &incx );
#else
	  dtrsv_( "L", "N", "U", &segsze, &lusup[luptr], 
		  &nsupr, TriTmp, &incx );
#endif
#else		
	  lsolve ( nsupr, segsze, &lusup[luptr], TriTmp, F );
#endif
	  
	  
	} /* else ... */
	
      }  /* for jj ... end tri-solves */
      
      /* Block row updates; push all the way into dense[*] block */
      for ( r_ind = 0; r_ind < nrow; r_ind += rowblk ) {
	
	r_hi = SUPERLU_MIN(nrow, r_ind + rowblk);
	block_nrow = SUPERLU_MIN(rowblk, r_hi - r_ind);
	luptr = xlusup[fsupc] + nsupc + r_ind;
	isub1 = lptr + nsupc + r_ind;
	
	repfnz_col = repfnz;
	TriTmp = tempv;
	dense_col = dense; // ?
	
	/* Sequence through each column in panel -- matrix-vector */
	for (jj = jcol; jj < jcol + w; jj++,
	       repfnz_col += m, dense_col += m, TriTmp += ldaTmp) {
	  
	  kfnz = repfnz_col[krep];
	  if ( kfnz == EMPTY ) continue; /* Skip any zero segment */
	  
	  segsze = krep - kfnz + 1;
	  if ( segsze <= 3 ) continue;   /* skip unrolled cases */
	  
	  /* Perform a block update, and scatter the result of
	     matrix-vector to dense[].		 */
	  no_zeros = kfnz - fsupc;
	  luptr1 = luptr + nsupr * no_zeros;
	  MatvecTmp = &TriTmp[maxsuper];
	  
#ifdef USE_VENDOR_BLAS
	  F.assign(alpha, one); 
	  F.assign(beta, zero);
#ifdef _CRAY
	  SGEMV(ftcs2, &block_nrow, &segsze, &alpha, &lusup[luptr1], 
		&nsupr, TriTmp, &incx, &beta, MatvecTmp, &incy);
#else
	  dgemv_("N", &block_nrow, &segsze, &alpha, &lusup[luptr1], 
		 &nsupr, TriTmp, &incx, &beta, MatvecTmp, &incy);
#endif
#else
	  matvec(nsupr, block_nrow, segsze, &lusup[luptr1], TriTmp, MatvecTmp, F);
#endif
	  
	  /* Scatter MatvecTmp[*] into SPA dense[*] temporarily
	   * such that MatvecTmp[*] can be re-used for the
	   * the next blok row update. dense[] will be copied into 
	   * global store after the whole panel has been finished.
	   */
	  isub = isub1;
	  for (i = 0; i < block_nrow; i++) {
	    irow = lsub[isub];
	    // dense_col[irow] -= MatvecTmp[i];
	    F.subin(dense_col[irow], MatvecTmp[i]); // A.Duran
	    F.assign(MatvecTmp[i], zero);
	    ++isub;
	  }
	  
	} /* for jj ... */
	
      } /* for each block row ... */
      
      /* Scatter the triangular solves into SPA dense[*] */
      repfnz_col = repfnz;
      TriTmp = tempv;
      dense_col = dense;
      
      for (jj = jcol; jj < jcol + w; jj++,
	     repfnz_col += m, dense_col += m, TriTmp += ldaTmp) {
	kfnz = repfnz_col[krep];
	if ( kfnz == EMPTY ) continue; /* Skip any zero segment */
	
	segsze = krep - kfnz + 1;
	if ( segsze <= 3 ) continue; /* skip unrolled cases */
	
	no_zeros = kfnz - fsupc;		
	isub = lptr + no_zeros;
	for (i = 0; i < segsze; i++) {
	  irow = lsub[isub];
	  F.assign(dense_col[irow], TriTmp[i]);
	  F.assign(TriTmp[i], zero);
	  ++isub;
	}
	
      } /* for jj ... */
      
    } else { /* 1-D block modification */
      
      
      /* Sequence through each column in the panel */
      for (jj = jcol; jj < jcol + w; jj++,
	     repfnz_col += m, dense_col += m) {
	
	kfnz = repfnz_col[krep];
	if ( kfnz == EMPTY ) continue;	/* Skip any zero segment */
	
	segsze = krep - kfnz + 1;
	luptr = xlusup[fsupc];
	
	ops[TRSV] += segsze * (segsze - 1);
	ops[GEMV] += 2 * nrow * segsze;
	
	/* Case 1: Update U-segment of size 1 -- col-col update */
	if ( segsze == 1 ) {
	  F.assign(ukj, dense_col[lsub[krep_ind]]);
	  luptr += nsupr*(nsupc-1) + nsupc;
	  
	  for (i = lptr + nsupc; i < xlsub[fsupc+1]; i++) {
	    irow = lsub[i];
	    // dense_col[irow] -= ukj * lusup[luptr];
	    F.mul(temp_mul, ukj, lusup[luptr]); // A.Duran
	    F.subin(dense_col[irow], temp_mul);
	    ++luptr;
	  }
	  
	} else if ( segsze <= 3 ) {
	  F.assign(ukj, dense_col[lsub[krep_ind]]);
	  luptr += nsupr*(nsupc-1) + nsupc-1;
	  F.assign(ukj1, dense_col[lsub[krep_ind - 1]]);
	  luptr1 = luptr - nsupr;
	  
	  if ( segsze == 2 ) {
	    // ukj -= ukj1 * lusup[luptr1];
	    F.mul(temp_mul, ukj1, lusup[luptr1]); // A.Duran
	    F.subin(ukj, temp_mul);
	    F.assign(dense_col[lsub[krep_ind]], ukj);
	    for (i = lptr + nsupc; i < xlsub[fsupc+1]; ++i) {
	      irow = lsub[i];
	      ++luptr;  ++luptr1;
	      // dense_col[irow] -= (ukj*lusup[luptr]
	      //	  + ukj1*lusup[luptr1]);
	      F.mul(temp_mul, ukj, lusup[luptr]); // A.Duran
	      F.mul(temp_mul1, ukj1, lusup[luptr1]);
	      F.add(temp_add, temp_mul, temp_mul1);
	      F.subin(dense_col[irow], temp_add);
	    }
	  } else {
	    F.assign(ukj2, dense_col[lsub[krep_ind - 2]]);
	    luptr2 = luptr1 - nsupr;
	    // ukj1 -= ukj2 * lusup[luptr2-1];
	    F.mul(temp_mul, ukj2, lusup[luptr2-1]); // A.Duran
	    F.subin(ukj1, temp_mul);
	    // ukj = ukj - ukj1*lusup[luptr1] - ukj2*lusup[luptr2];
	    F.mul(temp_mul, ukj1, lusup[luptr1]); // A.Duran
	    F.sub(temp_sub, ukj, temp_mul);
	    F.mul(temp_mul, ukj2, lusup[luptr2]);
	    F.sub(ukj, temp_sub, temp_mul);
	    F.assign(dense_col[lsub[krep_ind]], ukj);
	    F.assign(dense_col[lsub[krep_ind-1]], ukj1);
	    for (i = lptr + nsupc; i < xlsub[fsupc+1]; ++i) {
	      irow = lsub[i];
	      ++luptr; ++luptr1; ++luptr2;
	      // dense_col[irow] -= ( ukj*lusup[luptr]
	      //	   + ukj1*lusup[luptr1] + ukj2*lusup[luptr2] );
	      F.mul(temp_mul, ukj, lusup[luptr]);
	      F.mul(temp_mul1, ukj1, lusup[luptr1]);
	      F.add(temp_add, temp_mul, temp_mul1);
	      F.mul(temp_mul, ukj2, lusup[luptr2] );
	      F.addin(temp_add, temp_mul);
	      F.subin(dense_col[irow], temp_add);
	    }
	  }
	  
	} else  { /* segsze >= 4 */
	  /* 
	   * Perform a triangular solve and block update,
	   * then scatter the result of sup-col update to dense[].
	   */
	  no_zeros = kfnz - fsupc;
	  
	  /* Copy U[*,j] segment from dense[*] to tempv[*]: 
	   *    The result of triangular solve is in tempv[*];
	   *    The result of matrix vector update is in dense_col[*]
	   */
	  isub = lptr + no_zeros;
	  for (i = 0; i < segsze; ++i) {
	    irow = lsub[isub];
	    F.assign(tempv[i], dense_col[irow]); /* Gather */
	    ++isub;
	  }
	  
	  /* start effective triangle */
	  luptr += nsupr * no_zeros + no_zeros;
	  
#ifdef USE_VENDOR_BLAS
#ifdef _CRAY
	  STRSV( ftcs1, ftcs2, ftcs3, &segsze, &lusup[luptr], 
		 &nsupr, tempv, &incx );
#else
	  dtrsv_( "L", "N", "U", &segsze, &lusup[luptr], 
		  &nsupr, tempv, &incx );
#endif
	  
	  luptr += segsze;	/* Dense matrix-vector */
	  tempv1 = &tempv[segsze]; // ?
	  F.assign(alpha, one);
	  F.assign(beta, zero);
#ifdef _CRAY
	  SGEMV( ftcs2, &nrow, &segsze, &alpha, &lusup[luptr], 
		 &nsupr, tempv, &incx, &beta, tempv1, &incy );
#else
	  dgemv_( "N", &nrow, &segsze, &alpha, &lusup[luptr], 
		  &nsupr, tempv, &incx, &beta, tempv1, &incy );
#endif
#else
	  lsolve ( nsupr, segsze, &lusup[luptr], tempv, F );
	  
	  luptr += segsze;        /* Dense matrix-vector */
	  tempv1 = &tempv[segsze];
	  matvec (nsupr, nrow, segsze, &lusup[luptr], tempv, tempv1, F);
#endif
	  
	  /* Scatter tempv[*] into SPA dense[*] temporarily, such
	   * that tempv[*] can be used for the triangular solve of
	   * the next column of the panel. They will be copied into 
	   * ucol[*] after the whole panel has been finished.
	   */
	  isub = lptr + no_zeros;
	  for (i = 0; i < segsze; i++) {
	    irow = lsub[isub];
	    F.assign(dense_col[irow], tempv[i]);
	    F.assign(tempv[i], zero);
	    isub++;
	  }
	  
	  /* Scatter the update from tempv1[*] into SPA dense[*] */
	  /* Start dense rectangular L */
	  for (i = 0; i < nrow; i++) {
	    irow = lsub[isub];
	    // dense_col[irow] -= tempv1[i];
	    F.subin(dense_col[irow], tempv1[i]); // A.Duran
	    F.assign(tempv1[i], zero);
	    ++isub;	
	  }
	  
	} /* else segsze>=4 ... */
	
      } /* for each column in the panel... */
      
    } /* else 1-D update ... */
    
  } /* for each updating supernode ... */
  
}
