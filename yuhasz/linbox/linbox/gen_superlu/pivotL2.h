

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
#include <stdlib.h>
#include "linbox/gen_superlu/sp_defs.h"
#include "linbox/gen_superlu/util.h"

#define MIN_COL 0 // A.Duran

// #define DEBUGPERMR
// #define DEBUG // A.Duran

template <class Field>
int
pivotL(
        const int  jcol,     /* in */
        const typename Field::Element u, /* in - diagonal pivoting threshold */
        int        *usepr,   /* re-use the pivot sequence given by perm_r/iperm_r */
        int        *perm_r,  /* may be modified */
        int        *iperm_r, /* in - inverse of perm_r */
        int        *iperm_c, /* in - used to find diagonal of Pc*A*Pc' */
        int        *pivrow,  /* out */
        GlobalLU_t<Field> *Glu,      /* modified - global LU data structures */
	int        *nullity,  // nullity of matrix A  
	Field& F
       )
{
/*
 * Purpose
 * =======
 *   Performs the numerical pivoting on the current column of L,
 *   and the CDIV operation.
 *
 *   Pivot policy:
 *   (1) Compute thresh = u * max_(i>=j) abs(A_ij);
 *   (2) IF user specifies pivot row k and abs(A_kj) >= thresh THEN
 *           pivot row = k;
 *       ELSE IF abs(A_jj) >= thresh THEN
 *           pivot row = j;
 *       ELSE
 *           pivot row = m;
 * 
 *   Note: If you absolutely want to use a given pivot order, then set u=0.0.
 *
 *   Return value: 0      success;
 *                 i > 0  U(i,i) is exactly zero.
 *
 */
  int          fsupc;	/* first column in the supernode */
  int          nsupc;	/* no of columns in the supernode */
  int          nsupr;   /* no of rows in the supernode */
  int          lptr;	/* points to the starting subscript of the supernode */
  int          pivptr, old_pivptr, diag, diagind;
  typename Field::Element pivmax, rtemp, thresh; // ?? A.Duran
  typename Field::Element temp;
  typename Field::Element *lu_sup_ptr; 
  typename Field::Element *lu_col_ptr;
  int          *lsub_ptr;
  int          isub, icol, k, itemp;
  int          *lsub, *xlsub;
  typename Field::Element *lusup;
  int          *xlusup;
  extern SuperLUStat_t SuperLUStat;
  flops_t  *ops = SuperLUStat.ops;
  typename Field::Element zero; 
  F.init(zero, 0);
  typename Field::Element one; 
  F.init(one, 1);
  int found;
  /* Initialize pointers */
  lsub       = Glu->lsub;
  xlsub      = Glu->xlsub;
  lusup      = Glu->lusup;
  xlusup     = Glu->xlusup;
  fsupc      = (Glu->xsup)[(Glu->supno)[jcol]];
  nsupc      = jcol - fsupc;	        /* excluding jcol; nsupc >= 0 */
  lptr       = xlsub[fsupc];
  nsupr      = xlsub[fsupc+1] - lptr;
  lu_sup_ptr = &lusup[xlusup[fsupc]];	/* start of the current supernode */
  lu_col_ptr = &lusup[xlusup[jcol]];	/* start of jcol in the supernode */
  lsub_ptr   = &lsub[lptr];	/* start of row indices of the supernode */
  
#ifdef DEBUG
  if ( jcol == MIN_COL ) {
    printf("Before cdiv: col %d\n", jcol);
    for (k = nsupc; k < nsupr; k++) 
      printf("  lu[%d] %f\n", lsub_ptr[k], lu_col_ptr[k]);
  }
#endif

  // At genlinsol approach, if the field is NOT unparametric float or double,
  // we don't determine the largest abs numerical value for partial pivoting.
  // Instead, we pick only the first encountered nonzero value. A. Duran
 
  /* Determine the largest abs numerical value for partial pivoting;
     Also search for user-specified pivot, and diagonal element. */
  if ( *usepr ) *pivrow = iperm_r[jcol];
  diagind = iperm_c[jcol];
  F.init(pivmax, 0);
  pivptr = nsupc;
  diag = EMPTY;
  old_pivptr = nsupc;
#ifndef UNPARAMETRICFIELDFD // A. Duran added
  found = 0;
  for (isub = nsupc; isub < nsupr; ++isub) {
    F.assign(rtemp, (typename Field::Element)labs(lu_col_ptr[isub]));
    // if (!F.areEqual(rtemp, zero)) {
    if (rtemp > pivmax) {
      F.assign(pivmax, rtemp);
      pivptr = isub;
      // cout << rtemp << "visited\n";
      // break;
      // found = 1;
    }
    if ( *usepr && lsub_ptr[isub] == *pivrow ) old_pivptr = isub;
    if ( lsub_ptr[isub] == diagind ) diag = isub;
    // if (found) break;
  }
#else  
  for (isub = nsupc; isub < nsupr; ++isub) {
    F.assign(rtemp, (typename Field::Element)fabs (lu_col_ptr[isub]));
    if ( rtemp > pivmax ) { // ? Greater(rtemp, pivmax, F)
      F.assign(pivmax, rtemp);
      pivptr = isub;
      // cout << "visited\n" << rtemp;
    }
    if ( *usepr && lsub_ptr[isub] == *pivrow ) old_pivptr = isub;
    if ( lsub_ptr[isub] == diagind ) diag = isub;
  }
#endif
  
  /* Test for singularity */
  if ( F.areEqual(pivmax, zero) ) {
    *pivrow = lsub_ptr[pivptr];
    perm_r[*pivrow] = jcol;  // A. Duran
    // perm_r[pivptr] = jcol; // A. Duran added
    
#ifdef DEBUGPERMR
    cout<< "At pivotL .. test for singularity.. jcol "<< jcol <<"\n";
    cout<< "pivmax "<<pivmax<<"\n"; 
    cout<< "pivptr "<<pivptr<<"\n";
#endif
    ++(*nullity); // A. Duran
    *usepr = 0;
    return (jcol+1);
  }
  
  // thresh = u * pivmax;
  F.mul(thresh, u, pivmax);  // A.Duran 8/15/2002
  
  /* Choose appropriate pivotal element by our policy. */
  if ( *usepr ) {
    // do I need fabs or abs for Modular field ? I think no.
#ifndef UNPARAMETRICFIELDFD // A. Duran added
    F.assign(rtemp, (typename Field::Element)labs((long int)lu_col_ptr[old_pivptr]));
    if ( !F.areEqual(rtemp, zero)) // ?
      pivptr = old_pivptr;
    else
      *usepr = 0;
#else
    F.assign(rtemp, (typename Field::Element) fabs (lu_col_ptr[old_pivptr]));
    
    if ( (!F.areEqual(rtemp, zero)) && rtemp >= thresh ) // ?
      pivptr = old_pivptr;
    else
      *usepr = 0;
#endif
  }
  if ( *usepr == 0 ) {
    /* Use diagonal pivot? */
    if ( diag >= 0 ) { /* diagonal exists */
#ifndef UNPARAMETRICFIELDFD // A. Duran added
      F.assign(rtemp, (typename Field::Element)labs((long int)lu_col_ptr[diag]));
      if (!F.areEqual(rtemp, zero)) pivptr = diag; // ?
#else
      F.assign(rtemp, (typename Field::Element)fabs (lu_col_ptr[diag]));
      if ( (!F.areEqual(rtemp, zero)) && rtemp >= thresh ) pivptr = diag; // ?
#endif
    }
    *pivrow = lsub_ptr[pivptr];
  }
  
  /* Record pivot row */
  perm_r[*pivrow] = jcol;
  
  /* Interchange row subscripts */
  if ( pivptr != nsupc ) {
    itemp = lsub_ptr[pivptr];
    lsub_ptr[pivptr] = lsub_ptr[nsupc];
    lsub_ptr[nsupc] = itemp;
    
    /* Interchange numerical values as well, for the whole snode, such 
     * that L is indexed the same way as A.
     */
    for (icol = 0; icol <= nsupc; icol++) {
      itemp = pivptr + icol * nsupr;
      F.assign(temp, lu_sup_ptr[itemp]);
      F.assign(lu_sup_ptr[itemp], lu_sup_ptr[nsupc + icol*nsupr]);
      F.assign(lu_sup_ptr[nsupc + icol*nsupr], temp);
    }
  } /* if */
  
  /* cdiv operation */
  ops[FACT] += nsupr - nsupc;
  
  // temp = 1.0 / lu_col_ptr[nsupc];
  F.div(temp, one, lu_col_ptr[nsupc]); // A.Duran
  
  for (k = nsupc+1; k < nsupr; k++) 
    // lu_col_ptr[k] *= temp;
    F.mulin(lu_col_ptr[k], temp);
  
  
#ifdef DEBUG
  if ( jcol == MIN_COL ) {
    printf("At the end of pivotL: col %d\n", jcol);
    for (k = nsupc; k < nsupr; k++) 
      printf("  lu[%d] %f\n", lsub_ptr[k], lu_col_ptr[k]);
  }
#endif
  
  return 0;
}

  
