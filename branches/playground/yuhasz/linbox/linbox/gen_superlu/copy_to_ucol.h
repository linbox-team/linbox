

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

// #define DEBUG // A. Duran

template <class Field>
int
Fcopy_to_ucol(
	      int        jcol,	  /* in */
	      int        nseg,	  /* in */
	      int        *segrep,  /* in */
	      int        *repfnz,  /* in */
	      int        *perm_r,  /* in */
	      typename Field::Element *dense,   /* modified - reset to zero on return */
	      GlobalLU_t<Field> *Glu,      /* modified */
	      Field &F
	      )
{
/* 
 * Gather from SPA dense[*] to global ucol[*].
 */
    int ksub, krep, ksupno;
    int i, k, kfnz, segsze;
    int fsupc, isub, irow;
    int jsupno, nextu;
    int new_next, mem_error;
    int       *xsup, *supno;
    int       *lsub, *xlsub;
    typename Field::Element *ucol;
    int       *usub, *xusub;
    int       nzumax;

    typename Field::Element zero; // A.Duran
    F.init(zero,0);

    xsup    = Glu->xsup;
    supno   = Glu->supno;
    lsub    = Glu->lsub;
    xlsub   = Glu->xlsub;
    ucol    = Glu->ucol;
    usub    = Glu->usub;
    xusub   = Glu->xusub;
    nzumax  = Glu->nzumax;
    
    jsupno = supno[jcol];
    nextu  = xusub[jcol];
    k = nseg - 1;
    for (ksub = 0; ksub < nseg; ksub++) {
	krep = segrep[k--];
	ksupno = supno[krep];

	if ( ksupno != jsupno ) { /* Should go into ucol[] */
	    kfnz = repfnz[krep];
	    if ( kfnz != EMPTY ) {	/* Nonzero U-segment */

	    	fsupc = xsup[ksupno];
	        isub = xlsub[fsupc] + kfnz - fsupc;
	        segsze = krep - kfnz + 1;

		new_next = nextu + segsze;
		while ( new_next > nzumax ) {
		    if (mem_error = FLUMemXpand(jcol, nextu, UCOL, &nzumax, Glu, F))
			return (mem_error);
		    ucol = Glu->ucol;
		    if (mem_error = FLUMemXpand(jcol, nextu, USUB, &nzumax, Glu, F))
			return (mem_error);
		    usub = Glu->usub;
		    lsub = Glu->lsub;
		}
		
		for (i = 0; i < segsze; i++) {
		    irow = lsub[isub];
		    usub[nextu] = perm_r[irow];
		    // F.assign(ucol[nextu], dense[irow]);
		    ucol[nextu] = dense[irow];
		    F.assign(dense[irow], zero);
#ifdef DEBUG
		    cout << "At copy_to_ucol  ucol[nextu] "<<ucol[nextu] <<"\n";
#endif

		    nextu++;
		    isub++;
		} 

	    }

	}

    } /* for each segment... */

    xusub[jcol + 1] = nextu;      /* Close U[*,jcol] */
    return 0;
}

