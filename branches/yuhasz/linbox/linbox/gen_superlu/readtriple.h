#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#include "linbox/gen_superlu/sp_defs.h"
#include "linbox/gen_superlu/util.h"

template <class Field>
void
readtriple(int *m, int *n, int *nonz,
	    typename Field::Element **nzval, int **rowind, int **colptr, FILE *fp, Field& F )
{
/*
 * Output parameters
 * =================
 *   (a,asub,xa): asub[*] contains the row subscripts of nonzeros
 *	in columns of matrix A; a[*] the numerical values;
 *	row i of A is given by a[k],k=xa[i],...,xa[i+1]-1.
 *
 */
    int    i, j, k, jsize, lasta, nnz, nz;
    typename Field::Element *a, *val;
    int    *asub, *xa, *row, *col;
    
    /* 	Matrix format:
     *    First line:  #rows, #cols, #non-zero
     *    Triplet in the rest of lines:
     *                 row, col, value
     */
    float temp;
    fscanf(fp, "%d%d", n, nonz);
    *m = *n;
    printf("m %d, n %d, nonz %d\n", *m, *n, *nonz);
    FallocateA(*n, *nonz, nzval, rowind, colptr, F); /* Allocate storage */
    a    = *nzval;
    asub = *rowind;
    xa   = *colptr;

    val = (typename Field::Element *) SUPERLU_MALLOC(*nonz * sizeof(typename Field::Element));
    row = (int *) SUPERLU_MALLOC(*nonz * sizeof(int));
    col = (int *) SUPERLU_MALLOC(*nonz * sizeof(int));

    for (j = 0; j < *n; ++j) xa[j] = 0;

    /* Read into the triplet array from a file */
    for (nnz = 0, nz = 0; nnz < *nonz; ++nnz) {
	fscanf(fp, "%d%d%f", &row[nz], &col[nz], &temp);
	F.init(val[nz], (int)temp);
	/* Change to 0-based indexing. */
#if 0
	--row[nz];
	--col[nz];
#endif
	if (row[nz] < 0 || row[nz] >= *m || col[nz] < 0 || col[nz] >= *n
	    /*|| val[nz] == 0.*/) {
// #ifdef NTLZZP
//	    fprintf(stderr, "nz %d, (%d, %d) = %e out of bound, removed\n", 
//		    nz, row[nz], col[nz], rep(val[nz]));
//	    exit(-1);
//#else
	    fprintf(stderr, "nz %d, (%d, %d) =", nz, row[nz], col[nz]);
	     F.write(std::cerr, val[nz]);
	    fprintf(stderr, " out of bound, removed\n");
		
	    exit(-1);
//#endif    
	} else {
	    ++xa[col[nz]];
	    ++nz;
	}
    }

    *nonz = nz;

    /* Initialize the array of column pointers */
    k = 0;
    jsize = xa[0];
    xa[0] = 0;
    for (j = 1; j < *n; ++j) {
	k += jsize;
	jsize = xa[j];
	xa[j] = k;
    }
    
    /* Copy the triplets into the column oriented storage */
    for (nz = 0; nz < *nonz; ++nz) {
	j = col[nz];
	k = xa[j];
	asub[k] = row[nz];
	F.assign(a[k], val[nz]);
	++xa[j];
    }

    /* Reset the column pointers to the beginning of each column */
    for (j = *n; j > 0; --j)
	xa[j] = xa[j-1];
    xa[0] = 0;

    SUPERLU_FREE(val);
    SUPERLU_FREE(row);
    SUPERLU_FREE(col);

#ifdef CHK_INPUT
    for (i = 0; i < *n; i++) {
	printf("Col %d, xa %d\n", i, xa[i]);
	for (k = xa[i]; k < xa[i+1]; k++)
	    printf("%d\t%16.10f\n", asub[k], a[k]);
    }
#endif

}

/*
  void dreadrhs(int m, double *b)
  {
  FILE *fp, *fopen();
  int i, j;
  
  if ( !(fp = fopen("b.dat", "r")) ) {
  fprintf(stderr, "dreadrhs: file does not exist\n");
  exit(-1);
  }
  for (i = 0; i < m; ++i)
  fscanf(fp, "%lf\n", &b[i]);
  /*fscanf(fp, "%d%lf\n", &j, &b[i]);*/
/*        readpair_(j, &b[i]);*/
/* fclose(fp);
   }
*/

