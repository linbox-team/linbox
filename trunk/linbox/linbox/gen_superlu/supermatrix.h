#ifndef __SUPERLU_SUPERMATRIX /* allow multiple inclusions */
#define __SUPERLU_SUPERMATRIX
#include <stdio.h>
#include <stdlib.h>
#include "linbox/gen_superlu/sp_defs.h"

/********************************************
 * The matrix types are defined as follows. *
 ********************************************/
typedef enum {
  NC,        /* column-wise, no supernode */
  NR,        /* row-wize, no supernode */
  SC,        /* column-wise, supernode */
  SR,        /* row-wise, supernode */
  NCP,       /* column-wise, column-permuted, no supernode 
		(The consecutive columns of nonzeros, after permutation,
		may not be stored  contiguously.) */
  DN         /* Fortran style column-wise storage for dense matrix */
} Stype_t;

typedef enum {
  S_,         /* single */
  D_,         /* double */
  C_,         /* single complex */
  Z_          /* double complex */
} Dtype_t;

typedef enum {
  GE,        /* general */
  TRLU,      /* lower triangular, unit diagonal */
  TRUU,      /* upper triangular, unit diagonal */
  TRL,       /* lower triangular */
  TRU,       /* upper triangular */
  SYL,       /* symmetric, store lower half */
  SYU,       /* symmetric, store upper half */
  HEL,       /* Hermitian, store lower half */
  HEU        /* Hermitian, store upper half */
} Mtype_t;

template <class Field>
/*typedef*/ struct SuperMatrix {
  Stype_t Stype; /* Storage type: interprets the storage structure 
		    pointed to by *Store. */
  Dtype_t Dtype; /* Data type. */
  Mtype_t Mtype; /* Matrix type: describes the mathematical property of 
		    the matrix. */
  int  nrow;     /* number of rows */
  int  ncol;     /* number of columns */
  /* pointer to the actual storage of the matrix */
  typename Field::Element *Store;   
} ;//SuperMatrix;

/***********************************************
 * The storage schemes are defined as follows. *
 ***********************************************/

/* Stype == NC (Also known as Harwell-Boeing sparse matrix format (CCS)) */
template<class Field>
/*typedef*/ struct NCformat{
  int  nnz;	  /* number of nonzeros in the matrix */
  /* pointer to array of nonzero values, packed by column */
  typename Field::Element *nzval;  
  int  *rowind; /* pointer to array of row indices of the nonzeros */
  int  *colptr; /* pointer to array of beginning of columns in nzval[] 
		   and rowind[]  */
  /* Note:
     Zero-based indexing is used;
     colptr[] has ncol+1 entries, the last one pointing
     beyond the last column, so that colptr[ncol] = nnz. */
} ;//NCformat;

/* Stype == NR (Also known as row compressed storage (RCS). */
template<class Field>
/*typedef*/ struct NRformat{
  int  nnz;	  /* number of nonzeros in the matrix */
  /* pointer to array of nonzero values, packed by row */
  typename Field::Element *nzval;  
  int  *colind; /* pointer to array of column indices of the nonzeros */
  int  *rowptr; /* pointer to array of beginning of rows in nzval[] 
		   and colind[]  */
  /* Note:
     Zero-based indexing is used;
     rowptr[] has nrow+1 entries, the last one pointing
     beyond the last column, so that rowptr[nrow] = nnz. */
} ;//NRformat;

/* Stype == SC */
template<class Field>
/*typedef*/ struct SCformat{
  int  nnz;	     /* number of nonzeros in the matrix */
  int  nsuper;       /* number of supernodes, minus 1 */
  /* pointer to array of nonzero values, packed by column */
  typename Field::Element *nzval; 
  int  *nzval_colptr;/* pointer to array of beginning of columns in nzval[] */
  int  *rowind;      /* pointer to array of compressed row indices of 
			rectangular supernodes */
  int *rowind_colptr;/* pointer to array of beginning of columns in rowind[] */
  int *col_to_sup;   /* col_to_sup[j] is the supernode number to which column 
			j belongs; mapping from column to supernode number. */
  int *sup_to_col;   /* sup_to_col[s] points to the start of the s-th 
			supernode; mapping from supernode number to column.
		        e.g.: col_to_sup: 0 1 2 2 3 3 3 4 4 4 4 4 (ncol=12)
			sup_to_col: 0 1 2 4 7 12            (nsuper=4) */
  /* Note:
     Zero-based indexing is used;
     nzval_colptr[], rowind_colptr[], col_to_sup and
     sup_to_col[] have ncol+1 entries, the last one
     pointing beyond the last column.         */
} ;//SCformat;

/* Stype == NCP */
template<class Field>
/*typedef*/ struct NCPformat{
  int nnz;	  /* number of nonzeros in the matrix */
  /* pointer to array of nonzero values, packed by column */
  typename Field::Element *nzval;  
  int *rowind;  /* pointer to array of row indices of the nonzeros */
  /* Note: nzval[]/rowind[] always have the same length */
  int *colbeg;  /* colbeg[j] points to the beginning of column j in nzval[] 
		   and rowind[]  */
  int *colend;  /* colend[j] points to one past the last element of column
		   j in nzval[] and rowind[]  */
  /* Note:
     Zero-based indexing is used;
     The consecutive columns of the nonzeros may not be 
     contiguous in storage, because the matrix has been 
     postmultiplied by a column permutation matrix. */
} ;//NCPformat;

/* Stype == DN */
template<class Field>
/*typedef*/ struct DNformat {
  int lda;      /* leading dimension */
  /* array of size lda*ncol to represent a dense matrix */
  typename Field::Element *nzval;  
} ;//DNformat;



/*********************************************************
 * Macros used for easy access of sparse matrix entries. *
 *********************************************************/
#define L_SUB_START(col)     ( Lstore->rowind_colptr[col] )
#define L_SUB(ptr)           ( Lstore->rowind[ptr] )
#define L_NZ_START(col)      ( Lstore->nzval_colptr[col] )
#define L_FST_SUPC(superno)  ( Lstore->sup_to_col[superno] )
#define U_NZ_START(col)      ( Ustore->colptr[col] )
#define U_SUB(ptr)           ( Ustore->rowind[ptr] )


#endif  /* __SUPERLU_SUPERMATRIX */
