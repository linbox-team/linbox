/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/* linbox-sage.h
 * Copyright (C) 2007 Martin Albrecht
 *               2008 Clement Pernet
 *
 * Written by Martin Albrecht
 *            Clement Pernet <clement.pernet@gmail.com>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_SAGE_H
#define __LINBOX_SAGE_H
#ifdef __cplusplus
#define EXTERN extern "C"
#else 
#define EXTERN
#endif

#include<stddef.h>
typedef size_t mod_int; 

#include <cstdlib>
#include <vector>

/*****************************************************************

  Dense over Z/nZ.
 
*****************************************************************/

EXTERN int linbox_modn_dense_echelonize(mod_int modulus,  
					mod_int** matrix, size_t nrows, size_t ncols);


EXTERN void linbox_modn_dense_minpoly(mod_int modulus, mod_int **mp, size_t* degree, 
				      size_t n, mod_int **matrix, int do_minpoly);

EXTERN void linbox_modn_dense_delete_array(mod_int *f);

EXTERN void linbox_modn_dense_delete_dbl_array(double *f);

EXTERN int linbox_modn_dense_matrix_matrix_multiply(mod_int modulus, mod_int **ans,
						    mod_int **A, mod_int **B,
						    size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc); 

EXTERN int linbox_modn_dense_rank(mod_int modulus,  
				  mod_int** matrix, size_t nrows, size_t ncols);

EXTERN mod_int linbox_modn_dense_det(mod_int modulus,  
				     mod_int** matrix, size_t nrows, size_t ncols);

EXTERN int linbox_modn_dense_col_rankprofile_submatrix (mod_int modulus,
							mod_int** matrix,
							double ** outmatrix,
							size_t* rank,
							size_t nrows, size_t ncols);
EXTERN int linbox_modn_dense_col_rankprofile_submatrix_indices (mod_int modulus,
								mod_int** matrix,
								size_t ** row_idx,
								size_t ** col_idx,
								size_t * rank,
								size_t nrows,
								size_t ncols);
/*****************************************************************

  Dense over ZZ
 
*****************************************************************/
 
/* linbox_minpoly allocates space for minpoly, so you have to call linbox_delete_array
   to free it up afterwards. */
EXTERN void linbox_integer_dense_minpoly_hacked(mpz_t** minpoly, size_t* degree, 
                  size_t n, mpz_t** matrix, int do_minpoly);
EXTERN void linbox_integer_dense_minpoly(mpz_t** minpoly, size_t* degree, 
                  size_t n, mpz_t** matrix);
EXTERN void linbox_integer_dense_charpoly(mpz_t** charpoly, size_t* degree, 
                  size_t n, mpz_t** matrix);
EXTERN void linbox_integer_dense_delete_array(mpz_t* f);

/* ans must be a pre-allocated and pre-initialized array of GMP ints. */
EXTERN int linbox_integer_dense_matrix_matrix_multiply(mpz_t** ans, mpz_t **A, mpz_t **B,
			      size_t A_nr, size_t A_nc, size_t B_nr, size_t B_nc);

EXTERN unsigned long linbox_integer_dense_rank(mpz_t** matrix, size_t nrows,
					       size_t ncols);

EXTERN  void linbox_integer_dense_det(mpz_t ans, mpz_t** matrix, size_t nrows,
				      size_t ncols);

EXTERN void linbox_integer_dense_smithform(mpz_t **v, 
					   mpz_t **matrix, 
					   size_t nrows, size_t ncols);

EXTERN void linbox_integer_dense_double_det (mpz_t ans1, mpz_t ans2, mpz_t **a,
					     mpz_t ** b, mpz_t **c, size_t n, int proof);
/*****************************************************************

  Sparse over Z/nZ
 
*****************************************************************/


EXTERN unsigned long linbox_modn_sparse_matrix_rank(mod_int modulus, 
						    size_t numrows, 
						    size_t numcols, 
						    void *rows,
						    int reorder);

EXTERN std::vector<unsigned int> linbox_modn_sparse_matrix_solve(mod_int modulus, 
								 size_t numrows, 
								 size_t numcols,  
								 void *a, 
								 void *b,
								 int method);

#endif // __LINBOX_SAGE_H
