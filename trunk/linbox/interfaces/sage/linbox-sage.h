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

#include<stddef.h>

#include <cstdlib>
#include <vector>

/*****************************************************************

  Dense over Z/nZ.
 
*****************************************************************/

// Element could be either double or float

template<class Element>
unsigned long linbox_modn_dense_echelonize (Element modulus, Element* matrix,
				  size_t nrows, size_t ncols);

template <class Element>
void linbox_modn_dense_minpoly (Element modulus, Element * &mp, size_t& degree, 
				size_t n, Element * matrix);
template <class Element>
void linbox_modn_dense_charpoly (Element modulus, Element * &cp, size_t n, Element * matrix);

template <class Element>
void linbox_modn_dense_delete_array (Element *f);

template <class Element>
void linbox_modn_dense_matrix_matrix_multiply (Element modulus, Element * ans,
					      Element *A, Element *B,
					      size_t m, size_t n, size_t k); 

template <class Element>
void linbox_modn_dense_matrix_matrix_general_multiply (Element modulus,
						       Element * ans,
						       Element alpha, Element beta,
						       Element *A, Element *B,
						       size_t m, size_t n, size_t k); 
template <class Element>
unsigned long linbox_modn_dense_rank (Element modulus, Element* matrix, size_t nrows, size_t ncols);

template<class Element>
Element linbox_modn_dense_det (Element modulus, Element* matrix, size_t nrows, size_t ncols);

template<class Element>
unsigned long linbox_modn_dense_col_rankprofile_submatrix (Element modulus,
							   Element* matrix,
							   Element* outmatrix,
							   size_t& rank,
							   size_t nrows, size_t ncols);

template <class Element>
unsigned long linbox_modn_dense_col_rankprofile_submatrix_indices (Element modulus,
								   Element* matrix,
								   size_t *& row_idx,
								   size_t *& col_idx,
								   size_t & rank,
								   size_t nrows,
								   size_t ncols);

/*****************************************************************

  Dense over ZZ
 
*****************************************************************/
 
/* linbox_minpoly allocates space for minpoly, so you have to call linbox_delete_array
   to free it up afterwards. */

void linbox_integer_dense_minpoly (mpz_t*& minpoly, size_t& degree, 
				   size_t n, mpz_t* matrix);

void linbox_integer_dense_charpoly (mpz_t*& minpoly, size_t& degree, 
				    size_t n, mpz_t* matrix);

void linbox_integer_dense_delete_array (mpz_t* f);

/* ans must be a pre-allocated and pre-initialized array of GMP ints. */
int linbox_integer_dense_matrix_matrix_multiply (mpz_t* ans, mpz_t *A, mpz_t *B,
						 size_t m, size_t n, size_t k);

unsigned long linbox_integer_dense_rank(mpz_t* matrix, size_t nrows, size_t ncols);

void linbox_integer_dense_det(mpz_t ans, mpz_t* matrix, size_t nrows, size_t ncols);

void linbox_integer_dense_smithform(mpz_t *&v, mpz_t *matrix, size_t nrows, size_t ncols);

void linbox_integer_dense_double_det (mpz_t ans1, mpz_t ans2, mpz_t **a,
				      mpz_t ** b, mpz_t **c, size_t n, int proof);

/*****************************************************************

  Sparse over Z/nZ
 
*****************************************************************/

template <class Element>
unsigned long linbox_modn_sparse_matrix_rank(Element modulus, 
					     size_t numrows, 
					     size_t numcols, 
					     void *rows,
					     int reorder);
template<class Element>
std::vector<Element> linbox_modn_sparse_matrix_solve(Element modulus, 
						     size_t numrows, 
						     size_t numcols,  
						     void *a, 
						     void *b,
						     int method);

#endif // __LINBOX_SAGE_H
