/* bibd-matrix.h
 * Written by Bradford Hovinen <hovinen@ximian.com>
 * Create a matrix representing a BIBD problem
 */

#ifndef __BIBD_MATRIX_H
#define __BIBD_MATRIX_H

#include "matrix.h"
#include "vector.h"

matrix_t *make_bibd_matrix (int v, int k);
vector_t *make_bibd_vector (int v, int k, int lambda);
matrix_t *refine_bibd_matrix (matrix_t *bibd_matrix,
			      unsigned int *includes_list, int includes_len,
			      unsigned int *excludes_list, int excludes_len);
vector_t *refine_bibd_vector (vector_t *bibd_vector,
			      unsigned int *list, unsigned int len,
			      matrix_t *bibd_matrix);
unsigned int get_subset (unsigned int column);
unsigned int get_column_number (unsigned int subset);
unsigned int *get_subset_array (unsigned int n, unsigned int r, unsigned int *len);

/* Debugging only */
void print_subsets (int n, int r);

#endif /* __BIBD_MATRIX_H */
