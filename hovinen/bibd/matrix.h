/* matrix.h
 *
 * Prototypes for the matrix subsystem
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 */

#ifndef __MATRIX_H
#define __MATRIX_H

#include "vector.h"

typedef struct _matrix_t matrix_t;

matrix_t *matrix_new                 (vector_class_t *vector_class,
				      int rows, int cols);
matrix_t *matrix_new_identity        (vector_class_t *vector_class, int n);
matrix_t *matrix_clone               (matrix_t *matrix);
void      matrix_destroy             (matrix_t *matrix);

unsigned int matrix_get_rows         (matrix_t *matrix);
unsigned int matrix_get_cols         (matrix_t *matrix);

matrix_t *matrix_change_rep          (matrix_t *matrix, vector_class_t *klass);

matrix_t *matrix_get_block_submatrix (matrix_t *in, 
				      unsigned int r1, unsigned int c1,
				      unsigned int r2, unsigned int c2);

matrix_t *matrix_mult                (vector_class_t *vector_class, matrix_t *m1, matrix_t *m2);
matrix_t *matrix_mult_transpose      (vector_class_t *vector_class, matrix_t *m1, matrix_t *m2);
matrix_t *matrix_mult_permutation    (matrix_t *A, matrix_t *P, matrix_t *Q);
matrix_t *matrix_mult_permutation_transpose (matrix_t *A, matrix_t *Q);
vector_t *matrix_apply_to_vector     (matrix_t *in, vector_class_t *vector_class, vector_t *vector);

matrix_t *matrix_glue                (matrix_t *m1, matrix_t *m2);
matrix_t *matrix_inv                 (matrix_t *matrix, vector_class_t *vector_class);
void      matrix_row_echelon         (matrix_t *matrix, vector_t *vector);
matrix_t *matrix_nonsingular_minor   (matrix_t *in, unsigned int rank,
				      matrix_t **p, matrix_t **q);
matrix_t *matrix_transpose           (matrix_t *matrix);

matrix_t *matrix_read                (vector_class_t *klass, istream &input);
void      matrix_write               (matrix_t *matrix, ostream &output);
void      matrix_pretty_print        (matrix_t *matrix,
				      unsigned int col_width);

vector_t *matrix_get_row_vector      (matrix_t *matrix, unsigned int row);

int       matrix_equal               (matrix_t *m1, matrix_t *m2);

void      matrix_set_entry           (matrix_t *matrix, unsigned int row,
				      unsigned int col, Element &value);

#endif /* __MATRIX_H */
