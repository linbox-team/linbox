/* matrix.c
 *
 * Prototypes for the matrix subsystem
 *
 * Written by Bradford Hovinen <hovinen@ximian.com>
 */

#include "matrix.h"
#include "sparse-vector.h"
#include "sparse-compact-vector.h"

#include <stdlib.h>

#define BLOCK_SIZE 8

struct _matrix_t 
{
	vector_class_t  *vector_class;
	unsigned int     rows;
	unsigned int     cols;
	vector_t       **row_data;

	unsigned int     permutation;  /* Nonzero iff this matrix is a
					* permutation matrix */
};

typedef struct _pair_t pair_t;

struct _pair_t 
{
	void *a, *b;
};

static matrix_t *int_matrix_new    (vector_class_t *vector_class,
				    int rows, int cols);
static void matrix_swap_rows       (matrix_t *matrix, int i, int j);
static int matrix_elim_col         (matrix_t *matrix, unsigned int col,
				    unsigned int start_row,
				    matrix_t *p, int permute_only,
				    vector_t *vector);
static int matrix_find_pivot_row   (matrix_t *matrix, int col, int start_row);
static matrix_t *matrix_find_ind_cols (matrix_t *matrix, unsigned int rank);

/* matrix_new
 *
 * Allocate a new matrix of size rows x cols, using the given vector class
 */

matrix_t *
matrix_new (vector_class_t *vector_class, int rows, int cols) 
{
	matrix_t *matrix;
	int i;

	matrix = int_matrix_new (vector_class, rows, cols);

	for (i = 0; i < rows; i++)
		matrix->row_data[i] = vector_new (vector_class, matrix->cols);

	return matrix;
}

/* int_matrix_new
 *
 * Constructs a new matrix but does not allocate vectors
 */

static matrix_t *
int_matrix_new (vector_class_t *vector_class, int rows, int cols)
{
	matrix_t *matrix;

	matrix = (matrix_t *) calloc (sizeof (matrix_t), 1);
	matrix->vector_class = vector_class;
	matrix->rows = rows;
	matrix->cols = cols;

	matrix->row_data = 
		(vector_t **) calloc (sizeof (vector_t *), matrix->rows);

	return matrix;
}

/* matrix_new_identity
 *
 * Creates a new nxn identity matrix
 */

matrix_t *
matrix_new_identity (vector_class_t *vector_class, int n) 
{
	matrix_t *matrix;
	int i;

	matrix = matrix_new (vector_class, n, n);
	matrix->permutation = 1;

	for (i = 0; i < n; i++)
		vector_set_entry (matrix->row_data[i], i, field.one);

	return matrix;
}

/* matrix_clone
 *
 * Creates and returns a copy of the given matrix
 */

matrix_t *
matrix_clone (matrix_t *matrix) 
{
	matrix_t *new_matrix;
	unsigned i;

	new_matrix = int_matrix_new (matrix->vector_class, matrix->rows,
				     matrix->cols);

	for (i = 0; i < matrix->rows; i++)
		new_matrix->row_data[i] = vector_clone (matrix->row_data[i]);

	return new_matrix;
}

/* matrix_destroy
 *
 * Frees all data associated with a matrix
 */

void
matrix_destroy (matrix_t *matrix) 
{
	unsigned i;

	for (i = 0; i < matrix->rows; i++)
		if (matrix->row_data[i])
			vector_unref (matrix->row_data[i]);

	free (matrix->row_data);
	free (matrix);
}

/* matrix_get_rows
 *
 * Returns the number of rows in the matrix
 */

unsigned int
matrix_get_rows (matrix_t *matrix)
{
	return matrix->rows;
}

/* matrix_get_cols
 *
 * Returns the number of columns in the matrix
 */

unsigned int
matrix_get_cols (matrix_t *matrix)
{
	return matrix->cols;
}

/* matrix_change_rep
 *
 * Creates a new matrix with the given underlying vector representation
 */

matrix_t *matrix_change_rep (matrix_t *matrix, vector_class_t *klass)
{
	matrix_t *new_matrix;
	unsigned i;

	new_matrix = int_matrix_new (klass, matrix->rows, matrix->cols);

	for (i = 0; i < matrix->rows; i++)
		new_matrix->row_data[i] = 
			vector_construct (klass, matrix->row_data[i]);

	return new_matrix;
}

/* matrix_get_block_submatrix
 *
 * Construct a block submatrix of a matrix
 */

matrix_t *
matrix_get_block_submatrix (matrix_t *in, 
			    unsigned int r1, unsigned int c1,
			    unsigned int r2, unsigned int c2) 
{
	matrix_t *out;
	unsigned i;

	assert (r1 >= 0);
	assert (c1 >= 0);
	assert (r2 <= in->rows);
	assert (c2 <= in->cols);
	assert (r1 < r2);
	assert (c1 < c2);

	out = int_matrix_new (in->vector_class, r2 - r1, c2 - c1);

	for (i = r1; i < r2; i++)
		out->row_data[i - r1] = 
			vector_get_block_subvector (in->row_data[i], c1, c2);

	return out;
}

/* matrix_mult
 *
 * Multiply two matrices, forming a new matrix as the result
 *
 * Note: This is not meant to be efficient, only elegant
 */

matrix_t *matrix_mult (vector_class_t *vector_class, matrix_t *m1, matrix_t *m2) 
{
	matrix_t *out, *m2t;

	m2t = matrix_transpose (m2);
	out = matrix_mult_transpose (vector_class, m1, m2t);
	matrix_destroy (m2t);

	return out;
}

/* matrix_mult_transpose
 *
 * Finds the product of the matrix m1 and the transpose of the matrix
 * m2. Useful when the transpose of the multiplier is at least as easy to get
 * as the mulitplier itself, e.g. when multiplying a matrix and its own
 * transpose.
 */

matrix_t *matrix_mult_transpose (vector_class_t *vector_class, matrix_t *m1, matrix_t *m2) 
{
	matrix_t *out;
	unsigned i, j;
	Element newv;

	assert (m1->vector_class == m2->vector_class);
	assert (m1->cols == m2->cols);

	out = matrix_new (vector_class, m1->rows, m2->rows);

	for (i = 0; i < m1->rows; i++) {
		for (j = 0; j < m2->rows; j++) {
			vector_dot_product (m1->row_data[i], m2->row_data[j], newv);
			vector_set_entry (out->row_data[i], j, newv);
		}
	}

	return out;
}

typedef struct _left_mult_state_t left_mult_state_t;

struct _left_mult_state_t 
{
	matrix_t *PA, *PAQ;
	unsigned int cur_row;
};

static int
left_mult_cb (vector_t *vector, unsigned int col, Element &element,
	      left_mult_state_t *state) 
{
	if (field.isOne (element)) {
		state->PAQ->row_data[state->cur_row] = state->PA->row_data[col];
		vector_ref (state->PAQ->row_data[state->cur_row]);
	}

	return 1;
}

/* matrix_mult_permutation
 *
 * Computes and returns PAQ, where P and Q are permutation matrices. Either of
 * P or Q may be NULL, in which case they are treated as the identity. P and Q
 * should use the sparse representation for best effect. There is no
 * requirement that A use the same representation.
 */

matrix_t *
matrix_mult_permutation (matrix_t *A, matrix_t *P, matrix_t *Q) 
{
	matrix_t *A_T, *QA_T, *Q_T, *PA, *PAQ;
	left_mult_state_t state;

	assert (P == NULL || P->permutation);
	assert (Q == NULL || Q->permutation);
	assert (P == NULL || P->cols == A->rows);
	assert (Q == NULL || Q->rows == A->cols);

	if (Q != NULL) {
		A_T = matrix_transpose (A);
		Q_T = matrix_transpose (Q);
		QA_T = matrix_mult_permutation (A_T, Q_T, NULL);
		PA = matrix_transpose (QA_T);
		matrix_destroy (QA_T);
		matrix_destroy (A_T);
		matrix_destroy (Q_T);
	} else {
		PA = matrix_clone (A);
	}

	if (P != NULL) {
		state.PA = PA;
		state.PAQ = int_matrix_new (PA->vector_class, PA->rows,
					    PA->cols);

		for (state.cur_row = 0; state.cur_row < A->rows; state.cur_row++)
			vector_foreach_entry (P->row_data[state.cur_row],
					      (entry_cb_t) left_mult_cb, 0,
					      &state);
		matrix_destroy (PA);
		PAQ = state.PAQ;
	} else {
		PAQ = PA;
	}

	return PAQ;
}

/* matrix_mult_permutation_transpose
 *
 * Computes and returns AQ^T, where Q is a permutation matrix. Works best when 
 * Q is represented in the sparse manner
 */

matrix_t *
matrix_mult_permutation_transpose (matrix_t *A, matrix_t *Q) 
{
	matrix_t *A_T, *AQ_T, *AQ;

	assert (Q->permutation);
	assert (Q->cols == A->cols);

	A_T = matrix_transpose (A);
	AQ_T = matrix_mult_permutation (A_T, Q, NULL);
	AQ = matrix_transpose (AQ_T);
	matrix_destroy (AQ_T);

	return AQ;
}

/* matrix_apply_to_vector
 *
 * Apply the given sparse matrix to the given vector
 */

vector_t *matrix_apply_to_vector (matrix_t *matrix, vector_class_t *vector_class, vector_t *vector)
{
	unsigned i;
	vector_t *out;
	Element newv;

	assert (matrix->cols == vector_length (vector));

	out = vector_new (vector_class, matrix->rows);

	for (i = 0; i < matrix->rows; i++) {
		vector_dot_product (vector, matrix->row_data[i], newv);
		vector_set_entry (out, i, newv);
	}

	return out;
}

/* matrix_glue
 *
 * Appends one matrix onto another, forming a new matrix
 */

matrix_t *matrix_glue (matrix_t *m1, matrix_t *m2) 
{
	matrix_t *out;
	unsigned i;

	assert (m1->vector_class == m2->vector_class);
	assert (m1->rows == m2->rows);

	out = matrix_new (m1->vector_class, m1->rows, m1->cols + m2->cols);

	for (i = 0; i < m1->rows; i++)
		out->row_data[i] = 
			vector_glue (m1->row_data[i], m2->row_data[i]);

	return out;
}

/* matrix_inv
 *
 * Computes and returns the inverse of the given matrix. Converts the given
 * matrix to a row-equivalent dialonal matrix
 *
 * Returns the inverse if it exists, NULL otherwise
 */

matrix_t *matrix_inv (matrix_t *matrix, vector_class_t *vector_class)
{
	matrix_t *out;
	unsigned i;
	Element factor;

	assert (matrix->rows == matrix->cols);

	out = matrix_new_identity (vector_class, matrix->rows);
	out->permutation = 0;

	for (i = 0; i < matrix->cols; i++) {
		if (matrix_elim_col (matrix, i, i, out, 0, NULL) < 0) {
			matrix_destroy (out);
			return NULL;
		}
	}

	for (i = 0; i < matrix->cols; i++) {
		vector_get_entry (matrix->row_data[i], i, factor);
		field.inv (factor, factor);
		vector_scale (out->row_data[i], factor);
	}

	return out;
}

/* matrix_row_echelon
 *
 * Eliminates columns to transform the matrix, as much as possible, into row
 * echelon form; performs the same row operations on the accompanying vector
 */

void matrix_row_echelon (matrix_t *matrix, vector_t *vector) 
{
	unsigned int i, s = 0;

	assert (matrix->rows == vector_length (vector));

	for (i = 0; i < matrix->cols; i++)
		if (matrix_elim_col (matrix, i, s, NULL, 0, vector) != -1) s++;
}

/* matrix_nonsingular_minor
 *
 * Finds permutation matrices P and Q such that PAQ has a nonsingular leading
 * principle r x r minor (r <= rank A). Returns PAQ and sets *P and *Q to P
 * and Q respectively
 */

matrix_t *matrix_nonsingular_minor (matrix_t *A, unsigned int r,
				    matrix_t **P, matrix_t **Q)
{
	matrix_t *PAQ, *A_T, *Q_T;

	Q_T = matrix_find_ind_cols (A, r);

	if (matrix_get_rows (A) > r) {
		A_T = matrix_transpose (A);
		*P = matrix_find_ind_cols (A_T, r);
		*Q = matrix_transpose (Q_T);
		matrix_destroy (Q_T);
		PAQ = matrix_mult_permutation (A, *P, *Q);
	} else {
		PAQ = matrix_mult_permutation_transpose (A, Q_T);
		*P = matrix_new_identity (&sparse_compact_vector_class, r);
		*Q = matrix_transpose (Q_T);
		matrix_destroy (Q_T);
	}

	return PAQ;
}

static int transpose_cb (vector_t *v, unsigned int col,
			 Element &element, pair_t *pair) 
{
	matrix_t *out = (matrix_t *) pair->a;
	unsigned int cur_row = (unsigned int) pair->b;

	vector_set_entry (out->row_data[col], cur_row, element);
	return 0;
}

/* matrix_transpose
 *
 * Return the transpose of a given matrix
 */

matrix_t *matrix_transpose (matrix_t *in) 
{
	matrix_t *out;
	pair_t pair;
	unsigned int cur_row;

	out = matrix_new (in->vector_class, in->cols, in->rows);
	pair.a = out;
	if (in->permutation) out->permutation = 1;

	for (cur_row = 0; cur_row < in->rows; cur_row++) {
		pair.b = (void *) cur_row;
		vector_foreach_entry
			(in->row_data[cur_row], (entry_cb_t) transpose_cb, 0, &pair);
	}

	return out;
}

/* matrix_read
 *
 * Read a matrix in Linbox-style sparse format from a stream
 */

matrix_t *matrix_read (vector_class_t *klass, istream &input) 
{
	matrix_t *matrix;
	char c;
	int row = 0, col = 0, rows, cols;
	Element value;

	input >> rows >> cols >> c;

	if (c != 'M') return NULL;

	matrix = matrix_new (klass, rows, cols);

	while (row >= 0 && col >= 0 && input) {
		input >> row >> col;
		field.read (input, value);
		if (row == 0 || col == 0) break;
		vector_set_entry (matrix->row_data[row - 1], col - 1, value);
	}

	return matrix;
}

static int write_cb (vector_t *v, unsigned int col, 
		     Element &value, pair_t *pair) 
{
	ostream *out = (ostream *) pair->a;
	unsigned int cur_row = (unsigned int) pair->b;

	(*out) << cur_row + 1 << " " << col + 1 << " ";
	field.write (*out, value);
	(*out) << endl;
	return 0;
}

/* matrix_write
 *
 * Writes out a matrix in Linbox format to the given stream
 */

void matrix_write (matrix_t *matrix, ostream &out) 
{
	pair_t pair;
	unsigned int cur_row;

	pair.a = &out;
	out << matrix->rows << " " << matrix->cols << "M" << endl;;

	for (cur_row = 0; cur_row < matrix->rows; cur_row++) {
		pair.b = (void *) cur_row;
		vector_foreach_entry (matrix->row_data[cur_row],
				      (entry_cb_t) write_cb, 0,
				      &pair);
	}

	out << "0 0 0" << endl;
}

static int print_cb (vector_t *v, unsigned int col,
		     Element &value) 
{
	field.write (cout, value);
	cout << ", ";
	return 0;
}

/* matrix_pretty_print
 *
 * Prints the matrix onto the console in a human-readable fasion
 */

void matrix_pretty_print (matrix_t *matrix, unsigned int col_width) 
{
	unsigned i;

	for (i = 0; i < matrix->rows; i++) {
//  		cout.width (3);
		cout /* << i */ << " [ ";
		vector_foreach_entry (matrix->row_data[i], 
				      (entry_cb_t) print_cb, 1, NULL);
		cout << "]," << endl;
	}

	cout << endl;
	cout.flush ();
}

/* matrix_get_row_vector
 *
 * Retrieve the vector corresponding to the given row index
 */

vector_t *
matrix_get_row_vector (matrix_t *matrix, unsigned int row)
{
	assert (row < matrix->rows);

	return matrix->row_data[row];
}

/* matrix_equal
 *
 * Determines whether two matrices are equal and, if so, returns nonzero
 */

int
matrix_equal (matrix_t *m1, matrix_t *m2)
{
	unsigned int i;

	if (m1->rows != m2->rows || m1->cols != m2->cols)
		return 0;

	for (i = 0; i < m1->rows; i++) {
		if (!vector_equal (m1->row_data[i], m2->row_data[i]))
			return 0;
	}

	return 1;
}

/* matrix_set_entry
 *
 * Sets a specific (row,col) entry in the matrix to the value given
 */

void
matrix_set_entry (matrix_t *matrix, unsigned int row, unsigned int col,
		  Element &value) 
{
	assert (row < matrix->rows);
	assert (col < matrix->cols);

	vector_set_entry (matrix->row_data[row], col, value);
}

/* matrix_swap_rows
 *
 * Swaps rows i and j in the given matrix
 */

static void
matrix_swap_rows (matrix_t *matrix, int i, int j) 
{
	vector_t *v;

	v = matrix->row_data[i];
	matrix->row_data[i] = matrix->row_data[j];
	matrix->row_data[j] = v;
}

/* matrix_elim_col
 *
 * Fix the matrix through Gaussian elimination so that the leading prinicple
 * (col + 1)x(col + 1) minor is nonsingular
 *
 * Return pivot row on success, -1 if there are no entries in the given column
 */

static int
matrix_elim_col (matrix_t *matrix, unsigned int col, unsigned int start_row,
		 matrix_t *p, int permute_only, vector_t *vector) 
{
	int pivot_row = -1;
	unsigned i;
	Element factor, vinv;
	vector_t *new_row;

	if ((pivot_row = matrix_find_pivot_row (matrix, col, start_row)) == -1)
		return -1;

#ifdef VERBOSE_DEBUG
	printf ("Column: %d, Start row: %d, Pivot row: %d\n",
		col, start_row, pivot_row);
#endif

	if ((unsigned) pivot_row != start_row) {
		matrix_swap_rows (matrix, start_row, pivot_row);

		if (p != NULL)
			matrix_swap_rows (p, start_row, pivot_row);

		if (vector != NULL) {
			Element tmp, tmp1;

			vector_get_entry (vector, start_row, tmp);
			vector_get_entry (vector, pivot_row, tmp1);
			vector_set_entry (vector, start_row, tmp1);
			vector_set_entry (vector, pivot_row, tmp);
		}

		pivot_row = start_row;
	}

	vector_get_entry (matrix->row_data[pivot_row], col, vinv);
	field.inv (vinv, vinv);

	if (p == NULL)
		i = pivot_row + 1;
	else
		i = 0;

	for (; i < matrix->rows; i++) {
		if (i == (unsigned) pivot_row) continue;

		vector_get_entry (matrix->row_data[i], col, factor);
		field.mul (factor, factor, vinv);
		field.neg (factor, factor);

		if (field.isZero (factor)) continue;

#ifdef VERBOSE_DEBUG
		printf ("Need to eliminate row %d by adding factor %d\n",
			i, factor);
#endif

		new_row = vector_row_op (matrix->row_data[i],
					 matrix->row_data[pivot_row],
					 factor);
		vector_unref (matrix->row_data[i]);
		matrix->row_data[i] = new_row;

		if (!permute_only && p != NULL) {
			new_row = vector_row_op (p->row_data[i],
						 p->row_data[pivot_row],
						 factor);
			vector_unref (p->row_data[i]);
			p->row_data[i] = new_row;
		}

		if (vector != NULL) {
			Element tmp, tmp1;

			vector_get_entry (vector, pivot_row, tmp);
			vector_get_entry (vector, i, tmp1);
			field.mulin (tmp, factor);
			field.addin (tmp1, tmp);
			vector_set_entry (vector, i, tmp1);
		}
	}

	return pivot_row;
}

/* matrix_find_pivot_row
 *
 * Find a row index suitable for use as a pivot for the given column in
 * Gaussian elimination; return -1 if it could not be found
 */

static int
matrix_find_pivot_row (matrix_t *matrix, int col, int start_row) 
{
	unsigned i, n, min = matrix->cols;
	int min_idx = -1;
	Element value;

	for (i = start_row; i < matrix->rows; i++) {
		vector_get_entry (matrix->row_data[i], col, value);
		if (field.isZero (value))
			continue;

		n = vector_num_entries (matrix->row_data[i]);

		if (n < min || min_idx == -1) {
			min = n;
			min_idx = i;
		}
	}

	return min_idx;
}

/* matrix_find_ind_cols
 *
 * Find n independent columns in the matrix A and construct a permutation
 * matrix Q such that first n columns of AQ are independent. Returns Q^T
 *
 * This function is optimized for very long matrices, i.e., matrices with a
 * lot of columns, but relatively few rows. It will perform poorly on matrices
 * with many rows.
 */

static matrix_t *
matrix_find_ind_cols (matrix_t *A, unsigned int n) 
{
	matrix_t *A_block, *A_block_p, *tmp1, *tmp2, *Q, *Q_block, *AQ_block;
	int pivot;
	unsigned orig_found, found = 0, col, c = 0;

	Q = matrix_new_identity (&sparse_vector_class, A->cols);
	A_block = matrix_get_block_submatrix (A, 0, 0, A->rows,
					      MIN (BLOCK_SIZE, A->cols));

	while (found < n && c * BLOCK_SIZE < A->cols) {
		A_block_p = matrix_clone (A_block);
		Q_block = matrix_new_identity (&sparse_vector_class,
					       A_block->cols);

		for (col = 0; col < found; col++) {
#ifdef VERBOSE_DEBUG
			/* Debug */
			printf ("(Preprocess) Eliminating column %d, matrix is:\n", col);
			matrix_pretty_print (A_block_p, 2);
#endif
			pivot = matrix_elim_col (A_block_p, col, col, NULL, 0,
						 NULL);
			assert (pivot >= 0);
		}

		orig_found = found;

		for (col = found; col < A_block_p->cols && found < n; col++) {
#ifdef VERBOSE_DEBUG
			/* Debug */
			printf ("Eliminating column %d (%d %d), matrix is:\n",
				col, found, orig_found);
			matrix_pretty_print (A_block_p, 2);
#endif

			pivot = matrix_elim_col (A_block_p, col, found,
						 NULL, 0, NULL);

			if (pivot >= 0) {
				matrix_swap_rows (Q, found,
						  (col - orig_found) +
						  c * BLOCK_SIZE);
				matrix_swap_rows (Q_block, found, col);
				found++;
			}
		}

		c++;

		if (found < n && c * BLOCK_SIZE < A->cols) {
			AQ_block = matrix_mult_permutation_transpose
				(A_block, Q_block);
			tmp1 = matrix_get_block_submatrix (AQ_block, 0, 0, 
							   AQ_block->rows,
							   found);
			tmp2 = matrix_get_block_submatrix (A, 0, 
							   c * BLOCK_SIZE,
							   A->rows, 
							   MIN ((c + 1) *
								BLOCK_SIZE,
								A->cols));
			matrix_destroy (A_block);
			matrix_destroy (Q_block);
			A_block = matrix_glue (tmp1, tmp2);
			matrix_destroy (AQ_block);
			matrix_destroy (tmp1);
			matrix_destroy (tmp2);
		} else {
			matrix_destroy (A_block);
			matrix_destroy (Q_block);
		}
	}

	return Q;
}
