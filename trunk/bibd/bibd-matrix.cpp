/* bibd-matrix.cpp
 * Written by Bradford Hovinen <hovinen@ximian.com>
 * Create a matrix representing a BIBD problem
 */

#include <iostream>

#include "bibd-matrix.h"
#include "sparse-compact-vector.h"
#include "sparse-vector.h"
#include "dense-vector.h"

typedef struct _bibd_state_t bibd_state_t;

struct _bibd_state_t 
{
	int v, k;
	int row, col;
	int two_elem_subset;
	matrix_t *bibd_matrix;
	Element one;
};

typedef void (*subset_cb) (int subset_bitmap, void *data);

static int
ncr (int n, int r) 
{
	int i;
	int prod = 1;

	for (i = 1; i <= r; i++)
		prod = (prod * (i + n - r)) / i;

	return prod;
}

static void
all_subsets (int n, int k, subset_cb callback, int mask, void *data) 
{
	if (n == k) callback (mask | ((1 << n) - 1), data);
	else if (n == 0) return;
	else {
		all_subsets (n - 1, k, callback, mask, data);
		all_subsets (n - 1, k - 1, callback, mask | (1 << (n - 1)), data);
	}
}

static void
big_subset_cb (int bitmap, bibd_state_t *state) 
{
	if ((state->two_elem_subset & bitmap) == state->two_elem_subset)
		matrix_set_entry (state->bibd_matrix, state->row, state->col, state->one);

	state->col++;
}

static void
small_subset_cb (int bitmap, bibd_state_t *state) 
{
	state->col = 0;
	state->two_elem_subset = bitmap;
	all_subsets (state->v, state->k, (subset_cb) big_subset_cb, 0, state);
	state->row++;
}

matrix_t *
make_bibd_matrix (int v, int k) 
{
	bibd_state_t state;

	state.one = field.one;
	state.v = v;
	state.k = k;
	state.bibd_matrix = matrix_new (&sparse_vector_class,
					ncr (v, 2), ncr (v, k));

	state.row = 0;
	all_subsets (state.v, 2, (subset_cb) small_subset_cb, 0, &state);

	return state.bibd_matrix;
}

vector_t *
make_bibd_vector (int v, int k, int lambda) 
{
	int i, len;
	vector_t *vector;
	Element entry (lambda, 1);

	len = ncr (v, 2);
	vector = vector_new (&dense_vector_class, len);

	for (i = 0; i < len; i++)
		vector_set_entry (vector, i, entry);

	return vector;
}

/* Modify a BIBD matrix according to the given includes and excludes list;
 * assumes lists given are in sorted order */

typedef struct _refine_state_t refine_state_t;

struct _refine_state_t 
{
	unsigned int *includes_list;
	unsigned int includes_len;
	unsigned int includes_idx;

	unsigned int *excludes_list;
	unsigned int excludes_len;
	unsigned int excludes_idx;

	vector_t *out;
};

static int
refine_bibd_cb (vector_t *vector, unsigned int col, Element &element,
		refine_state_t *state) 
{
	while (state->includes_idx < state->includes_len &&
	       col > state->includes_list[state->includes_idx])
		state->includes_idx++;
	while (state->excludes_idx < state->excludes_len &&
	       col > state->excludes_list[state->excludes_idx])
		state->excludes_idx++;

	if (state->includes_idx < state->includes_len &&
	    col == state->includes_list[state->includes_idx])
	{
		state->includes_idx++;
		return 0;
	}

	if (state->excludes_idx < state->excludes_len &&
	    col == state->excludes_list[state->excludes_idx])
	{
		state->excludes_idx++;
		return 0;
	}

	vector_set_entry (state->out, col - (state->includes_idx + state->excludes_idx), element);
	return 0;
}

matrix_t *
refine_bibd_matrix (matrix_t *bibd_matrix,
		    unsigned int *includes_list, int includes_len,
		    unsigned int *excludes_list, int excludes_len) 
{
	unsigned int row, rows, cols;
	refine_state_t state;
	matrix_t *out;

	rows = matrix_get_rows (bibd_matrix);
	cols = matrix_get_cols (bibd_matrix);
	out = matrix_new (&sparse_vector_class, rows, cols - (includes_len + excludes_len));

	state.includes_list = includes_list;
	state.includes_len = includes_len;
	state.excludes_list = excludes_list;
	state.excludes_len = excludes_len;

	for (row = 0; row < rows; row++) {
		state.out = matrix_get_row_vector (out, row);
		state.includes_idx = state.excludes_idx = 0;
		vector_foreach_entry (matrix_get_row_vector (bibd_matrix, row),
				      (entry_cb_t) refine_bibd_cb, 0,
				      &state);
	}

	return out;
}

/* Modify a BIBD vector by subtracting out the given column vectors from the
 * BIBD matrix */

vector_t *
refine_bibd_vector (vector_t *bibd_vector, unsigned int *list, unsigned int len,
		    matrix_t *bibd_matrix) 
{
	vector_t *out, *tmp = NULL;
	matrix_t *bibd_t;
	unsigned int i;

	bibd_t = matrix_transpose (bibd_matrix);
	out = bibd_vector;
	vector_ref (out);

	for (i = 0; i < len; i++) {
		tmp = out;
		out = vector_row_op (tmp, matrix_get_row_vector (bibd_t, list[i]),
				     field.neg_one);
		vector_unref (tmp);
	}

	matrix_destroy (bibd_t);

	return out;
}

/* Return the subset bitmap corresponding to the given column number */

unsigned int
get_subset (unsigned int column) 
{
	return 0;
}

/* Return the column number for a given subset bitmap */

unsigned int
get_column_number (unsigned int subset) 
{
	return 0;
}

/* Fill an array that maps column numbers to subsets */

typedef struct _subset_array_state_t subset_array_state_t;

struct _subset_array_state_t 
{
	unsigned int *array;
	unsigned int i;
};

static void
subset_array_cb (unsigned int bitmap, subset_array_state_t *state)
{
	state->array[state->i++] = bitmap;
}

unsigned int *
get_subset_array (unsigned int n, unsigned int r, unsigned int *len) 
{
	subset_array_state_t state;

	*len = ncr (n, r);
	state.array = new unsigned int[*len];
	state.i = 0;
	all_subsets (n, r, (subset_cb) subset_array_cb, 0, &state);
	return state.array;
}

/* Debugging function: print out all the subset bitmaps with their associated column/row numbers */

static void
print_subset_bitmap_cb (int bitmap, int *num) 
{
	cout.width (2);
	cout << *num << " ";
	cout.setf (ios::hex, ios::basefield);
	cout << bitmap << endl;
	cout.setf (ios::dec);
	(*num)++;
}

void
print_subsets (int n, int r) 
{
	int num = 0;

	all_subsets (n, r, (subset_cb) print_subset_bitmap_cb, 0, &num);
}
