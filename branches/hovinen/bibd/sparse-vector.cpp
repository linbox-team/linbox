/* sparse-vector.c
 *
 * Sparse vector subsystem
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 */

#include <stdlib.h>
#include <string.h>

#include "sparse-vector.h"

#define INIT_ALLOC 16

typedef struct _sparse_vector_t sparse_vector_t;
typedef struct _sparse_entry_t sparse_entry_t;

struct _sparse_entry_t
{
	Element     value;
	unsigned int            column;
};

struct _sparse_vector_t 
{
	vector_t        parent;
	unsigned int    length;
	unsigned int    used;
	unsigned int    alloc_len;
	sparse_entry_t *entries;
};

static void sparse_vector_set_entry          (vector_t *v,
					      unsigned int column,
					      const Element &value);
static unsigned sparse_vector_get_column_idx (sparse_vector_t *vector,
					      unsigned col);
static unsigned int_sparse_vector_get_column_idx (sparse_vector_t *vector,
						  unsigned col, 
						  unsigned low,
						  unsigned high);
static void sparse_vector_check_alloc        (sparse_vector_t *vector,
					      unsigned int n);

/* sparse_vector_new
 *
 * Construct a new vector structure with room for n entries
 */

static vector_t *
sparse_vector_new (unsigned int length) 
{
	sparse_vector_t *vector;
	unsigned i;

	assert (length > 0);

	vector = new sparse_vector_t;
	vector->length = length;
	vector->used = 0;
	vector->alloc_len = INIT_ALLOC;
	vector->entries = new sparse_entry_t[INIT_ALLOC];

	for (i = 0; i < vector->alloc_len; i++) {
		vector->entries[i].value = field.zero;
		vector->entries[i].column = (unsigned int) -1;
	}

	return (vector_t *) vector;
}

/* sparse_vector_clone
 *
 * Clone a sparse vector
 */

static vector_t *
sparse_vector_clone (vector_t *v) 
{
	sparse_vector_t *vector = (sparse_vector_t *) v;
	sparse_vector_t *out;
	unsigned i;

	out = (sparse_vector_t *) vector_new (&sparse_vector_class,
					      vector->length);
	out->used = vector->used;
	out->alloc_len = vector->alloc_len;
	delete[] out->entries;
	out->entries = new sparse_entry_t[out->alloc_len];

	for (i = 0; i < vector->used; i++)
		out->entries[i] = vector->entries[i];

	return (vector_t *) out;
}

static int
construct_cb (vector_t *vector, unsigned int col, Element &value,
	      sparse_vector_t *new_vector)
{
	sparse_vector_set_entry ((vector_t *) new_vector, col, value);
	return 0;
}

/* sparse_vector_construct
 *
 * Construct a new vector from an existing one that is not necessarily of the
 * same class
 */

static vector_t *
sparse_vector_construct (vector_t *vector) 
{
	sparse_vector_t *new_vector;

	new_vector = (sparse_vector_t *)
		vector_new (&sparse_vector_class, vector_length (vector));
	vector_foreach_entry (vector, (entry_cb_t) construct_cb,
			      0, new_vector);
	return (vector_t *) new_vector;
}

/* sparse_vector_destroy
 *
 * Destroy a vector structure */

static void
sparse_vector_destroy (vector_t *vector)
{
	sparse_vector_t *sparse_vector = (sparse_vector_t *) vector;

	delete[] sparse_vector->entries;
	delete sparse_vector;
}

/* sparse_vector_length
 *
 * Return the length of a sparse vector
 */

static unsigned int
sparse_vector_length (vector_t *vector)
{
	return ((sparse_vector_t *) vector)->length;
}

/* sparse_vector_num_entries
 *
 * Return the number of entries in a sparse vector
 */

static unsigned int
sparse_vector_num_entries (vector_t *vector)
{
	return ((sparse_vector_t *) vector)->used;
}

/* sparse_vector_set_entry
 *
 * Set the value for the given column, optimized for the case where data are
 * being appended to the end of the column
 */

static void
sparse_vector_set_entry (vector_t *v, unsigned int column, 
			 const Element &value)
{
	sparse_vector_t *vector = (sparse_vector_t *) v;
	int idx;

	if (vector->used == 0 ||
	    vector->entries[vector->used - 1].column < column)
	{
		if (field.isZero (value)) return;
		sparse_vector_check_alloc (vector, vector->used + 1);
		vector->entries[vector->used].column = column;
		vector->entries[vector->used].value = value;
		vector->used++;
		return;
	}

	idx = sparse_vector_get_column_idx (vector, column);

	if (vector->entries[idx].column == column) {
		vector->entries[idx].value = value;
	} else {
		if (field.isZero (value)) return;
		sparse_vector_check_alloc (vector, vector->used + 1);
		memmove (&(vector->entries[idx + 1]), &(vector->entries[idx]),
			 sizeof (sparse_entry_t) * vector->used - idx);
		vector->entries[idx].column = column;
		vector->entries[idx].value = value;
		vector->used++;
	}
}

/* sparse_vector_get_entry
 *
 * Returns the element at the given column position
 */

static Element &
sparse_vector_get_entry (vector_t *v, unsigned int column,
			 Element &res) 
{
	sparse_vector_t *vector = (sparse_vector_t *) v;
	int idx;

	idx = sparse_vector_get_column_idx (vector, column);

	if (vector->entries[idx].column == column)
		res = vector->entries[idx].value;
	else
		res = field.zero;

	return res;
}

/* sparse_vector_foreach_entry
 *
 * Issues the given callback for each nonzero entry in the vector. Terminates
 * if any callback returns a nonzero value.
 *
 * Returns the number of entries iterated
 */

static int
sparse_vector_foreach_entry (vector_t *v, entry_cb_t callback, int iterate_all,
			     void *data) 
{
	sparse_vector_t *vector = (sparse_vector_t *) v;
	Element zero = field.zero;
	unsigned i, col;

	if (iterate_all) {
		for (col = 0, i = 0; col < vector->length; col++) {
			if (vector->entries[i].column == col) {
				if (callback (v, col, vector->entries[i].value,
					      data))
					return col;
				i++;
			} else {
				if (callback (v, col, zero, data))
					return col;
			}
		}

		return vector->length;
	} else {
		for (i = 0; i < vector->used; i++)
			if (callback (v, vector->entries[i].column,
				      vector->entries[i].value, data))
				return i;

		return vector->used;
	}
}

/* sparse_vector_get_block_subvector
 *
 * Construct a new vector that contains elements a1 up to but not including a2
 * of the given vector
 */

static vector_t *
sparse_vector_get_block_subvector (vector_t *inv, unsigned int a1,
				   unsigned int a2) 
{
	sparse_vector_t *in = (sparse_vector_t *) inv;
	sparse_vector_t *out;
	unsigned idx, i;

	idx = sparse_vector_get_column_idx (in, a1);

	out = (sparse_vector_t *) vector_new (&sparse_vector_class, a2 - a1);

	for (i = 0; in->entries[i + idx].column < a2 && i + idx < in->used;
	     i++)
		sparse_vector_set_entry ((vector_t *) out,
					 in->entries[i + idx].column - a1,
					 in->entries[i + idx].value);

	return (vector_t *) out;
}

/* sparse_vector_glue
 *
 * Append one vector onto another vector, forming a new vector
 */

static vector_t *
sparse_vector_glue (vector_t *v1, vector_t *v2) 
{
	sparse_vector_t *sv1 = (sparse_vector_t *) v1;
	sparse_vector_t *sv2 = (sparse_vector_t *) v2;
	sparse_vector_t *out;
	unsigned i;

	out = (sparse_vector_t *)
		vector_new (&sparse_vector_class, sv1->length + sv2->length);
	sparse_vector_check_alloc (out, sv1->used + sv2->used);

	for (i = 0; i < sv1->used; i++) {
		out->entries[i].column = sv1->entries[i].column;
		out->entries[i].value = sv1->entries[i].value;
	}

	for (i = 0; i < sv2->used; i++) {
		out->entries[i + sv1->used].column =
			sv2->entries[i].column + sv1->length;
		out->entries[i + sv1->used].value = sv2->entries[i].value;
	}

	out->used = sv1->used + sv2->used;

	return (vector_t *) out;
}

typedef struct _dot_product_state_t dot_product_state_t;

struct _dot_product_state_t 
{
	sparse_vector_t *v1;
	Element *tmp_res;
	unsigned int c_idx;
};

static int
dot_product_cb (vector_t *v2, unsigned int col, Element &value,
		dot_product_state_t *state) 
{
	Element tmp;

	while (state->c_idx < state->v1->used && state->v1->entries[state->c_idx].column < col) 
		state->c_idx++;

	if (state->c_idx >= state->v1->used || state->v1->entries[state->c_idx].column > col) {
		return 0;
	} else {
		field.mul (tmp, value, state->v1->entries[state->c_idx].value);
		field.add (*state->tmp_res, *state->tmp_res, tmp);
		state->c_idx++;
		return 0;
	}
}

/* sparse_vector_dot_product
 *
 * Compute the dot product of two sparse vectors
 */

static Element &
sparse_vector_dot_product (vector_t *v1, vector_t *v2, Element &res) 
{
	dot_product_state_t state;

	res = field.zero;
	state.v1 = (sparse_vector_t *) v1;
	state.c_idx = 0;
	state.tmp_res = &res;
	vector_foreach_entry (v2, (entry_cb_t) dot_product_cb, 0, &state);

	return res;
}

/* sparse_vector_row_op
 *
 * Perform a row operation on the given sparse vector by scaling the second
 * vector by factor and adding it to te vector; return the newly created
 * vector
 *
 * We support sparse/sparse and sparse/-- separately here. The sparse/sparse 
 * case is O(n), while the sparse/-- case is O(n lg n)
 */

typedef struct _row_op_state_t row_op_state_t;

struct _row_op_state_t 
{
	sparse_vector_t *out;
	Element *prod;
	const Element *factor;
};

static int
row_op_cb (vector_t *v2, unsigned int col, Element &entry, row_op_state_t *state) 
{
	unsigned int idx;

	idx = sparse_vector_get_column_idx (state->out, col);
	field.mul (*state->prod, entry, *state->factor);
	field.addin (state->out->entries[idx].value, *state->prod);

	return 0;
}

static vector_t *
sparse_vector_row_op (vector_t *v1, vector_t *v2,
		      const Element &factor) 
{
	sparse_vector_t *vector = (sparse_vector_t *) v1;
	sparse_vector_t *xform;
	sparse_vector_t *out;
	unsigned int vidx = 0, xidx = 0;
	Element newv;

	assert (vector->length == vector_length (v2));

	if (v2->klass == &sparse_vector_class) {
		out = (sparse_vector_t *) vector_new (&sparse_vector_class,
						      vector->length);
		xform = (sparse_vector_t *) v2;

		sparse_vector_check_alloc (out, vector->used + xform->used);

		while (vidx < vector->used || xidx < xform->used) {
			while (vidx < vector->used && (xidx >= xform->used ||
						       vector->entries[vidx].column <
						       xform->entries[xidx].column))
			{
				out->entries[out->used].column =
					vector->entries[vidx].column;
				out->entries[out->used].value =
					vector->entries[vidx].value;
				vidx++; out->used++;
			}

			while (xidx < xform->used && (vidx == vector->used ||
						      vector->entries[vidx].column >
						      xform->entries[xidx].column))
			{
				field.mul (newv, factor,
					   xform->entries[xidx].value);

				if (!field.isZero (newv)) {
					out->entries[out->used].column =
						xform->entries[xidx].column;
					out->entries[out->used].value = newv;
					out->used++;
				}

				xidx++;
			}

			while (vidx < vector->used &&
			       xidx < xform->used &&
			       vector->entries[vidx].column ==
			       xform->entries[xidx].column) 
			{
				Element tmp;

				field.mul (tmp, factor, xform->entries[xidx].value);
				field.add (newv, vector->entries[vidx].value, tmp);

				if (!field.isZero (newv)) {
					out->entries[out->used].column =
						vector->entries[vidx].column;
					out->entries[out->used].value = newv;
					out->used++;
				}

				vidx++; xidx++;
			}
		}
	} else {
		row_op_state_t state = { NULL, NULL, &factor };

		out = (sparse_vector_t *) sparse_vector_clone (v1);
		state.out = out;
		state.prod = &newv;
		vector_foreach_entry (v2, (entry_cb_t) row_op_cb, 0, &state);
	}

	return (vector_t *) out;
}

/* sparse_vector_scale
 *
 * Scales the given vector by the given factor
 */

static void
sparse_vector_scale (vector_t *v, Element &factor) 
{
	sparse_vector_t *vector = (sparse_vector_t *) v;
	unsigned i;

	for (i = 0; i < vector->used; i++)
		field.mul (vector->entries[i].value,
			   vector->entries[i].value, factor);
}

/* sparse_vector_read
 *
 * Read a vector from a stream
 */

static void sparse_vector_read (vector_t *v, istream &file) 
{
	sparse_vector_t *vector = (sparse_vector_t *) v;
	unsigned i;
	Element elem;

	for (i = 0; i < vector->length; i++) {
		field.read (file, elem);
		sparse_vector_set_entry (v, i, elem);
	}
}

/* sparse_vector_write
 *
 * Write a vector to a stream
 */

static void sparse_vector_write (vector_t *v, ostream &file)
{
	unsigned i;
	sparse_vector_t *vector = (sparse_vector_t *) v;

	file << vector->length << " V" << endl;

	for (i = 0; i < vector->used; i++) {
		file << vector->entries[i].column << " ";
		field.write (file, vector->entries[i].value);
		file << endl;
	}
}

typedef struct _equal_state_t equal_state_t;

struct _equal_state_t 
{
	sparse_vector_t *v1;
	unsigned int c_idx;
};

static int equal_cb (vector_t *v2, unsigned int col, Element &value,
		     equal_state_t *state) 
{
	if (state->v1->entries[state->c_idx].column == col && 
	    field.areEqual (value, state->v1->entries[state->c_idx].value))
		state->c_idx++;
	else if (!field.isZero (value))
		return 1;

	return 0;
}

/* sparse_vector_equal
 *
 * Test whether two vectors are equal; return 1 if so and 0 if not
 */

static int sparse_vector_equal (vector_t *v1, vector_t *v2)
{
	equal_state_t state;

	state.v1 = (sparse_vector_t *) v1;
	state.c_idx = 0;
	vector_foreach_entry (v2, (entry_cb_t) equal_cb, 0, &state);

	if (state.c_idx == state.v1->used)
		return 1;
	else
		return 0;
}

/* int_sparse_vector_get_column_idx
 *
 * Internal function. Returns the index with the smallest column number not
 * less than col.
 */

static unsigned
int_sparse_vector_get_column_idx (sparse_vector_t *vector,
				  unsigned col, unsigned low, unsigned high) 
{
	unsigned mid;

	if (low == high) {
		if (vector->entries[low].column >= col)
			return low;
		else
			return low + 1;
	}

	mid = (high + low) / 2;

	if (vector->entries[mid].column == col) return mid;
	else if (vector->entries[mid].column > col) 
		return int_sparse_vector_get_column_idx (vector, col,
							 low, mid);
	else if (high - low == 1)
		return high;
	else
		return int_sparse_vector_get_column_idx (vector, col,
							 mid, high);
}

/* sparse_vector_get_column_idx
 *
 * Get the index in the entries array associated with the given column; or the
 * index with the smallest column number not less than the given number if no
 * such column exists.
 */

static unsigned
sparse_vector_get_column_idx (sparse_vector_t *vector, unsigned col) 
{
	return int_sparse_vector_get_column_idx (vector, col, 0, vector->used);
}

/* sparse_vector_check_alloc
 *
 * Check the allocation of the sparse vector to make sure there's room for the
 * given number of elements
 */

static void
sparse_vector_check_alloc (sparse_vector_t *vector, unsigned int n)
{
	unsigned i;
	int factor = 1;
	sparse_entry_t *array;

	while (vector->alloc_len * factor < n + 1) factor <<= 1;

	if (factor > 1) {
		vector->alloc_len *= factor;
		array = new sparse_entry_t[vector->alloc_len];

		for (i = 0; i < vector->used; i++)
			array[i] = vector->entries[i];

		delete[] vector->entries;
		vector->entries = array;

		for (i = vector->used + 1; i < vector->alloc_len; i++) {
			vector->entries[i].value = field.zero;
			vector->entries[i].column = (unsigned int) -1;
		}
	}
}

vector_class_t sparse_vector_class = {
	sparse_vector_new,
	sparse_vector_clone,
	sparse_vector_construct,
	sparse_vector_destroy,

	sparse_vector_length,
	sparse_vector_num_entries,

	sparse_vector_set_entry,
	sparse_vector_get_entry,
	sparse_vector_foreach_entry,

	sparse_vector_get_block_subvector,

	sparse_vector_glue,
	sparse_vector_dot_product,
	sparse_vector_row_op,
	sparse_vector_scale,

	sparse_vector_read,
	sparse_vector_write,

	sparse_vector_equal
};
