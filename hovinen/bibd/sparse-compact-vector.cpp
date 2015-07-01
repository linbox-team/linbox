/* sparse-compact-vector.c
 *
 * "Valueless" sparse vector subsystem
 *
 * Written by Bradford Hovinen <hovinen@ximian.com>
 */

#include <stdlib.h>
#include <string.h>

#include "sparse-compact-vector.h"

#define INIT_ALLOC 16

typedef struct _sparse_compact_vector_t sparse_compact_vector_t;

#define ENTRY_COLUMN(e) ((unsigned int) ((e < 0) ? (-1 - e) : (e - 1)))
#define ENTRY_VALUE(e) ((e > 0) ? 1 : ((e < 0) ? -1 : 0))

struct _sparse_compact_vector_t 
{
	vector_t        parent;
	unsigned int    length;
	unsigned int    used;
	unsigned int    alloc_len;
	int            *entries;
};

static void sparse_compact_vector_set_entry (vector_t *v,
					     unsigned int column,
					     const Element &value);
static unsigned sparse_compact_vector_get_column_idx (sparse_compact_vector_t *vector,
						      unsigned col);
static unsigned int_sparse_compact_vector_get_column_idx (sparse_compact_vector_t *vector,
							  unsigned col, 
							  unsigned low,
							  unsigned high);
static void sparse_compact_vector_check_alloc (sparse_compact_vector_t *vector,
					       unsigned int n);

/* sparse_compact_vector_new
 *
 * Construct a new vector structure with room for n entries
 */

static vector_t *
sparse_compact_vector_new (unsigned int length) 
{
	sparse_compact_vector_t *vector;
	unsigned i;

	assert (length > 0);

	vector = new sparse_compact_vector_t;
	vector->length = length;
	vector->used = 0;
	vector->alloc_len = INIT_ALLOC;
	vector->entries = new int[INIT_ALLOC];

	for (i = 0; i < vector->alloc_len; i++)
		vector->entries[i] = 0;

	return (vector_t *) vector;
}

/* sparse_compact_vector_clone
 *
 * Clone a sparse vector
 */

static vector_t *
sparse_compact_vector_clone (vector_t *v) 
{
	sparse_compact_vector_t *vector = (sparse_compact_vector_t *) v;
	sparse_compact_vector_t *out;
	unsigned i;

	out = (sparse_compact_vector_t *) vector_new (&sparse_compact_vector_class,
						      vector->length);
	out->used = vector->used;
	out->alloc_len = vector->alloc_len;
	delete[] out->entries;
	out->entries = new int[out->alloc_len];

	for (i = 0; i < vector->used; i++)
		out->entries[i] = vector->entries[i];

	return (vector_t *) out;
}

static int
construct_cb (vector_t *vector, unsigned int col, Element &value,
	      sparse_compact_vector_t *new_vector)
{
	sparse_compact_vector_set_entry ((vector_t *) new_vector, col, value);
	return 0;
}

/* sparse_compact_vector_construct
 *
 * Construct a new vector from an existing one that is not necessarily of the
 * same class
 */

static vector_t *
sparse_compact_vector_construct (vector_t *vector) 
{
	sparse_compact_vector_t *new_vector;

	new_vector = (sparse_compact_vector_t *)
		vector_new (&sparse_compact_vector_class, vector_length (vector));
	vector_foreach_entry (vector, (entry_cb_t) construct_cb,
			      0, new_vector);
	return (vector_t *) new_vector;
}

/* sparse_compact_vector_destroy
 *
 * Destroy a vector structure */

static void
sparse_compact_vector_destroy (vector_t *vector)
{
	sparse_compact_vector_t *sparse_compact_vector = (sparse_compact_vector_t *) vector;

	delete[] sparse_compact_vector->entries;
	delete sparse_compact_vector;
}

/* sparse_compact_vector_length
 *
 * Return the length of a sparse vector
 */

static unsigned int
sparse_compact_vector_length (vector_t *vector)
{
	return ((sparse_compact_vector_t *) vector)->length;
}

/* sparse_compact_vector_num_entries
 *
 * Return the number of entries in a sparse vector
 */

static unsigned int
sparse_compact_vector_num_entries (vector_t *vector)
{
	return ((sparse_compact_vector_t *) vector)->used;
}

/* sparse_compact_vector_set_entry
 *
 * Set the value for the given column, optimized for the case where data are
 * being appended to the end of the column
 */

static void
sparse_compact_vector_set_entry (vector_t *v, unsigned int column, 
			 const Element &value)
{
	sparse_compact_vector_t *vector = (sparse_compact_vector_t *) v;
	int idx, svalue;

	if (field.isZero (value))
		svalue = 0;
	else if (field.isOne (value))
		svalue = 1;
	else if (field.areEqual (value, field.neg_one))
		svalue = -1;
	else
		return;

	if (vector->used == 0 ||
	    ENTRY_COLUMN (vector->entries[vector->used - 1]) < column)
	{
		if (svalue == 0) return;
		sparse_compact_vector_check_alloc (vector, vector->used + 1);
		vector->entries[vector->used] = (column + 1) * svalue;
		vector->used++;
		return;
	}

	idx = sparse_compact_vector_get_column_idx (vector, column);

	if (ENTRY_COLUMN (vector->entries[idx]) == column && svalue != 0) {
		vector->entries[idx] =
			(ENTRY_COLUMN (vector->entries[idx]) + 1) * svalue;
	}
	else if (ENTRY_COLUMN (vector->entries[idx]) == column) {
		memmove (&(vector->entries[idx]),
			 &(vector->entries[idx + 1]),
			 sizeof (unsigned int) * vector->used - idx - 1);
		vector->used--;
	} else {
		if (svalue == 0) return;
		sparse_compact_vector_check_alloc (vector, vector->used + 1);
		memmove (&(vector->entries[idx + 1]), &(vector->entries[idx]),
			 sizeof (unsigned int) * vector->used - idx);
		vector->entries[idx] = (column + 1) * svalue;
		vector->used++;
	}
}

/* sparse_compact_vector_get_entry
 *
 * Returns the element at the given column position
 */

static Element &
sparse_compact_vector_get_entry (vector_t *v, unsigned int column,
				 Element &res) 
{
	sparse_compact_vector_t *vector = (sparse_compact_vector_t *) v;
	int idx;

	idx = sparse_compact_vector_get_column_idx (vector, column);

	if (ENTRY_COLUMN (vector->entries[idx]) != column)
		res = field.zero;
	else if (ENTRY_VALUE (vector->entries[idx]) == 1)
		res = field.one;
	else if (ENTRY_VALUE (vector->entries[idx]) == -1)
		res = field.neg_one;

	return res;
}

/* sparse_compact_vector_foreach_entry
 *
 * Issues the given callback for each nonzero entry in the vector. Terminates
 * if any callback returns a nonzero value.
 *
 * Returns the number of entries iterated
 */

static int
sparse_compact_vector_foreach_entry (vector_t *v, entry_cb_t callback, int iterate_all,
				     void *data) 
{
	sparse_compact_vector_t *vector = (sparse_compact_vector_t *) v;
	Element zero = field.zero;
	Element one = field.one;
	Element neg_one = field.neg_one;
	unsigned i, col;

	if (iterate_all) {
		for (col = 0, i = 0; col < vector->length; col++) {
			if (ENTRY_COLUMN (vector->entries[i]) != col) {
				if (callback (v, col, zero, data))
					return col;
			} else if (ENTRY_VALUE (vector->entries[i]) == 1) {
				if (callback (v, col, one, data))
					return col;
				i++;
			} else {
				if (callback (v, col, neg_one, data))
					return col;
				i++;
			}
		}

		return vector->length;
	} else {
		for (i = 0; i < vector->used; i++) {
			if (ENTRY_VALUE (vector->entries[i]) == 1) {
				if (callback (v, ENTRY_COLUMN (vector->entries[i]),
					      one, data))
					return i;
			}
			else if (ENTRY_VALUE (vector->entries[i]) == -1) {
				if (callback (v, ENTRY_COLUMN (vector->entries[i]),
					      neg_one, data))
					return i;
			}
		}

		return vector->used;
	}
}

/* sparse_compact_vector_get_block_subvector
 *
 * Construct a new vector that contains elements a1 up to but not including a2
 * of the given vector
 */

static vector_t *
sparse_compact_vector_get_block_subvector (vector_t *inv, unsigned int a1,
				   unsigned int a2) 
{
	sparse_compact_vector_t *in = (sparse_compact_vector_t *) inv;
	sparse_compact_vector_t *out;
	unsigned idx, end, i;

	idx = sparse_compact_vector_get_column_idx (in, a1);
	end = sparse_compact_vector_get_column_idx (in, a2);

	out = (sparse_compact_vector_t *) vector_new (&sparse_compact_vector_class, a2 - a1);

	sparse_compact_vector_check_alloc (out, end - idx);

	for (i = 0; i < end; i++)
		out->entries[out->used++] =
			(ENTRY_COLUMN (in->entries[i + idx]) - a1 + 1) *
			ENTRY_VALUE (in->entries[i + idx]);

	return (vector_t *) out;
}

/* sparse_compact_vector_glue
 *
 * Append one vector onto another vector, forming a new vector
 */

static vector_t *
sparse_compact_vector_glue (vector_t *v1, vector_t *v2) 
{
	sparse_compact_vector_t *sv1 = (sparse_compact_vector_t *) v1;
	sparse_compact_vector_t *sv2 = (sparse_compact_vector_t *) v2;
	sparse_compact_vector_t *out;
	unsigned i;

	out = (sparse_compact_vector_t *)
		vector_new (&sparse_compact_vector_class, sv1->length + sv2->length);
	sparse_compact_vector_check_alloc (out, sv1->used + sv2->used);

	for (i = 0; i < sv1->used; i++) {
		out->entries[i] = sv1->entries[i];
	}

	for (i = 0; i < sv2->used; i++) {
		out->entries[i + sv1->used] =
			(ENTRY_COLUMN (sv2->entries[i]) + sv1->length + 1) *
			ENTRY_VALUE (sv2->entries[i]);
	}

	out->used = sv1->used + sv2->used;

	return (vector_t *) out;
}

typedef struct _dot_product_state_t dot_product_state_t;

struct _dot_product_state_t 
{
	sparse_compact_vector_t *v1;
	unsigned int c_idx;
	Element *tmp_res;
};

static int
dot_product_cb (vector_t *v2, unsigned int col, Element &value,
		dot_product_state_t *state) 
{
	while (state->c_idx < state->v1->used &&
	       ENTRY_COLUMN (state->v1->entries[state->c_idx]) < col) 
		state->c_idx++;

	if (state->c_idx >= state->v1->used || ENTRY_COLUMN (state->v1->entries[state->c_idx]) > col) {
		return 0;
	} else {
		if (ENTRY_VALUE (state->v1->entries[state->c_idx]) == 1)
			field.add (*state->tmp_res, *state->tmp_res, value);
		else
			field.sub (*state->tmp_res, *state->tmp_res, value);
		state->c_idx++;
		return 0;
	}
}

/* sparse_compact_vector_dot_product
 *
 * Compute the dot product of two sparse vectors
 */

static Element &
sparse_compact_vector_dot_product (vector_t *v1, vector_t *v2,
				   Element &res) 
{
	dot_product_state_t state;

	res = field.zero;
	state.v1 = (sparse_compact_vector_t *) v1;
	state.c_idx = 0;
	state.tmp_res = &res;
	vector_foreach_entry (v2, (entry_cb_t) dot_product_cb, 0, &state);

	return res;
}

/* sparse_compact_vector_row_op
 *
 * Perform a row operation on the given sparse vector by scaling the second
 * vector by factor and adding it to te vector; return the newly created
 * vector
 */

static vector_t *
sparse_compact_vector_row_op (vector_t *v1, vector_t *v2,
			      const Element &factor) 
{
	cerr << "Cannot do a row operation when one of the vectors is of the"
	     << endl << "sparse compact representation." << endl;

	return (vector_t *) NULL;
}

/* sparse_compact_vector_scale
 *
 * Scales the given vector by the given factor
 */

static void
sparse_compact_vector_scale (vector_t *v, Element &factor) 
{
	cerr << "Cannot scale sparse compact vectors" << endl;
}

/* sparse_compact_vector_read
 *
 * Read a vector from a stream
 */

static void sparse_compact_vector_read (vector_t *v, istream &file) 
{
	sparse_compact_vector_t *vector = (sparse_compact_vector_t *) v;
	unsigned i;
	Element value;

	for (i = 0; i < vector->length; i++) {
		field.read (file, value);
		sparse_compact_vector_set_entry (v, i, value);
	}
}

/* sparse_compact_vector_write
 *
 * Write a vector to a stream
 */

static void sparse_compact_vector_write (vector_t *v, ostream &file)
{
	unsigned i;
	sparse_compact_vector_t *vector = (sparse_compact_vector_t *) v;

	file << vector->length << " V" << endl;

	for (i = 0; i < vector->used; i++) {
		file << ENTRY_COLUMN (vector->entries[i]) << " "
		     << ENTRY_VALUE (vector->entries[i])
		     << endl;
	}
}

typedef struct _equal_state_t equal_state_t;

struct _equal_state_t 
{
	sparse_compact_vector_t *v1;
	unsigned int c_idx;
};

static int equal_cb (vector_t *v2, unsigned int col,
		     Element &value,
		     equal_state_t *state) 
{
	if (ENTRY_COLUMN (state->v1->entries[state->c_idx]) == col) {
		if ((ENTRY_VALUE (state->v1->entries[state->c_idx]) == 1 &&
		     !field.isOne (value)) ||
		    (ENTRY_VALUE (state->v1->entries[state->c_idx]) == -1 &&
		     !field.areEqual (value, field.neg_one)))
			return 1;
		else
			state->c_idx++;
	}
	else if (!field.isZero (value))
		return 1;

	return 0;
}

/* sparse_compact_vector_equal
 *
 * Test whether two vectors are equal; return 1 if so and 0 if not
 */

static int sparse_compact_vector_equal (vector_t *v1, vector_t *v2)
{
	equal_state_t state;

	state.v1 = (sparse_compact_vector_t *) v1;
	state.c_idx = 0;
	vector_foreach_entry (v2, (entry_cb_t) equal_cb, 0, &state);

	if (state.c_idx == state.v1->used)
		return 1;
	else
		return 0;
}

/* int_sparse_compact_vector_get_column_idx
 *
 * Internal function. Returns the index with the smallest column number not
 * less than col.
 */

static unsigned
int_sparse_compact_vector_get_column_idx (sparse_compact_vector_t *vector,
					  unsigned col, unsigned low, unsigned high) 
{
	unsigned mid;

	if (low == high) {
		if (ENTRY_COLUMN (vector->entries[low]) >= col)
			return low;
		else
			return low + 1;
	}

	mid = (high + low) / 2;

	if (ENTRY_COLUMN (vector->entries[mid]) == col) return mid;
	else if (ENTRY_COLUMN (vector->entries[mid]) > col) 
		return int_sparse_compact_vector_get_column_idx (vector, col, low, mid);
	else if (high - low == 1)
		return high;
	else
		return int_sparse_compact_vector_get_column_idx (vector, col, mid, high);
}

/* sparse_compact_vector_get_column_idx
 *
 * Get the index in the entries array associated with the given column; or the
 * index with the smallest column number not less than the given number if no
 * such column exists.
 */

static unsigned
sparse_compact_vector_get_column_idx (sparse_compact_vector_t *vector, unsigned col) 
{
	return int_sparse_compact_vector_get_column_idx (vector, col, 0, vector->used);
}

/* sparse_compact_vector_check_alloc
 *
 * Check the allocation of the sparse vector to make sure there's room for the
 * given number of elements
 */

static void
sparse_compact_vector_check_alloc (sparse_compact_vector_t *vector, unsigned int n)
{
	unsigned i;
	int factor = 1;
	int *array;

	while (vector->alloc_len * factor < n + 1) factor <<= 1;

	if (factor > 1) {
		vector->alloc_len *= factor;
		array = new int[vector->alloc_len];

		for (i = 0; i < vector->used; i++)
			array[i] = vector->entries[i];

		delete[] vector->entries;
		vector->entries = array;

		for (i = vector->used + 1; i < vector->alloc_len; i++)
			vector->entries[i] = 0;
	}
}

vector_class_t sparse_compact_vector_class = {
	sparse_compact_vector_new,
	sparse_compact_vector_clone,
	sparse_compact_vector_construct,
	sparse_compact_vector_destroy,

	sparse_compact_vector_length,
	sparse_compact_vector_num_entries,

	sparse_compact_vector_set_entry,
	sparse_compact_vector_get_entry,
	sparse_compact_vector_foreach_entry,

	sparse_compact_vector_get_block_subvector,

	sparse_compact_vector_glue,
	sparse_compact_vector_dot_product,
	sparse_compact_vector_row_op,
	sparse_compact_vector_scale,

	sparse_compact_vector_read,
	sparse_compact_vector_write,

	sparse_compact_vector_equal
};
