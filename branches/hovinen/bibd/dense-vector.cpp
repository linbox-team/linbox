/* dense-vector.c
 *
 * Dense vector subsystem
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 */

#include <stdlib.h>

#include "dense-vector.h"

typedef struct _dense_vector_t dense_vector_t;

struct _dense_vector_t 
{
	vector_t             parent;
	unsigned int         length;
	Element *entries;
};

/* dense_vector_new
 *
 * Construct a new dense_vector structure with room for n entries
 */

static vector_t *
dense_vector_new (unsigned int length) 
{
	dense_vector_t *dense_vector;
	unsigned int i;

	assert (length > 0);

	dense_vector = (dense_vector_t *) calloc (sizeof (dense_vector_t), 1);
	dense_vector->length = length;
	dense_vector->entries =
		(Element *) calloc (sizeof (Element), length);

	for (i = 0; i < dense_vector->length; i++)
		dense_vector->entries[i] = field.zero;

	return (vector_t *) dense_vector;
}

/* dense_vector_clone
 *
 * Copy a dense_vector
 */

static vector_t *
dense_vector_clone (vector_t *v) 
{
	dense_vector_t *vector = (dense_vector_t *) v;
	dense_vector_t *out;
	unsigned int i;

	out = (dense_vector_t *) vector_new (&dense_vector_class,
					     vector->length);

	for (i = 0; i < vector->length; i++)
		out->entries[i] = vector->entries[i];

	return (vector_t *) out;
}

static int
construct_cb (vector_t *vector, unsigned int col, Element &value,
	      dense_vector_t *new_vector)
{
	new_vector->entries[col] = value;
	return 0;
}

/* dense_vector_construct
 *
 * Construct a new vector from an existing one that is not necessarily of the
 * same class
 */

static vector_t *
dense_vector_construct (vector_t *vector) 
{
	dense_vector_t *new_vector;

	new_vector = (dense_vector_t *) 
		vector_new (&dense_vector_class, vector_length (vector));
	vector_foreach_entry (vector, (entry_cb_t) construct_cb,
			      0, new_vector);
	return (vector_t *) new_vector;
}

/* dense_vector_destroy
 *
 * Destroy a dense_vector structure
 */

static void
dense_vector_destroy (vector_t *v) 
{
	dense_vector_t *vector = (dense_vector_t *) v;

	free (vector->entries);
	free (vector);
}

/* dense_vector_length
 *
 * Returns the length of the vector
 */

static unsigned int
dense_vector_length (vector_t *vector) 
{
	return ((dense_vector_t *) vector)->length;
}

/* dense_vector_num_entries
 *
 * Returns the length of the vector
 */

static unsigned int
dense_vector_num_entries (vector_t *vector) 
{
	return ((dense_vector_t *) vector)->length;
}

/* dense_vector_set_entry
 *
 * Sets the value at the given column
 */

static void
dense_vector_set_entry (vector_t *v, unsigned int column,
			const Element &value) 
{
	dense_vector_t *vector = (dense_vector_t *) v;

	vector->entries[column] = value;
}

/* dense_vector_get_entry
 *
 * Returns the element at the given column position
 */

static Element &
dense_vector_get_entry (vector_t *v, unsigned int column,
			Element &res) 
{
	dense_vector_t *vector = (dense_vector_t *) v;

	res = vector->entries[column];
	return res;
}

/* dense_vector_foreach_entry
 *
 * Issues the given callback for each entry; terminates if the entry returns
 * nonzero
 *
 * Returns number of entries iterated
 */

static int
dense_vector_foreach_entry (vector_t *v, entry_cb_t callback, int iterate_all,
			    void *data) 
{
	dense_vector_t *vector = (dense_vector_t *) v;
	unsigned i;

	for (i = 0; i < vector->length; i++)
		if (callback (v, i, vector->entries[i], data))
			return i;

	return vector->length;
}

/* dense_vector_get_block_subvector
 *
 * Construct a new dense_vector that contains elements a1 through a2 of the
 * given vector
 */

static vector_t *
dense_vector_get_block_subvector (vector_t *v, unsigned int a1,
				  unsigned int a2) 
{
	dense_vector_t *vector = (dense_vector_t *) v;
	dense_vector_t *out;
	unsigned int i;

	out = (dense_vector_t *) vector_new (&dense_vector_class, a2 - a1);

	for (i = 0; i < a2 - a1; i++)
		out->entries[i] = vector->entries[i + a1];

	return (vector_t *) out;
}

/* dense_vector_glue
 *
 * Paste two dense vectors together
 */

static vector_t *
dense_vector_glue (vector_t *v1, vector_t *v2) 
{
	dense_vector_t *dv1 = (dense_vector_t *) v1;
	dense_vector_t *dv2 = (dense_vector_t *) v2;
	dense_vector_t *out;

	out = (dense_vector_t *) vector_new (&dense_vector_class,
					     dv1->length + dv2->length);

	memcpy (out->entries, dv1->entries, dv1->length * sizeof (Element));
	memcpy (&(out->entries[dv1->length]),
		dv2->entries, dv2->length * sizeof (Element));

	return (vector_t *) out;
}

typedef struct _dot_product_state_t dot_product_state_t;

struct _dot_product_state_t 
{
	Element *tmp_res;
	dense_vector_t *v1;
};

static int
dot_product_cb (vector_t *v2, unsigned int col, Element &value,
		dot_product_state_t *state) 
{
	Element tmp;

	field.mul (tmp, value, state->v1->entries[col]);
	field.add (*state->tmp_res, *state->tmp_res, tmp);
	return 0;
}

/* dense_vector_dot_product
 *
 * Find the dot product of two dense vectors
 */

static Element &
dense_vector_dot_product (vector_t *v1, vector_t *v2,
			  Element &res) 
{
	dot_product_state_t state;

	res = field.zero;
	state.v1 = (dense_vector_t *) v1;
	state.tmp_res = &res;
	vector_foreach_entry (v2, (entry_cb_t) dot_product_cb, 0, &state);

	return res;
}

/* dense_vector_row_op
 *
 * Perform a row operation on the given dense vector by scaling the second
 * vector by factor and adding it to te vector; return the newly created
 * vector. May affect first vector.
 *
 * We support two flavors here: dense<->dense (which is very slightly faster)
 * and dense<->sparse/sparse-compact. The latter uses the foreach mechanism
 */

typedef struct _row_op_state_t row_op_state_t;

struct _row_op_state_t 
{
	dense_vector_t *vector;
	Element *prod;
	const Element *factor;
};

static int
row_op_cb (vector_t *v2, unsigned int col, Element &entry, row_op_state_t *state) 
{
	field.mul (*state->prod, entry, *state->factor);
	field.add (state->vector->entries[col], state->vector->entries[col],
		   *state->prod);
	return 0;
}

static vector_t *
dense_vector_row_op (vector_t *v1, vector_t *v2,
		     const Element &factor) 
{
	dense_vector_t *vector = (dense_vector_t *) v1;
	dense_vector_t *xform;
	Element prod;
	unsigned i;
	row_op_state_t state = { NULL, NULL, &factor };

	assert (vector->length == vector_length (v2));

	if (v2->klass == &dense_vector_class) {
		xform = (dense_vector_t *) v2;

		for (i = 0; i < vector->length; i++) {
			field.mul (prod, xform->entries[i], factor);
			field.add (vector->entries[i], vector->entries[i], prod);
		}
	} else {
		state.vector = vector;
		state.factor = &factor;
		state.prod = &prod;
		vector_foreach_entry (v2, (entry_cb_t) row_op_cb, 1, &state);
	}

	vector_ref (v1);

	return v1;
}

/* dense_vector_scale
 *
 * Scales the given vector by the given factor
 */

static void
dense_vector_scale (vector_t *v, Element &factor) 
{
	dense_vector_t *vector = (dense_vector_t *) v;
	unsigned i;

	for (i = 0; i < vector->length; i++)
		field.mul (vector->entries[i], vector->entries[i], factor);
}

/* dense_vector_read
 *
 * Read a dense_vector from a stream
 */

static void
dense_vector_read (vector_t *v, istream &input) 
{
	dense_vector_t *vector = (dense_vector_t *) v;
	unsigned i;

	for (i = 0; i < vector->length; i++)
		field.read (input, vector->entries[i]);
}

/* dense_vector_write
 *
 * Write a dense_vector to a stream
 */

static void
dense_vector_write (vector_t *v, ostream &output)
{
	dense_vector_t *vector = (dense_vector_t *) v;
	unsigned i;

	for (i = 0; i < vector->length; i++) {
		field.write (output, vector->entries[i]);
		output << endl;
	}
}

/* dense_vector_equal
 *
 * Test whether two vectors are equal
 */

static int dense_vector_equal (vector_t *v1, vector_t *v2)
{
	unsigned int idx;
	dense_vector_t *dv1 = (dense_vector_t *) v1;
	dense_vector_t *dv2 = (dense_vector_t *) v2;

	if (v2->klass != &dense_vector_class)
		return vector_equal (v2, v1);

	if (dv1->length != dv2->length)
		return 0;

	for (idx = 0; idx < dv1->length; idx++) {
		if (!field.areEqual (dv1->entries[idx], dv2->entries[idx]))
			return 0;
	}

	return 1;
}

vector_class_t dense_vector_class = {
	dense_vector_new,
	dense_vector_clone,
	dense_vector_construct,
	dense_vector_destroy,

	dense_vector_length,
	dense_vector_num_entries,

	dense_vector_set_entry,
	dense_vector_get_entry,
	dense_vector_foreach_entry,

	dense_vector_get_block_subvector,

	dense_vector_glue,
	dense_vector_dot_product,
	dense_vector_row_op,
	dense_vector_scale,

	dense_vector_read,
	dense_vector_write,

	dense_vector_equal
};
