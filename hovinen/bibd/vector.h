/* vector.h
 *
 * Prototypes for generic vector operations
 *
 * Written by Bradford Hovinen <hovinen@helixcode.com>
 */

#ifndef _VECTOR_H
#define _VECTOR_H

#include <assert.h>
#include <iostream>

#include "domain.h"

#undef __inline__

#ifdef __GNUC__
#  define __inline__ static inline
#else
#  define __inline__ static
#endif

typedef struct _vector_t vector_t;
typedef struct _vector_class_t vector_class_t;

typedef int (*entry_cb_t) (vector_t *, unsigned int,
			   Element &, void *);

struct _vector_t
{
	vector_class_t *klass;
	unsigned int    refcount;
};

struct _vector_class_t 
{
	vector_t     *(*v_new)                 (unsigned int);
	vector_t     *(*v_clone)               (vector_t *);
	vector_t     *(*v_construct)           (vector_t *);
	void          (*v_destroy)             (vector_t *);

	unsigned int  (*v_length)              (vector_t *);
	unsigned int  (*v_num_entries)         (vector_t *);

	void          (*v_set_entry)           (vector_t *, unsigned int,
						const Element &);
	Element &(*v_get_entry)    (vector_t *, unsigned int,
						Element &);

	int           (*v_foreach_entry)       (vector_t *, entry_cb_t,
						int, void *);

	vector_t     *(*v_get_block_subvector) (vector_t *, 
						unsigned int, unsigned int);

	vector_t     *(*v_glue)                (vector_t *, vector_t *);
	Element &(*v_dot_product)  (vector_t *, vector_t *,
						Element &);
	vector_t     *(*v_row_op)              (vector_t *, vector_t *, 
						const Element &);
	void          (*v_scale)               (vector_t *,
						Element &);

	void          (*v_read)                (vector_t *, istream &);
	void          (*v_write)               (vector_t *, ostream &);

	int           (*v_equal)               (vector_t *, vector_t *);
};

__inline__ vector_t *
vector_new (vector_class_t *klass, unsigned int n) 
{
	vector_t *vector;

	vector = klass->v_new (n);
	vector->klass = klass;
	vector->refcount = 1;
	return vector;
}

__inline__ vector_t *
vector_clone (vector_t *vector) 
{
	return vector->klass->v_clone (vector);
}

__inline__ vector_t *
vector_construct (vector_class_t *klass, vector_t *vector) 
{
	vector_t *new_vector;

	new_vector = klass->v_construct (vector);
	new_vector->refcount = 1;
	return new_vector;
}

__inline__ void
vector_destroy (vector_t *vector) 
{
	vector->klass->v_destroy (vector);
}

__inline__ void
vector_ref (vector_t *vector) 
{
	vector->refcount++;
}

__inline__ void
vector_unref (vector_t *vector) 
{
	vector->refcount--;

	if (vector->refcount == 0)
		vector_destroy (vector);
}

__inline__ unsigned int
vector_length (vector_t *vector) 
{
	return vector->klass->v_length (vector);
}

__inline__ unsigned int
vector_num_entries (vector_t *vector) 
{
	return vector->klass->v_num_entries (vector);
}

__inline__ void
vector_set_entry (vector_t *vector, unsigned int idx,
		  const Element &element) 
{
	vector->klass->v_set_entry (vector, idx, element);
}

__inline__ Element &
vector_get_entry (vector_t *vector, unsigned int idx,
		  Element &res) 
{
	return vector->klass->v_get_entry (vector, idx, res);
}

__inline__ int
vector_foreach_entry (vector_t *vector, entry_cb_t callback, int iterate_all,
		      void *data) 
{
	return vector->klass->v_foreach_entry (vector, callback,
					       iterate_all, data);
}

__inline__ vector_t *
vector_get_block_subvector (vector_t *vector, unsigned int start,
			    unsigned int end) 
{
	return vector->klass->v_get_block_subvector (vector, start, end);
}

__inline__ vector_t *
vector_glue (vector_t *v1, vector_t *v2) 
{
	assert (v1->klass == v2->klass);

	return v1->klass->v_glue (v1, v2);
}

__inline__ Element &
vector_dot_product (vector_t *v1, vector_t *v2, Element &res) 
{
	return v1->klass->v_dot_product (v1, v2, res);
}

__inline__ vector_t *
vector_row_op (vector_t *v1, vector_t *v2, const Element &factor) 
{
	return v1->klass->v_row_op (v1, v2, factor);
}

__inline__ void
vector_scale (vector_t *vector, Element &factor) 
{
	vector->klass->v_scale (vector, factor);
}

__inline__ void
vector_read (vector_t *vector, istream &file) 
{
	vector->klass->v_read (vector, file);
}

__inline__ void
vector_write (vector_t *vector, ostream &file) 
{
	vector->klass->v_write (vector, file);
}

__inline__ int
vector_equal (vector_t *v1, vector_t *v2)
{
	return v1->klass->v_equal (v1, v2);
}


#endif /* _VECTOR_H */
