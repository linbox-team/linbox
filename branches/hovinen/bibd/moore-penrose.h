/* moore-penrose.h
 *
 * Header file for Moore Penrose algorithms
 *
 * Written by Bradford Hovinen <hovinen@ximian.com>
 */

#ifndef __MOORE_PENROSE_H
#define __MOORE_PENROSE_H

#include "matrix.h"
#include "vector.h"

vector_t *moore_penrose (matrix_t *A, vector_t *b, unsigned int rank);
vector_t *full_moore_penrose (matrix_t *A, vector_t *b, unsigned int rank);
vector_t *simple_moore_penrose (matrix_t *A, vector_t *b, unsigned int rank);

#endif /* __MOORE_PENROSE_H */
