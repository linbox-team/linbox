/* moore-penrose.cpp
 *
 * Implementation of an algorithm for computing the Moore-Penrose inverse of
 * a sparse matrix
 *
 * Written by Bradford Hovinen <hovinen@ximian.com>
 */

#include <assert.h>

#include "moore-penrose.h"
#include "dense-vector.h"

/* Compute the result of applying the Moore-Penrose inverse of a matrix A to
 * the vector b. Assume the leading principle rank x rank minor of the matrix
 * is nonsingular
 */

vector_t *
moore_penrose (matrix_t *A, vector_t *b, unsigned int rank) 
{
	matrix_t *AP, *ABA, *AAC, *ABAT, *AACT, *ABATP, *AACTP, *ABATINV, *AACTINV;
	vector_t *bp, *tmp;

	AP = matrix_get_block_submatrix (A, 0, 0, rank, rank);

	ABA = matrix_get_block_submatrix (A, 0, 0, rank, matrix_get_cols (A));
	ABATP = matrix_mult_transpose (&dense_vector_class, ABA, ABA);
	ABATINV = matrix_inv (ABATP, &dense_vector_class);
	assert (ABATINV != NULL);
	matrix_destroy (ABATP);
	ABAT = matrix_transpose (ABA);
	matrix_destroy (ABA);

	AAC = matrix_get_block_submatrix (A, 0, 0, matrix_get_rows (A), rank);
	AACT = matrix_transpose (AAC);
	AACTP = matrix_mult_transpose (&dense_vector_class, AACT, AACT);
	AACTINV = matrix_inv (AACTP, &dense_vector_class);
	assert (AACTINV != NULL);
	matrix_destroy (AACTP);

	bp = matrix_apply_to_vector (AACT, &dense_vector_class, b);
	matrix_destroy (AACT);
	tmp = matrix_apply_to_vector (AACTINV, &dense_vector_class, bp);
	matrix_destroy (AACTINV);
	vector_unref (bp);
	bp = tmp;

	tmp = matrix_apply_to_vector (AP, &dense_vector_class, bp);
	matrix_destroy (AP);
	vector_unref (bp);
	bp = tmp;

	tmp = matrix_apply_to_vector (ABATINV, &dense_vector_class, bp);
	matrix_destroy (ABATINV);
	vector_unref (bp);
	bp = tmp;

	tmp = matrix_apply_to_vector (ABAT, &dense_vector_class, bp);
	matrix_destroy (ABAT);
	vector_unref (bp);

	return tmp;
}

/* Compute the result of applying the Moore-Penrose inverse of a matrix A to
 * the vector b; first reorder the rows and columns of the matrix so that the
 * leading principle rank x rank minor is nonsingular
 */

vector_t *
full_moore_penrose (matrix_t *A, vector_t *b, unsigned int rank) 
{
	vector_t *result, *tmp;
	matrix_t *PAQ, *P, *Q;

	PAQ = matrix_nonsingular_minor (A, rank, &P, &Q);
	tmp = matrix_apply_to_vector (P, &dense_vector_class, b);
	matrix_destroy (P);

	result = moore_penrose (PAQ, tmp, rank);
	matrix_destroy (PAQ);
	vector_unref (tmp);
	tmp = result;

	result = matrix_apply_to_vector (Q, &dense_vector_class, tmp);
	matrix_destroy (Q);
	vector_unref (tmp);
	return result;
}

/* Compute the result of applying the Moore-Penrose inverse of a matrix A to
 * the vector b in the special case that the first rank rows of A are
 * independent.
 */

vector_t *
simple_moore_penrose (matrix_t *A, vector_t *b, unsigned int rank) 
{
	matrix_t *AP, *APT, *AAT, *AATINV;
	vector_t *bp, *tmp;

	AP = matrix_get_block_submatrix (A, 0, 0, rank, matrix_get_cols (A));
	APT = matrix_transpose (AP);
	AAT = matrix_mult_transpose (&dense_vector_class, AP, AP);
	AATINV = matrix_inv (AAT, &dense_vector_class);

	matrix_destroy (AP);

	if (AATINV == NULL) {
		rank = MIN (rank, matrix_get_cols (A));
		matrix_destroy (APT);
		matrix_destroy (AAT);
		bp = vector_clone (b);
		matrix_row_echelon (A, bp);
		AP = matrix_get_block_submatrix (A, 0, 0, rank, rank);
		tmp = vector_get_block_subvector (bp, 0, rank);
		vector_unref (bp);
		bp = tmp;
		AATINV = matrix_inv (AP, &dense_vector_class);
		matrix_destroy (AP);
		tmp = matrix_apply_to_vector (AATINV, &dense_vector_class, bp);
		matrix_destroy (AATINV);
		vector_unref (bp);
		return tmp;
	}

	matrix_destroy (AAT);

	bp = vector_get_block_subvector (b, 0, rank);
	tmp = matrix_apply_to_vector (AATINV, &dense_vector_class, bp);
	matrix_destroy (AATINV);
	vector_unref (bp);
	bp = tmp;
	tmp = matrix_apply_to_vector (APT, &dense_vector_class, bp);
	matrix_destroy (APT);
	vector_unref (bp);
	return tmp;
}
