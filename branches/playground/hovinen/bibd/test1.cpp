/* test1.cpp
 *
 * Test suite
 * Written by Bradford Hovinen <hovinen@ximian.com>
 */

/* Tests in this suite:
 *
 * Black box tests:
 * - Reflexive, symmetric, transitive properties of equality
 * - Multiply (sparse, dense) vector by its inverse and test for identity
 *   property
 * - Apply the identity matrix to a vector and test for equality
 * - Find the dot product of two orthogonal (sparse, dense) vectors and see if
 *   the result is 0
 * - Take a matrix known to be singular and permute it to have a nonsingular
 *   r x r minor; compute the inverse of the minor
 * - Take a matrix and vector known to be in its range, test the least-norm
 *   solution to see if |Ax-b| = 0
 * - Find the transpose of a transpose and test for equality
 * - Create a Vandermonde matrix and invert
 * - Perform computations (e.g., multiplication, inverse, transpose) on
 *   matrices using sparse and dense representations and compare the results
 *   for equality
 *
 * Glass box tests:
 * - Multiply a matrix by a permutation matrix using mult_tranpose and
 *   mult_permutation_transpose; test if the two resulting matrices are the
 *   same
 * - Inequality: sparse vs. dense vs. compact sparse, positions of 0's
 */

#include <iostream>
#include <fstream>

#include "LinBox/gmp-rational-random.h"

#include "matrix.h"
#include "dense-vector.h"
#include "sparse-vector.h"
#include "sparse-compact-vector.h"

/* test1: Equality of sparse and dense vectors */

void test1 (void) {
	vector_t *v1, *v2, *v3, *v4;
	ifstream input;
	bool bad;

	bad = false;
	cerr << "Testing properties of equality...";
	v1 = vector_new (&dense_vector_class, 200);
	v2 = vector_new (&sparse_vector_class, 200);
	v3 = vector_new (&dense_vector_class, 200);
	v4 = vector_new (&sparse_vector_class, 200);

	input.open ("test/test-1.vector");
	vector_read (v1, input);
	input.close ();

	input.open ("test/test-1.vector");
	vector_read (v2, input);
	input.close ();

	if (!vector_equal (v1, v1)) {
		cerr << endl << " *** Test failed on vector_equal (v1, v1)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v2, v2)) {
		cerr << endl << " *** Test failed on vector_equal (v2, v2)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v1, v2)) {
		cerr << endl << " *** Test failed on vector_equal (v1, v2)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v2, v1)) {
		cerr << endl << " *** Test failed on vector_equal (v2, v1)"
		     << endl;
		bad = true;
	}

	input.open ("test/test-2.vector");
	vector_read (v3, input);
	input.close ();

	input.open ("test/test-2.vector");
	vector_read (v4, input);
	input.close ();

	if (!vector_equal (v3, v3)) {
		cerr << endl << " *** Test failed on vector_equal (v3, v3)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v4, v4)) {
		cerr << endl << " *** Test failed on vector_equal (v4, v4)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v3, v4)) {
		cerr << endl << " *** Test failed on vector_equal (v3, v4)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v4, v3)) {
		cerr << endl << " *** Test failed on vector_equal (v4, v3)"
		     << endl;
		bad = true;
	}

	if (!bad)
		cout << "good" << endl;

	vector_destroy (v1);
	vector_destroy (v2);
	vector_destroy (v3);
	vector_destroy (v4);
}

/* test2: Inequality of sparse and dense vectors */

void test2 (void) {
	vector_t *v1, *v2, *v3, *v4;
	ifstream input;
	bool bad;

	bad = false;
	cerr << "Testing inequality...";
	v1 = vector_new (&sparse_vector_class, 200);
	v2 = vector_new (&sparse_vector_class, 200);
	v3 = vector_new (&sparse_vector_class, 200);
	v4 = vector_new (&sparse_vector_class, 200);

	input.open ("test/test-3-1.vector");
	vector_read (v1, input);
	input.close ();

	input.open ("test/test-3-2.vector");
	vector_read (v2, input);
	input.close ();

	if (vector_equal (v1, v2)) {
		cerr << endl << " *** Test failed on vector_equal (v1, v2)"
		     << endl;
		bad = true;
	}

	if (vector_equal (v2, v1)) {
		cerr << endl << " *** Test failed on vector_equal (v2, v1)"
		     << endl;
		bad = true;
	}

	input.open ("test/test-4-1.vector");
	vector_read (v3, input);
	input.close ();

	input.open ("test/test-4-2.vector");
	vector_read (v4, input);
	input.close ();

	if (vector_equal (v3, v4)) {
		cerr << endl << " *** Test failed on vector_equal (v3, v4)"
		     << endl;
		bad = true;
	}

	if (vector_equal (v4, v3)) {
		cerr << endl << " *** Test failed on vector_equal (v4, v3)"
		     << endl;
		bad = true;
	}

	if (!bad)
		cout << "good" << endl;

	vector_destroy (v1);
	vector_destroy (v2);
	vector_destroy (v3);
	vector_destroy (v4);
}

/* test3: Equality and inequality with a 1, 0, -1 -vector */

void test3 (void) {
	vector_t *v1, *v2, *v3;
	ifstream input;
	bool bad;

	bad = false;
	cerr << "Testing equality against 1, 0, -1 -vector...";
	v1 = vector_new (&dense_vector_class, 200);
	v2 = vector_new (&sparse_vector_class, 200);
	v3 = vector_new (&sparse_compact_vector_class, 200);

	input.open ("test/test-5.vector");
	vector_read (v1, input);
	input.close ();

	input.open ("test/test-5.vector");
	vector_read (v2, input);
	input.close ();

	input.open ("test/test-5.vector");
	vector_read (v3, input);
	input.close ();

	if (!vector_equal (v3, v3)) {
		cerr << endl << " *** Test failed on vector_equal (v3, v3)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v1, v3)) {
		cerr << endl << " *** Test failed on vector_equal (v1, v3)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v3, v1)) {
		cerr << endl << " *** Test failed on vector_equal (v3, v1)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v2, v3)) {
		cerr << endl << " *** Test failed on vector_equal (v2, v3)"
		     << endl;
		bad = true;
	}

	if (!vector_equal (v3, v2)) {
		cerr << endl << " *** Test failed on vector_equal (v3, v2)"
		     << endl;
		bad = true;
	}

	if (!bad)
		cout << "good" << endl;

	vector_destroy (v1);
	vector_destroy (v2);
	vector_destroy (v3);
}

/* test4: Take the transpose of a transpose and test for equality */

void test4 (vector_class_t *klass, char *file) 
{
	matrix_t *m1, *mt, *mtt;
	ifstream input;
	bool bad = false;

	cout << "Testing whether (A^T)^T=A (";

	if (klass == &dense_vector_class)
		cout << "dense";
	else if (klass == &sparse_vector_class)
		cout << "sparse";
	else
		cout << "???";

	cout << ")...";

	input.open (file);
	m1 = matrix_read (klass, input);
	input.close ();

	mt = matrix_transpose (m1);

	if (mt == NULL) {
		cout << endl << " *** Could not compute transpose." << endl;
		bad = true;
	}

	mtt = matrix_transpose (mt);

	if (mtt == NULL) {
		cout << endl << " *** Could not compute transpose." << endl;
		bad = true;
	}

	if (!matrix_equal (m1, mtt)) {
		cout << endl << " *** Matrices are not equal." << endl;
		bad = true;
	}

	if (!bad)
		cout << "good" << endl;

	matrix_destroy (m1);
	matrix_destroy (mt);
	matrix_destroy (mtt);
}

/* test5: Apply the identity matrix to a vector and test for equality */

void test5 (vector_class_t *klass) 
{
	matrix_t *id;
	vector_t *v1, *v2;
	ifstream input;
	bool bad = false;

	cout << "Testing whether I*x=x (";

	if (klass == &dense_vector_class)
		cout << "dense";
	else if (klass == &sparse_vector_class)
		cout << "sparse";
	else if (klass == &sparse_compact_vector_class)
		cout << "sparse compact";
	else
		cout << "???";

	cout << ")...";

	id = matrix_new_identity (klass, 200);

	if (id == NULL) {
		cout << endl << " *** Could not create identity matrix." << endl;
		bad = true;
	}

	v1 = vector_new (klass, 200);

	input.open ("test/test-1.vector");
	vector_read (v1, input);
	input.close ();

	v2 = matrix_apply_to_vector (id, &dense_vector_class, v1);

	if (!vector_equal (v1, v2)) {
		cout << endl << " *** Vectors are not equal" << endl;
		bad = true;
	}

	if (!bad)
		cout << "good" << endl;

	matrix_destroy (id);

	vector_destroy (v1);
	vector_destroy (v2);
}

/* test6: Multiply matrix by its inverse, using different representations */

void test6 (vector_class_t *klass) 
{
	matrix_t *m1, *m2, *minv, *prod, *id;
	ifstream input;
	bool bad = false;

	cout << "Testing A*A^-1 for identity property (";

	if (klass == &dense_vector_class)
		cout << "dense";
	else if (klass == &sparse_vector_class)
		cout << "sparse";
	else
		cout << "???";

	cout << ")...";
	cout.flush ();

	input.open ("test/test-1.matrix");
	m1 = matrix_read (klass, input);
	input.close ();

	m2 = matrix_clone (m1);
	minv = matrix_inv (m2, klass);

	if (minv == NULL) {
		cout << endl << " *** Could not invert matrix; since it was a dense matrix"
		     << endl << "     this is probably an error." << endl;
		bad = true;
	}

	prod = matrix_mult (&dense_vector_class, m1, minv);

	if (prod == NULL) {
		cout << endl << " *** Could not get matrix product" << endl;
		bad = true;
	}

	id = matrix_new_identity (klass, 200);

	if (!matrix_equal (id, prod)) {
		cout << endl << " *** A*A^-1 is not the identity matrix" << endl;
		bad = true;
	}

	if (!bad)
		cout << "good" << endl;

	matrix_destroy (m1);
	matrix_destroy (minv);
	matrix_destroy (prod);
}

/* test7: Take inner product of two known orthogonal sparse vectors and check
 * if the result is 0
 */

void test7 (vector_class_t *klass) 
{
	vector_t *v1, *v2;
	bool bad = false;
	ifstream input; 
	Element prod;

	cout << "Testing inner product of orthogonal vectors (";
	
	if (klass == &dense_vector_class)
		cout << "dense";
	else if (klass == &sparse_vector_class)
		cout << "sparse";
	else
		cout << "???";

	cout << ")...";
	cout.flush ();

	v1 = vector_new (klass, 200);
	v2 = vector_new (klass, 200);

	input.open ("test/test-5-1.vector");
	vector_read (v1, input);
	input.close ();

	input.open ("test/test-5-2.vector");
	vector_read (v2, input);
	input.close ();

	vector_dot_product (v1, v2, prod);

	if (!field.isZero (prod)) {
		cout << endl << " *** Dot product is nonzero" << endl;
		bad = true;
	}

	if (!bad)
		cout << "good" << endl;

	vector_destroy (v1);
	vector_destroy (v2);
}

/* test8: Apply a 1,0,-1 -matrix to a vector using all different
 * representations and test for equality
 */

void test8 (void) 
{
	matrix_t *m[3];
	vector_t *v[2], *res[2][3];

	static const char *names[] = { "dense", "sparse", "sparse compact" };
	static vector_class_t *classes[] = {
		&dense_vector_class, &sparse_vector_class, &sparse_compact_vector_class
	};

	ifstream input;
	bool bad = false;
	int i, j;

	cout << "Testing application of 1,0,-1 -matrix to a vector...";
	cout.flush ();

	for (i = 0; i < 3; i++) {
		input.open ("test/test-3.matrix");
		m[i] = matrix_read (classes[i], input);
		input.close ();
	}

	for (i = 0; i < 2; i++) {
		v[i] = vector_new (classes[i], 200);
		input.open ("test/test-2.vector");
		vector_read (v[i], input);
		input.close ();
	}

	for (i = 0; i < 2; i++)
		for (j = 0; j < 3; j++)
			res[i][j] = matrix_apply_to_vector (m[j], &sparse_vector_class, v[i]);

	for (i = 0; i < 2; i++) {
		for (j = 0; j < 3; j++) {
			if (i == 0 && j == 0) continue;
			if (!vector_equal (res[i][j], res[0][0])) {
				cout << endl << "Vectors are not equal (matrix is " << names[j] << ", vector is " << names[i] << ")" << endl;
				bad = true;
			}
		}
	}

	if (!bad)
		cout << "good" << endl;

	for (i = 0; i < 2; i++)
		for (j = 0; j < 3; j++)
			vector_destroy (res[i][j]);

	for (i = 0; i < 2; i++)
		vector_destroy (v[i]);

	for (i = 0; i < 3; i++)
		matrix_destroy (m[i]);
}

/* test9: Test row operations */

void test9 (void) 
{
	Random random (field);
	vector_t *a, *alpha_a, *b, *c, *alpha_c, *d, *e;
	Element alpha, neg_alpha, beta;
	ifstream input;
	bool bad = false;

	cout << "Testing row operations...";

	alpha = random ();
	beta = random ();

	a = vector_new (&sparse_vector_class, 200);
	input.open ("test/test-5-1.vector");
	vector_read (a, input);
	input.close ();

	b = vector_new (&sparse_vector_class, 200);
	input.open ("test/test-5-2.vector");
	vector_read (b, input);
	input.close ();

	field.neg (neg_alpha, alpha);

	alpha_a = vector_clone (a);
	vector_scale (alpha_a, alpha);

	c = vector_row_op (alpha_a, b, beta);

	alpha_c = vector_clone (c);
	vector_scale (alpha_c, alpha);

	d = vector_row_op (c, a, alpha);
	e = vector_row_op (alpha_c, d, neg_alpha);

	vector_scale (a, alpha);
	vector_scale (a, neg_alpha);

	if (!vector_equal (e, a)) {
		cout << endl << " *** Vectors are not equal" << endl;
		bad = true;
	}

	if (!bad)
		cout << "good" << endl;

	vector_destroy (a);
	vector_destroy (b);
	vector_destroy (c);
	vector_destroy (d);
	vector_destroy (e);
	vector_destroy (alpha_a);
	vector_destroy (alpha_c);
}

int main (int argc, char **argv)
{
	test1 ();
	test2 ();
	test3 ();
	test4 (&dense_vector_class, "test/test-1.matrix");
	test4 (&dense_vector_class, "test/test-2.matrix");
	test4 (&sparse_vector_class, "test/test-1.matrix");
	test4 (&sparse_vector_class, "test/test-2.matrix");
	test5 (&dense_vector_class);
	test5 (&sparse_vector_class);
	test5 (&sparse_compact_vector_class);
	test7 (&dense_vector_class);
	test7 (&sparse_vector_class);
	test8 ();
	test9 ();
	test6 (&dense_vector_class);
	test6 (&sparse_vector_class);
}
