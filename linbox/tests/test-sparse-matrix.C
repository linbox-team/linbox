/* -*- mode: c; style: linux -*- */

/* tests/test-sparse-matrix.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/param-modular.h"

#include "linbox/blackbox/sparse-matrix.h"
#include "linbox/solutions/minpoly.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Application of identity matrix onto random vectors
 *
 * Construct the identity matrix and a series of randomly-generated
 * vectors. Apply the identity to each vector and test whether the input and
 * output are equal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityApply (Field &F, size_t n, ostream &report, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	cout << "Testing identity apply...";
	cout.flush ();
	report << "Testing identity apply:" << endl;

	bool ret = true;
	bool iter_passed = true;
	Blackbox A (F, n, n);

	int i, j;
	typename Field::Element e;
	F.init (e, 1);

	for (i = 0; i < n; i++)
		A.put_value (pair<size_t, size_t> (i, i), e);

	Vector v(n), w(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;
		iter_passed = true;

		for (j = 0; j < n; j++)
			r.random (v[j]);

		report << "    Input vector:  ";
		printVector<Field> (F, report, v);

		A.apply (w, v);

		report << "    Output vector: ";
		printVector<Field> (F, report, w);

		for (j = 0; j < n; j++)
			if (!F.areEqual (w[j], v[j]))
				ret = iter_passed = false;

		if (!iter_passed)
			report << "    ERROR: Vectors are not equal" << endl;
	}

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();

	return ret;
}

/* Test 2: Application of nilpotent linear map to random vectors
 *
 * Generates an index-n nilpotent linear map and applies it to randomly
 * generated vectors n times. Generates the vectors so that they are guaranteed
 * to be outside the matrix's index- n-1 subspace; tests to make sure the
 * vectors do not prematurely reach zero.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testNilpotentApply (Field &F, size_t n, ostream &report, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	cout << "Testing nilpotent apply...";
	cout.flush ();
	report << "Testing nilpotent apply:" << endl;

	bool ret = true;
	bool iter_passed;
	bool all_zero;
	Blackbox A (F, n, n);

	int i, j;
	typename Field::Element e;
	F.init (e, 1);
	bool even = false;

	for (i = 1; i < n; i++)
		A.put_value (pair<size_t, size_t> (i - 1, i), e);

	Vector v(n), w(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;
		iter_passed = true;

		for (j = 0; j < n; j++)
			do r.random (v[j]); while (F.isZero (v[j]));

		report << "    Input vector:  ";
		printVector<Field> (F, report, v);

		for (j = 0; j < n - 1; j++, even = !even)
			if (even)
				A.apply (v, w);
			else
				A.apply (w, v);

		report << "    A^(n-1) v:     ";
		printVector<Field> (F, report, even ? w : v);

		for (all_zero = true, j = 0; j < n; j++)
			if (!F.isZero (even ? w[j] : v[j]))
				all_zero = false;

		if (all_zero) {
			ret = false;
			report << "    ERROR: A^(n-1) v is prematurely zero" << endl;
		}

		if (even)
			A.apply (v, w);
		else
			A.apply (w, v);

		report << "    A^n v:         ";
		printVector<Field> (F, report, even ? v : w);

		for (j = 0; j < n; j++)
			if (!F.isZero (even ? v[j] : w[j]))
				ret = iter_passed = false;

		if (!iter_passed)
			report << "    ERROR: A^n v is non-zero" << endl;
	}

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();

	return ret;
}

/* Test 3: Random apply to sparse matrix of K nonzero elements per row
 *
 * Generates a random sparse matrix with K nonzero elements per row and applies
 * it to the vectors {e_i | i=1..n} to test whether the output matches the ith
 * column of the matrix.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 * K - Number of nonzero elements per row
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testRandomApply1 (Field &F, size_t n, ostream &report, int iterations, int K) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	cout << "Testing sparse random apply (1)...";
	cout.flush ();
	report << "Testing sparse random apply (1):" << endl;

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	typename Field::RandIter r (F);
	typename Field::Element x;

	Integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (K > n) K = n;

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;
		iter_passed = true;

		Blackbox A (F, n, n);

		for (j = 0; j < n; j++) {
			for (k = 0; k < K; k++) {
				pair<size_t, size_t> p (j, 0);

				do
					p.second = rand () % n;
				while (!F.isZero (A[p]));

				r.random (x);
				A.put_value (p, x);
			}
		}

		report << "    Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		for (j = 0; j < n; j++) {
			Vector e_j(n), w(n);

			F.init (e_j[j], 1);

			A.apply (w, e_j);

			for (k = 0; k < n; k++)
				if (!F.areEqual (A[pair<size_t, size_t>(k, j)], w[k]))
					ret = iter_passed = false;

			report << "    Output vector " << j << ": ";
			printVector<Field> (F, report, w);
		}

		if (!iter_passed)
			report << "    ERROR: Output vectors were incorrect" << endl;
	}

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();

	return ret;
}

/* Test 4: Random apply to sparse matrix of N nonzero elements
 *
 * Generates a random sparse matrix with N nonzero elements and applies it to
 * the vectors {e_i | i=1..n} to test whether the output matches the ith column
 * of the matrix.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 * N - Number of nonzero elements
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testRandomApply2 (Field &F, size_t n, ostream &report, int iterations, int N) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	cout << "Testing sparse random apply (2)...";
	cout.flush ();
	report << "Testing sparse random apply (2):" << endl;

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	typename Field::RandIter r (F);
	typename Field::Element x;

	Integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (N > n * n) N = n * n;

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;
		iter_passed = true;

		Blackbox A (F, n, n);

		for (k = 0; k < N; k++) {
			pair<size_t, size_t> p (j, 0);

			do {
				p.first = rand () % n;
				p.second = rand () % n;
			} while (!F.isZero (A[p]));

			r.random (x);
			A.put_value (p, x);
		}

		report << "    Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		for (j = 0; j < n; j++) {
			Vector e_j(n), w(n);

			F.init (e_j[j], 1);

			A.apply (w, e_j);

			for (k = 0; k < n; k++)
				if (!F.areEqual (A[pair<size_t, size_t>(k, j)], w[k]))
					ret = iter_passed = false;

			report << "    Output vector " << j << ": ";
			printVector<Field> (F, report, w);
		}

		if (!iter_passed)
			report << "    ERROR: Output vectors were incorrect" << endl;
	}

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();

	return ret;
}

/* Test 5: Random apply to sparse matrix of K nonzero elements per row
 *
 * Generates a random sparse matrix with K nonzero elements per row and applies
 * it to the vectors (1,1,...,1) to test whether the output matches the sum of
 * the input's columns
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 * K - Number of nonzero elements per row
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testRandomApply3 (Field &F, size_t n, ostream &report, int iterations, int K) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	cout << "Testing sparse random apply (3)...";
	cout.flush ();
	report << "Testing sparse random apply (3):" << endl;

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	typename Field::RandIter r (F);
	typename Field::Element x, sum;

	Integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (K > n) K = n;

	Vector v(n), w(n);

	for (k = 0; k < n; k++)
		F.init (v[k], 1);

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;
		iter_passed = true;

		Blackbox A (F, n, n);

		for (j = 0; j < n; j++) {
			for (k = 0; k < K; k++) {
				pair<size_t, size_t> p (j, 0);

				do
					p.second = rand () % n;
				while (!F.isZero (A[p]));

				r.random (x);
				A.put_value (p, x);
			}
		}

		report << "    Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		A.apply (w, v);

		for (j = 0; j < n; j++) {
			F.init (sum, 0);

			for (k = 0; k < n; k++)
				F.addin (sum, A[pair<size_t, size_t>(j, k)]);

			if (!F.areEqual (sum, w[j]))
				ret = iter_passed = false;
		}

		report << "    Output vector: ";
		printVector<Field> (F, report, w);

		if (!iter_passed)
			report << "    ERROR: Output vector was incorrect" << endl;
	}

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();

	return ret;
}

/* Test 6: Minimal polynomial of the identity matrix
 *
 * Construct the identity matrix and compute its minimal polynomial. Confirm
 * that the result is 1-x
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityMinpoly (Field &F, size_t n, ostream &report) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	cout << "Testing identity minpoly...";
	cout.flush ();
	report << "Testing identity minpoly:" << endl;

	bool ret = true;
	Blackbox A (F, n, n);

	int i;
	typename Field::Element e, c0, c1;
	F.init (e, 1);

	for (i = 0; i < n; i++)
		A.put_value (pair<size_t, size_t> (i, i), e);

	Polynomial phi;

	minpoly<Field, Polynomial, Vector> (phi, A, F);

	report << "  Minimal polynomial is: ";
	printPolynomial<Field, Polynomial> (F, report, phi);

	F.init (c0, 1);
	F.init (c1, -1);

	if (phi.size () != 2 ||
	    !F.areEqual (phi[0], c0) ||
	    !F.areEqual (phi[1], c1)) {
		ret = false;
		report << "  ERROR: Minimal polynomial is incorrect" << endl;
	}

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();

	return ret;
}

/* Test 7: Minimal polynomial of a nilpotent matrix
 *
 * Construct an index-n nilpotent matrix and compute its minimal
 * polynomial. Confirm that the result is x^n
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testNilpotentMinpoly (Field &F, size_t n, ostream &report) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	cout << "Testing nilpotent minpoly...";
	cout.flush ();
	report << "Testing nilpotent minpoly:" << endl;

	bool ret = true;
	bool lowerTermsCorrect = true;
	Blackbox A (F, n, n);

	int i;
	typename Field::Element e;
	F.init (e, 1);

	for (i = 1; i < n; i++)
		A.put_value (pair<size_t, size_t> (i - 1, i), e);

	Polynomial phi;

	minpoly<Field, Polynomial, Vector> (phi, A, F);

	report << "  Minimal polynomial is: ";
	printPolynomial<Field, Polynomial> (F, report, phi);

	for (i = 0; i < n - 1; i++)
		if (!F.isZero (phi[i]))
			lowerTermsCorrect = false;

	if (phi.size () != n + 1 || !F.isOne (phi[0]) || !lowerTermsCorrect) {
		ret = false;
		report << "  ERROR: Minimal polynomial is incorrect" << endl;
	}

	if (ret) {
		cout << "passed" << endl;
		report << "Test passed" << endl << endl;
	} else {
		cout << "FAILED" << endl;
		report << "Test FAILED" << endl << endl;
	}

	cout.flush ();

	return ret;
}

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 10;
	static Integer q = 101;
	static int iterations = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",                 TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)",                 TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",                   TYPE_INT,     &iterations },
		{ 'k', "-k K", "K nonzero elements per row in sparse random apply test (default 3)", TYPE_INT,     &k },
		{ 'N', "-N N", "N nonzero elements in sparse random apply test (default 20)",        TYPE_INT,     &N }
	};

	parseArguments (argc, argv, report, args);
	ParamModular F (q);

	srand (time (NULL));

	cout << "Sparse matrix black box test suite" << endl << endl;

	if (!testIdentityApply<ParamModular>    (F, n, report, iterations)) pass = false;
	if (!testNilpotentApply<ParamModular>   (F, n, report, iterations)) pass = false;
	if (!testRandomApply1<ParamModular>     (F, n, report, iterations, k)) pass = false;
	if (!testRandomApply2<ParamModular>     (F, n, report, iterations, N)) pass = false;
	if (!testRandomApply3<ParamModular>     (F, n, report, iterations, k)) pass = false;
	if (!testIdentityMinpoly<ParamModular>  (F, n, report)) pass = false;
	if (!testNilpotentMinpoly<ParamModular> (F, n, report)) pass = false;

	return pass ? 0 : -1;
}
