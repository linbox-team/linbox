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
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/large-modular.h"
#include "linbox/blackbox/sparse-matrix.h"

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
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityApply (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	commentator.start ("Testing identity apply...", "testIdentityApply", iterations);

	bool ret = true;
	bool iter_passed = true;
	Blackbox A (F, n, n);

	int i, j;
	typename Field::element e;
	F.init (e, 1);

	for (i = 0; i < n; i++)
		A.put_value (pair<size_t, size_t> (i, i), e);

	Vector v(n), w(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		for (j = 0; j < n; j++)
			r.random (v[j]);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		printVector<Field> (F, report, v);

		A.apply (w, v);

		commentator.indent (report);
		report << "Output vector: ";
		printVector<Field> (F, report, w);

		for (j = 0; j < n; j++)
			if (!F.areEqual (w[j], v[j]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityApply");

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
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testNilpotentApply (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	commentator.start ("Testing nilpotent apply...", "testNilpotentApply", iterations);

	bool ret = true;
	bool iter_passed;
	bool all_zero;
	Blackbox A (F, n, n);

	int i, j;
	typename Field::element e;
	F.init (e, 1);
	bool even = false;

	for (i = 1; i < n; i++)
		A.put_value (pair<size_t, size_t> (i - 1, i), e);

	Vector v(n), w(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;
		even = false;

		for (j = 0; j < n; j++)
			do r.random (v[j]); while (F.isZero (v[j]));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		printVector<Field> (F, report, v);

		for (j = 0; j < n - 1; j++, even = !even)
			if (even)
				A.apply (v, w);
			else
				A.apply (w, v);

		commentator.indent (report);
		report << "A^(n-1) v:     ";
		printVector<Field> (F, report, even ? w : v);

		for (all_zero = true, j = 0; j < n; j++)
			if (!F.isZero (even ? w[j] : v[j]))
				all_zero = false;

		if (all_zero) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A^(n-1) v is prematurely zero" << endl;
		}

		if (even)
			A.apply (v, w);
		else
			A.apply (w, v);

		commentator.indent (report);
		report << "A^n v:         ";
		printVector<Field> (F, report, even ? v : w);

		for (j = 0; j < n; j++)
			if (!F.isZero (even ? v[j] : w[j]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A^n v is non-zero" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNilpotentApply");

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
 * K - Number of nonzero elements per row
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testRandomApply1 (Field &F, size_t n, int iterations, int K) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	commentator.start ("Testing sparse random apply (1)...", "testRandomApply1", iterations);

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	typename Field::RandIter r (F);
	typename Field::element x;

	integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (K > n) K = n;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

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

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		for (j = 0; j < n; j++) {
			Vector e_j(n), w(n);

			F.init (e_j[j], 1);

			A.apply (w, e_j);

			for (k = 0; k < n; k++)
				if (!F.areEqual (A[pair<size_t, size_t>(k, j)], w[k]))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << j << ": ";
			printVector<Field> (F, report, w);
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vectors were incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply1");

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
 * N - Number of nonzero elements
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testRandomApply2 (Field &F, size_t n, int iterations, int N) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	commentator.start ("Testing sparse random apply (2)...", "testRandomApply2", iterations);

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	typename Field::RandIter r (F);
	typename Field::element x;

	integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (N > n * n) N = n * n;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

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

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		for (j = 0; j < n; j++) {
			Vector e_j(n), w(n);

			F.init (e_j[j], 1);

			A.apply (w, e_j);

			for (k = 0; k < n; k++)
				if (!F.areEqual (A[pair<size_t, size_t>(k, j)], w[k]))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << j << ": ";
			printVector<Field> (F, report, w);
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vectors were incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply2");

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
 * K - Number of nonzero elements per row
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testRandomApply3 (Field &F, size_t n, int iterations, int K) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef SparseMatrix <Field, Row, Vector> Blackbox;

	commentator.start ("Testing sparse random apply (3)...", "testRandomApply3", iterations);

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	typename Field::RandIter r (F);
	typename Field::element x, sum;

	integer c;
	long width;

	F.characteristic (c);
	width = logp (c, 10) + 1;

	if (K > n) K = n;

	Vector v(n), w(n);

	for (k = 0; k < n; k++)
		F.init (v[k], 1);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

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

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.prettyPrint (report, 6, width);

		A.apply (w, v);

		for (j = 0; j < n; j++) {
			F.init (sum, 0);

			for (k = 0; k < n; k++)
				F.addin (sum, A[pair<size_t, size_t>(j, k)]);

			if (!F.areEqual (sum, w[j]))
				ret = iter_passed = false;
		}

		commentator.indent (report);
		report << "Output vector: ";
		printVector<Field> (F, report, w);

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector was incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply3");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;
	static int iterations = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",                 TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 65521)",               TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",                   TYPE_INT,     &iterations },
		{ 'k', "-k K", "K nonzero elements per row in sparse random apply test (default 3)", TYPE_INT,     &k },
		{ 'N', "-N N", "N nonzero elements in sparse random apply test (default 20)",        TYPE_INT,     &N }
	};

	parseArguments (argc, argv, args);
	LargeModular F (q);

	srand (time (NULL));

	cout << "Sparse matrix black box test suite" << endl << endl;

	if (!testIdentityApply<LargeModular>    (F, n, iterations)) pass = false;
	if (!testNilpotentApply<LargeModular>   (F, n, iterations)) pass = false;
	if (!testRandomApply1<LargeModular>     (F, n, iterations, k)) pass = false;
	if (!testRandomApply2<LargeModular>     (F, n, iterations, N)) pass = false;
	if (!testRandomApply3<LargeModular>     (F, n, iterations, k)) pass = false;

	return pass ? 0 : -1;
}
