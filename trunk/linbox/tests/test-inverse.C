/* -*- mode: c; style: linux -*- */

/* tests/test-inverse.C
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

#include "linbox/field/large-modular.h"

#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/hilbert.h"
#include "linbox/blackbox/inverse.h"

#include "linbox/solutions/minpoly.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Inverse of the identity matrix
 *
 * Constructs a black box for the inverse of an n x n identity matrix and checks
 * that that inverse is itself the identity operator by applying it to random
 * vectors.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 * iterations - Number of random vectors to which to apply identity inverse
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityInverse (Field &F, size_t n, ostream &report, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	cout << "Testing identity inverse...";
	cout.flush ();
	report << "Testing identity inverse:" << endl;

	bool ret = true;
	bool iter_passed = true;

	Vector d(n);

	int i, j;

	for (i = 0; i < n; i++)
		F.init (d[i], 1);

	Blackbox D (F, d);
	Inverse<Field, Vector> DT (F, &D);

	Vector v(n), w(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;
		iter_passed = true;

		for (j = 0; j < n; j++)
			r.random (v[j]);

		report << "    Input vector:  ";
		printVector<Field> (F, report, v);

		DT.apply (w, v);

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

/* Test 2: Inverse of hilbert matrix
 *
 * Constructs an n x n Hilbert matrix and a black box for its inverse. Applies
 * each to random vectors and checks that the results are equal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testHilbertInverse (Field &F, size_t n, ostream &report, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <typename Field::element> Polynomial;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Hilbert <Field, Vector> Blackbox;

	test_header("Hilbert inverse", report);
	//cout << "Testing Hilbert inverse...";
	//cout.flush ();
	//report << "Testing Hilbert inverse:" << endl;

	bool ret = true;
	bool iter_passed;

	int i, j;

	Blackbox H (F, n);
	Inverse<Field, Vector> HT (F, &H);

	Vector v(n), w(n), z(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;
		iter_passed = true;

		for (j = 0; j < n; j++)
			r.random (v[j]);

		report << "    Input vector: ";
		printVector<Field> (F, report, v);

		H.apply (z, v);
		HT.apply (w, z);

		report << "    Output vector: ";
		printVector<Field> (F, report, w);

		for (j = 0; j < n; j++)
			if (!F.areEqual (w[j], v[j]))
				ret = iter_passed = false;

		if (!iter_passed)
			report << "    ERROR: Vectors are not equal" << endl;
	}

	//if (ret) {
	//	cout << "passed" << endl;
	//	report << "Test passed" << endl << endl;
	//} else {
	//	cout << "FAILED" << endl;
	//	report << "Test FAILED" << endl << endl;
	//}
	//
	//cout.flush ();

	return test_trailer(ret, report);
}

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 10;
	static integer q = 4294967291U;
	static int iterations = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 4294967291)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, report, args);
	LargeModular F (q);

	srand (time (NULL));

	cout << "Black box inverse test suite" << endl << endl;

	if (!testIdentityInverse<LargeModular>   (F, n, report, iterations)) pass = false;
	if (!testHilbertInverse<LargeModular>    (F, n, report, iterations)) pass = false;

	return pass ? 0 : -1;
}
