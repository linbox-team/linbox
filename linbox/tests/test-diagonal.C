/* -*- mode: c; style: linux -*- */

/* tests/test-diagonal.C
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

#include "linbox/field/archetype.h"
#include "linbox/field/large-modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/minpoly.h"

#include "test-common.h"
#include "test-generic.h"

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
	typedef vector <typename Field::element> Vector;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	cout << "Testing identity apply...";
	cout.flush ();
	report << "Testing identity apply:" << endl;

	bool ret = true;
	bool iter_passed = true;

	Vector d(n);

	int i, j;

	for (i = 0; i < n; i++)
		F.init (d[i], 1);

	Blackbox D (F, d);

	Vector v(n), w(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;
		iter_passed = true;

		for (j = 0; j < n; j++)
			r.random (v[j]);

		report << "    Input vector:  ";
		printVector<Field> (F, report, v);

		D.apply (w, v);

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

/* Test 2: Constant term in minimal polynomial of diagonal map
 *
 * Generates a random diagonal nonsingular matrix and computes its minimal
 * polynomial. Checks that the constant term thereof equals the product of the
 * entries on the diagonal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomMinpoly (Field &F, size_t n, ostream &report, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <typename Field::element> Polynomial;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	cout << "Testing random minpoly...";
	cout.flush ();
	report << "Testing random minpoly:" << endl;

	bool ret = true;

	int i, j;
	typename Field::element pi;
	Polynomial m_D;

	Vector d(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;

		F.init (pi, 1);
		for (j = 0; j < n; j++) {
			do r.random (d[j]); while (F.isZero (d[j]));
			F.mulin (pi, d[j]);
		}

		report << "    Diagonal vector: ";
		printVector<Field> (F, report, d);

		report << "    Product: ";
		F.write (report, pi);
		report << endl;

		Blackbox D (F, d);
		minpoly<Field, Polynomial, Vector> (m_D, D, F);

		report << "    Minimal polynomial: ";
		printPolynomial<Field, Polynomial> (F, report, m_D);

		if (!F.areEqual (m_D[0], pi)) {
			report << "    ERROR: m_D(0) != det(D)" << endl;
			ret = false;
		}
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

/* Test 3: Random transpose
 *
 * Compute a random diagonal matrix and use the transpose test in test-generic.h
 * to check consistency of transpose apply.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * report - Stream to which to output detailed report of failures, if any
 * iterations - Number of random vectors to which to apply matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomTranspose (Field &F, size_t n, ostream &report, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef Diagonal <Field, Vector> Blackbox;

	Vector d(n);
	typename Field::RandIter r (F);

	for (int i = 0; i < n; i++)
		r.random (d[i]);

	Blackbox D (F, d);

	return testTranpose<Field> (F, D, report, iterations);
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

	cout << "Diagonal matrix black box test suite" << endl << endl;

	if (!testIdentityApply<LargeModular>    (F, n, report, iterations)) pass = false;
	if (!testRandomMinpoly<LargeModular>    (F, n, report, iterations)) pass = false;
	if (!testRandomTranspose<LargeModular>  (F, n, report, iterations)) pass = false;

	return pass ? 0 : -1;
}
