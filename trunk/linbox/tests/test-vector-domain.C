/* -*- mode: c; style: linux -*- */

/* tests/test-matrix-domain.C
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
#include "linbox/field/matrix-domain.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Dot product of dense vectors
 *
 * Construct two random dense vectors and compute their doc product
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * report - Stream to which to output detailed report of failures, if any
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDenseDotProduct (Field &F, size_t n, ostream &report, int iterations) 
{
	typedef vector <typename Field::element> Vector;

	cout << "Testing dense/dense dot product...";
	cout.flush ();
	report << "Testing dense/dense dot product:" << endl;

	bool ret = true;

	Vector v1 (n), v2 (n);
	typename Field::element sigma, rho;
	typename Field::RandIter r (F);

	MatrixDomain<Field, Vector, Vector> MD (F);

	int i, j;

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;

		F.init (sigma, 0);

		for (j = 0; j < n; j++) {
			r.random (v1[j]);
			r.random (v2[j]);
			F.axpyin (sigma, v1[j], v2[j]);
		}

		report << "    Input vector 1:  ";
		printVector<Field> (F, report, v1);

		report << "    Input vector 2:  ";
		printVector<Field> (F, report, v2);

		MD.dotprod (rho, v1, v2);

		report << "    True dot product: ";
		F.write (report, sigma);
		report << endl;

		report << "    Dot product from matrix domain: ";
		F.write (report, rho);
		report << endl;

		if (!F.areEqual (sigma, rho)) {
			ret = false;
			report << "    ERROR: Dot products are not equal" << endl;
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

/* Test 1: Dot product dense vector and sparse vector
 *
 * Construct a random dense vector and a random sparse vector and compute their
 * doc product
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * report - Stream to which to output detailed report of failures, if any
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDenseSparseDotProduct (Field &F, size_t n, ostream &report, int iterations) 
{
	typedef vector <pair <size_t, typename Field::element> > Vector1;
	typedef vector <typename Field::element> Vector2;

	cout << "Testing dense/sparse dot product...";
	cout.flush ();
	report << "Testing dense/sparse dot product:" << endl;

	bool ret = true;

	Vector1 v1 (n);
	Vector2 v2 (n);
	typename Field::element sigma, rho, tmp;
	typename Field::RandIter r (F);

	MatrixDomain<Field, Vector1, Vector2> MD (F);

	int i, j;

	for (i = 0; i < iterations; i++) {
		report << "  Iteration " << i << ": " << endl;

		F.init (sigma, 0);
		v1.clear ();

		for (j = 0; j < n; j++) {
			r.random (v2[j]);

			// Give the sparse vector an entry about 10% of the time
			if (rand () % 100 < 10) {
				r.random (tmp);
				v1.push_back (pair <size_t, typename Field::element> (j, tmp));
				F.axpyin (sigma, v2[j], tmp);
			}
		}

		report << "    Input vector 1:  ";
		printSparseSeqVector<Field> (F, report, v1);

		report << "    Input vector 2:  ";
		printVector<Field> (F, report, v2);

		MD.dotprod (rho, v1, v2);

		report << "    True dot product: ";
		F.write (report, sigma);
		report << endl;

		report << "    Dot product from matrix domain: ";
		F.write (report, rho);
		report << endl;

		if (!F.areEqual (sigma, rho)) {
			ret = false;
			report << "    ERROR: Dot products are not equal" << endl;
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

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 100;
	static integer q = 101;
	static int iterations = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to N (default 100)",          TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 4294967291)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, report, args);
	ParamModular F (q);

	srand (time (NULL));

	cout << "Matrix domain test suite" << endl << endl;

	if (!testDenseDotProduct<ParamModular>       (F, n, report, iterations)) pass = false;
	if (!testDenseSparseDotProduct<ParamModular> (F, n, report, iterations)) pass = false;

	return pass ? 0 : -1;
}
