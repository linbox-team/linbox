/* -*- mode: C++; style: linux -*- */

/* tests/test-dense-matrix.C
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

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/dense-matrix1.h"
#include "linbox/solutions/minpoly.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Identity matrix in dense representation
 *
 * Construct a dense representation of an n x n identity matrix and check
 * whether the output of its application to a series of random vectors is equal
 * to the input.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply identity inverse
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentity (Field &F, long n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing identity apply", "testIdentity", iterations);

	bool ret = true;
	bool iter_passed = true;

	int i, j;

	Blackbox I (F, n, n);
	typename Field::Element one;

	F.init (one, 1);

	for (i = 0; i < n; i++)
		I.setEntry (i, i, one);

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
		report << "Input vector: ";
		printVector<Field> (F, report, v);

		I.apply (w, v);

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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentity");

	return ret;
}

/* Test 2: Application of Vandermonde matrix in dense representation
 *
 * Computes a random Vandermonde matrix and applies it to a series of random
 * vectors. The random vectors contain the coefficients of polynomials over the
 * ground field. The output of the application is the result of evaluating these
 * polynomials at the points given by the second column of the matrix. This
 * function interpolates (using Lagrange interpolants) the evaluation points to
 * get the original polynomials and checks whether the coefficients match the
 * original vectors.
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 * N - Number of random vectors to which to apply random Vandermonde matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testVandermonde (Field &F, long n, int iterations, int N) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing Vandermonde apply", "testVandermonde", iterations);

	bool ret = true;
	bool inner_iter_passed;

	int i, j, k;

	Blackbox V (F, n, n);

	Vector x(n), v(n), y(n), f(n);
	typename Field::RandIter r (F);
	typename Field::Element t;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		/* Evaluation points */
		for (j = 0; j < n; j++) {
			bool flag = true;

			// Make sure points are all distinct
			while (flag) {
				r.random (x[j]);
				flag = false;
				for (k = 0; k < j; k++)
					if (F.areEqual (x[j], x[k]))
						flag = true;
			}
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Evaluation points: ";
		printVector<Field> (F, report, x);

		/* Build the Vandermonde matrix */
		for (j = 0; j < n; j++) {
			F.init (t, 1);

			for (k = 0; k < n; k++) {
				V.setEntry (j, k, t);
				F.mulin (t, x[j]);
			}
		}

		for (j = 0; j < N; j++) {
			inner_iter_passed = true;

			/* Random vector of evaluation results */
			for (k = 0; k < n; k++)
				r.random (v[k]);

			commentator.indent (report);
			report << "Input vector: ";
			printVector<Field> (F, report, v);

			/* w should now be a vector of polynomial evaluations */
			V.apply (y, v);

			commentator.indent (report);
			report << "Output vector: ";
			printVector<Field> (F, report, y);

			/* Polynomial interpolation to check whether w is correct */
			interpolatePoly (F, f, x, y);

			commentator.indent (report);
			report << "Interpolation results: ";
			printVector<Field> (F, report, f);

			for (k = 0; k < n; k++)
				if (!F.areEqual (f[k], v[k]))
					ret = inner_iter_passed = false;

			if (!inner_iter_passed)
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Vectors are not equal" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testVandermonde");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;
	static int iterations = 100;
	static int N = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",   TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Dense matrix black box test suite" << endl << endl;

	if (!testIdentity<Modular<long> >    (F, n, iterations)) pass = false;
	if (!testVandermonde<Modular<long> > (F, n, iterations, N)) pass = false;

	return pass ? 0 : -1;
}
