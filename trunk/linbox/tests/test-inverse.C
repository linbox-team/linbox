/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/hilbert.h"
#include "linbox/blackbox/dense-matrix.h"
#include "linbox/blackbox/inverse.h"
#include "linbox/util/vector-factory.h"

#include "linbox/solutions/minpoly.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Inverse of the identity matrix
 *
 * Constructs a black box for the inverse of an n x n identity matrix and checks
 * that that inverse is itself the identity operator by applying it to random
 * vectors.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply identity inverse
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityInverse (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing identity inverse", "testIdentityInverse", iterations);

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
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		for (j = 0; j < n; j++)
			r.random (v[j]);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		printVector<Field> (F, report, v);

		DT.apply (w, v);

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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityInverse");

	return ret;
}

/* Test 2: Inverse of Hilbert matrix
 *
 * Constructs an n x n Hilbert matrix and a black box for its inverse. Applies
 * each to random vectors and checks that the results are equal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testHilbertInverse (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Hilbert <Field, Vector> Blackbox;

	commentator.start ("Testing Hilbert inverse", "testHilbertInverse", iterations);

	bool ret = true;
	bool iter_passed;

	int i, j;

	Blackbox H (F, n);
	Inverse<Field, Vector> HT (F, &H);

	Vector v(n), w(n), z(n);
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

		H.apply (z, v);
		HT.apply (w, z);

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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testHilbertInverse");

	return ret;
}

/* Test 3: Inverse of Vandermonde matrix
 *
 * Computes a random Vandermonde matrix and its inverse. This inverse is a
 * linear operator that interpolates the values given in the input vector to
 * produce a polynomial whose coefficients are the elements of the output
 * vector. We then evaluate the polynomial in Horner fashion at each of the
 * evaluation points generated above to check whether the result is the original
 * input vector.
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 * N - Number of random vectors to which to apply random Vandermonde matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testVandermondeInverse (Field &F, size_t n, int iterations, int N) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing Vandermonde inverse", "testVandermondeInverse", iterations);

	bool ret = true;
	bool inner_iter_passed;

	int i, j, k;

	Blackbox V (F, n, n);

	Vector x(n), v(n), w(n), z(n);
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
		report.flush ();

		/* Build the Vandermonde matrix */
		for (j = 0; j < n; j++) {
			F.init (t, 1);

			for (k = 0; k < n; k++) {
				V.setEntry (j, k, t);
				F.mulin (t, x[j]);
			}
		}

		Inverse<Field, Vector> VT (F, &V);

		for (j = 0; j < N; j++) {
			inner_iter_passed = true;

			/* Random vector of evaluation results */
			for (k = 0; k < n; k++)
				r.random (v[k]);

			commentator.indent (report);
			report << "Input vector: ";
			printVector<Field> (F, report, v);

			/* w should now be the requisite polynomial */
			VT.apply (w, v);

			commentator.indent (report);
			report << "Output vector: ";
			printVector<Field> (F, report, w);

			/* Multipoint evaluation to check whether w is correct */
			multiEvalPoly (F, z, w, x);

			commentator.indent (report);
			report << "Evaluation results: ";
			printVector<Field> (F, report, z);

			for (k = 0; k < n; k++)
				if (!F.areEqual (z[k], v[k]))
					ret = inner_iter_passed = false;

			if (!inner_iter_passed)
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Vectors are not equal" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testVandermondeInverse");

	return ret;
}

/* Test 3: Inverse of diagonal inverse
 *
 * Constructs a random nonsingular diagonal matrix and its inverse, and extracts
 * the values along the diagonal of the inverse. Checks that those values are in
 * fact the inverses of the diagonal elements in the original.
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalInverse (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing diagonal inverse", "testDiagonalInverse", iterations);

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	Vector d(n), di(n), dt(n), e(n), DTe(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		for (j = 0; j < n; j++) {
			do r.random (d[j]); while (F.isZero (d[j]));
			F.inv (di[j], d[j]);
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		printVector<Field> (F, report, d);

		commentator.indent (report);
		report << "Expeted diagonal entries of inverse: ";
		printVector<Field> (F, report, di);

		Blackbox D (F, d);
		Inverse <Field, Vector> DT (F, &D);

		for (j = 0; j < n; j++) {
			F.init (e[j], 1);
			DT.apply (DTe, e);
			dt[j] = DTe[j];
		}

		commentator.indent (report);
		report << "Diagonal entries of computed inverse: ";
		printVector<Field> (F, report, dt);

		for (j = 0; j < n; j++)
			if (!F.areEqual (di[j], dt[j]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed inverse does not match expected inverse" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalInverse");

	return ret;
}

/* Test 3: Random transpose
 *
 * Compute the inverse of a random dense matrix and apply its transpose to
 * random vectors
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply matrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testRandomTranspose (Field &F,
				 VectorFactory<Vector> &factory1,
				 VectorFactory<Vector> &factory2) 
{
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing random transpose", "testRandomTranspose", factory1.m ());

	int i, j;
	typename Field::Element x;
	typename Field::RandIter r (F);

	Blackbox A (F, factory1.n (), factory2.n ());
	Inverse<Field, Vector> Ainv (F, &A);

	for (i = 0; i < factory1.n (); i++) {
		for (j = 0; j < factory2.n (); j++) {
			r.random (x);
			A.setEntry (i, j, x);
		}
	}

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	bool ret = testTranspose (F, Ainv, factory1, factory2);

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Black box inverse test suite" << endl << endl;

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	RandomDenseVectorFactory<Modular<long> > factory1 (F, n, iterations), factory2 (F, n, iterations);

	if (!testIdentityInverse<Modular<long> >    (F, n, iterations)) pass = false;
	if (!testVandermondeInverse<Modular<long> > (F, n, iterations, N)) pass = false;
	if (!testDiagonalInverse<Modular<long> >    (F, n, iterations)) pass = false;
	if (!testRandomTranspose<Modular<long> >    (F, factory1, factory2)) pass = false;

	return pass ? 0 : -1;
}
