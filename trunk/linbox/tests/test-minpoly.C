/* -*- mode: c; style: linux -*- */

/* tests/test-minpoly.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * 2002-04-03: William J. Turner <wjturner@acm.org>
 *
 * changed name of sparse-matrix file.
 *
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse0.h"	// wjt: was sparse-matrix.h
#include "linbox/solutions/minpoly.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Minimal polynomial of the identity matrix
 *
 * Construct the identity matrix and compute its minimal polynomial. Confirm
 * that the result is 1-x
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityMinpoly (Field &F, size_t n) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing identity minpoly", "testIdentityMinpoly");

	bool ret = true;
	Blackbox A (F, n, n);

	int i;
	typename Field::Element e, c0, c1;
	F.init (e, 1);

	for (i = 0; i < n; i++)
		A.put_value (pair<size_t, size_t> (i, i), e);

	Polynomial phi;

	minpoly<Field, Polynomial, Vector> (phi, A, F);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Minimal polynomial is: ";
	printPolynomial<Field, Polynomial> (F, report, phi);

	F.init (c0, -1);
	F.init (c1, 1);

	if (phi.size () != 2 ||
	    !F.areEqual (phi[0], c0) ||
	    !F.areEqual (phi[1], c1)) {
		ret = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Minimal polynomial is incorrect" << endl;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityMinpoly");

	return ret;
}

/* Test 2: Minimal polynomial of a nilpotent matrix
 *
 * Construct an index-n nilpotent matrix and compute its minimal
 * polynomial. Confirm that the result is x^n
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testNilpotentMinpoly (Field &F, size_t n) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing nilpotent minpoly", "testNilpotentMinpoly");

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

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Minimal polynomial is: ";
	printPolynomial<Field, Polynomial> (F, report, phi);

	for (i = 0; i < n - 1; i++)
		if (!F.isZero (phi[i]))
			lowerTermsCorrect = false;

	if (phi.size () != n + 1 || !F.isOne (phi[n]) || !lowerTermsCorrect) {
		ret = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Minimal polynomial is incorrect" << endl;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNilpotentMinpoly");

	return ret;
}

/* Test 3: Random minpoly of sparse matrix of K nonzero elements per row
 *
 * Generates a random sparse matrix with K nonzero elements per row and computes
 * its minimal polynomial. Then computes random vectors and applies the
 * polynomial to them in Horner style, checking whether the result is 0.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * K - Number of nonzero elements per row
 * numVectors - Number of random vectors to which to apply the minimal polynomial
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testRandomMinpoly1 (Field &F, size_t n, int iterations, int K, int numVectors) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing sparse random minpoly (1)", "testRandomMinpoly1", iterations);

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	typename Field::RandIter r (F);
	typename Field::Element x;

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

		Polynomial phi;

		minpoly<Field, Polynomial, Vector> (phi, A, F);

		commentator.indent (report);
		report << "Minimal polynomial is: ";
		printPolynomial<Field, Polynomial> (F, report, phi);

		for (j = 0; j < numVectors; j++) {
			Vector v(n), w(n);
			bool even = false;

			for (k = 0; k < n; k++)
				r.random (v[k]);

			applyPoly<Field, Polynomial> (F, w, A, phi, v);

			for (k = 0; k < n; k++)
				if (!F.isZero (w[k]))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << j << ": ";
			printVector<Field> (F, report, w);
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector was incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomMinpoly1");

	return ret;
}

/* Test 4: Random minpoly of sparse matrix of N nonzero elements
 *
 * Generates a random sparse matrix with K nonzero elements per row and computes
 * its minimal polynomial. Then computes random vectors and applies the
 * polynomial to them in Horner style, checking whether the result is 0.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * N - Number of nonzero elements
 * numVectors - Number of random vectors to which to apply the minimal polynomial
 *
 * Return true on success and false on failure
 */

template <class Field>
bool testRandomMinpoly2 (Field &F, size_t n, int iterations, int N, int numVectors) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing sparse random minpoly (2)", "testRandomMinpoly2", iterations);

	bool ret = true;
	bool iter_passed;

	int i, j, k;

	typename Field::RandIter r (F);
	typename Field::Element x;

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

		Polynomial phi;

		minpoly<Field, Polynomial, Vector> (phi, A, F);

		commentator.indent (report);
		report << "Minimal polynomial is: ";
		printPolynomial<Field, Polynomial> (F, report, phi);

		for (j = 0; j < numVectors; j++) {
			Vector v(n), w(n);
			bool even = false;

			for (k = 0; k < n; k++)
				r.random (v[k]);

			applyPoly<Field, Polynomial> (F, w, A, phi, v);

			for (k = 0; k < n; k++)
				if (!F.isZero (w[k]))
					ret = iter_passed = false;

			commentator.indent (report);
			report << "Output vector " << j << ": ";
			printVector<Field> (F, report, w);
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector was incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomMinpoly2");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 10;
	static int numVectors = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",                 TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)",          TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",                    TYPE_INT,     &iterations },
		{ 'v', "-v V", "Use V test vectors for the random minpoly tests (default 100)",      TYPE_INT,     &numVectors },
		{ 'k', "-k K", "K nonzero Elements per row in sparse random apply test (default 3)", TYPE_INT,     &k },
		{ 'N', "-N N", "N nonzero Elements in sparse random apply test (default 20)",        TYPE_INT,     &N }
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);

	cout << "Black box minimal polynomial test suite" << endl << endl;

	if (!testIdentityMinpoly<Modular<long> >  (F, n)) pass = false;
	if (!testNilpotentMinpoly<Modular<long> > (F, n)) pass = false;
	if (!testRandomMinpoly1<Modular<long> >   (F, n, iterations, k, numVectors)) pass = false;
	if (!testRandomMinpoly2<Modular<long> >   (F, n, iterations, N, numVectors)) pass = false;

	return pass ? 0 : -1;
}
