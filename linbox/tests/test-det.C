/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-det.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/det.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Determinant of nonsingular diagonal matrix with distinct entries
 *
 * Construct a random diagonal matrix with distinct entries and see that its
 * computed determinant equals the product of diagonal entries
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalDet1 (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing nonsingular diagonal determinant (1)", "testDiagonalDet1", iterations);

	bool ret = true;
	bool done;
	int i;
	size_t j, k;

	VectorDomain<Field> VD (F);

	Vector d(n);
	typename Field::Element pi, phi;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		F.init (pi, 1);

		for (j = 0; j < n; j++) {
			do {
				do r.random (d[j]); while (F.isZero (d[j]));

				done = true;
				for (k = 0; k < j; k++) {
					if (F.areEqual (d[j], d[k])) {
						done = false;
						break;
					}
				}
			} while (!done);

			F.mulin (pi, d[j]);
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		commentator.indent (report);
		report << "True determinant: ";
		F.write (report, pi);
		report << endl;

		Blackbox D (F, d);

		det (phi, D, F);

		commentator.indent (report);
		report << "Computed determinant: ";
		F.write (report, phi);
		report << endl;

		if (!F.areEqual (pi, phi)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed determinant is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalDet1");

	return ret;
}

/* Test 2: Determinant of nonsingular diagonal matrix with nondistinct entries
 *
 * Construct a random diagonal matrix with nondistinct entries and see that its
 * computed determinant equals the product of diagonal entries
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalDet2 (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing nonsingular diagonal determinant (2)", "testDiagonalDet2", iterations);

	bool ret = true;
	int i, k;
	size_t j;

	Vector d(n);
	typename Field::Element pi, phi;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		F.init (pi, 1);

		for (j = 0; j < n / 2; j++) {
			do r.random (d[j]); while (F.isZero (d[j]));
			F.mulin (pi, d[j]);
		}

		for (j = n / 2; j < n; j++) {
			k = rand () % (n / 2);
			d[j] = d[k];
			F.mulin (pi, d[j]);
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		printVector<Field> (F, report, d);

		commentator.indent (report);
		report << "True determinant: ";
		F.write (report, pi);
		report << endl;

		Blackbox D (F, d);

		det <Field, Vector> (phi, D, F);

		commentator.indent (report);
		report << "Computed determinant: ";
		F.write (report, phi);
		report << endl;

		if (!F.areEqual (pi, phi)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed determinant is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalDet2");

	return ret;
}

/* Test 3: Determinant of singular diagonal matrix
 *
 * Construct a random diagonal matrix with one zero entry
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testSingularDiagonalDet (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing singular diagonal determinant", "testSingularDiagonalDet", iterations);

	bool ret = true;
	int i;
	size_t j;

	Vector d(n);
	typename Field::Element phi;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		for (j = 0; j < n; j++)
			r.random (d[j]);

		F.init (d[rand () % n], 0);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		printVector<Field> (F, report, d);

		Blackbox D (F, d);

		det <Field, Vector> (phi, D, F);

		commentator.indent (report);
		report << "Computed determinant: ";
		F.write (report, phi);
		report << endl;

		if (!F.isZero (phi)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed determinant is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSingularDiagonalDet");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101U;
	static int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",    TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<uint32> F (q);

	srand (time (NULL));

	cout << "Black box determinant test suite" << endl << endl;

	if (!testDiagonalDet1        (F, n, iterations)) pass = false;
	if (!testDiagonalDet2        (F, n, iterations)) pass = false;
	if (!testSingularDiagonalDet (F, n, iterations)) pass = false;

	return pass ? 0 : -1;
}
