/* -*- mode: c; style: linux -*- */

/* tests/test-det.C
 * Copyright (C) 2002 Bradford Hovinen
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
	typedef vector <typename Field::element> Vector;
	typedef vector <typename Field::element> Polynomial;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing nonsingular diagonal determinant (1)", "testDiagonalDet1", iterations);

	bool ret = true;
	bool done;
	int i, j, k;

	Vector d(n);
	typename Field::element pi, phi;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

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
	typedef vector <typename Field::element> Vector;
	typedef vector <typename Field::element> Polynomial;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing nonsingular diagonal determinant (2)", "testDiagonalDet2", iterations);

	bool ret = true;
	int i, j, k;

	Vector d(n);
	typename Field::element pi, phi;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

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
	typedef vector <typename Field::element> Vector;
	typedef vector <typename Field::element> Polynomial;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing singular diagonal determinant", "testSingularDiagonalDet", iterations);

	bool ret = true;
	int i, j, k;

	Vector d(n);
	typename Field::element phi;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

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
	static int numVectors = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",    TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Black box determinant test suite" << endl << endl;

	if (!testDiagonalDet1<Modular<long> >        (F, n, iterations)) pass = false;
	if (!testDiagonalDet2<Modular<long> >        (F, n, iterations)) pass = false;
	if (!testSingularDiagonalDet<Modular<long> > (F, n, iterations)) pass = false;

	return pass ? 0 : -1;
}
