/* -*- mode: c; style: linux -*- */

/* tests/test-vector-domain.C
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
#include "linbox/field/modular.h"
#include "linbox/field/vector-domain.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Dot product of dense vectors
 *
 * Construct two random dense vectors and compute their doc product
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDenseDotProduct (Field &F, long n, int iterations) 
{
	typedef vector <typename Field::element> Vector;

	commentator.start ("Testing dense/dense dot product", "testDenseDotProduct", iterations);

	bool ret = true;

	Vector v1 (n), v2 (n);
	typename Field::element sigma, rho;
	typename Field::RandIter r (F);

	VectorDomain<Field, Vector, Vector> VD (F);

	int i, j;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		F.init (sigma, 0);

		for (j = 0; j < n; j++) {
			r.random (v1[j]);
			r.random (v2[j]);
			F.axpyin (sigma, v1[j], v2[j]);
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		VD.dotprod (rho, v1, v2);

		commentator.indent (report);
		report << "True dot product: ";
		F.write (report, sigma);
		report << endl;

		commentator.indent (report);
		report << "Dot product from vector domain: ";
		F.write (report, rho);
		report << endl;

		if (!F.areEqual (sigma, rho)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Dot products are not equal" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDenseDotProduct");

	return ret;
}

/* Test 2: Dot product dense vector and sparse vector
 *
 * Construct a random dense vector and a random sparse vector and compute their
 * doc product
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDenseSparseDotProduct (Field &F, long n, int iterations) 
{
	typedef vector <pair <size_t, typename Field::element> > Vector1;
	typedef vector <typename Field::element> Vector2;

	commentator.start ("Testing dense/sparse dot product", "testDenseSparseDotProduct", iterations);

	bool ret = true;

	Vector1 v1 (n);
	Vector2 v2 (n);
	typename Field::element sigma, rho, tmp;
	typename Field::RandIter r (F);

	VectorDomain<Field, Vector1, Vector2> VD (F);

	int i, j;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

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

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		VD.dotprod (rho, v1, v2);

		commentator.indent (report);
		report << "True dot product: ";
		F.write (report, sigma);
		report << endl;

		commentator.indent (report);
		report << "Dot product from vector domain: ";
		F.write (report, rho);
		report << endl;

		if (!F.areEqual (sigma, rho)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Dot products are not equal" << endl;
		}

		commentator.stop ("done");

		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDenseDotProduct");

	return ret;
}

/* Test 3: Vector-vector axpy, dense vectors
 *
 * Construct two random dense vectors x and y and a random element a and compute
 * (x + a*y) - a*(y + a^-1*x). Check whether the result is 0.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDenseAXPY (Field &F, long n, int iterations) 
{
	typedef vector <typename Field::element> Vector;

	commentator.start ("Testing dense vector axpy", "testDenseAXPY", iterations);

	bool ret = true;
	bool iter_passed;

	Vector v1 (n);
	Vector v2 (n);
	Vector v3 (n);
	Vector v4 (n);
	typename Field::element a;
	typename Field::element ainv;
	typename Field::element aneg;
	typename Field::RandIter r (F);

	VectorDomain<Field, Vector, Vector> VD (F);

	int i, j;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		for (j = 0; j < n; j++) {
			r.random (v2[j]);
			r.random (v1[j]);
		}

		do r.random (a); while (F.isZero (a));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		commentator.indent (report);
		report << "Element a:  ";
		F.write (report, a);
		report << endl;

		F.inv (ainv, a);
		F.neg (aneg, a);
		VD.axpy (v3, v1, a, v2);
		VD.axpy (v4, v2, ainv, v1);
		VD.axpyin (v3, aneg, v4);

		commentator.indent (report);
		report << "Output vector:  ";
		printVector<Field> (F, report, v3);

		for (j = 0; j < n; j++)
			if (!F.isZero (v3[j]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (x + a*y) - a*(y + a^-1*x) != 0" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDenseAXPY");

	return ret;
}

/* Test 4: Vector-vector axpy, sparse vectors
 *
 * Construct two random dense vectors x and y and a random element a and compute
 * (x + a*y) - a*(y + a^-1*x). Check whether the result is 0.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * iterations - Number of iterations over which to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testSparseAXPY (Field &F, long n, int iterations) 
{
	typedef vector <pair <size_t, typename Field::element> > Vector;

	commentator.start ("Testing sparse vector axpy", "testSparseAXPY", iterations);

	bool ret = true;
	bool iter_passed;

	Vector v1;
	Vector v2;
	Vector v3;
	Vector v4;
	typename Field::element a;
	typename Field::element ainv;
	typename Field::element aneg;
	typename Field::RandIter r (F);

	VectorDomain<Field, Vector, vector <typename Field::element> > VD (F);

	int i, j;
	Vector::iterator k;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		v1.clear ();
		v2.clear ();

		for (j = 0; j < n; j++) {
			// Give each sparse vector an entry about 10% of the time
			if (rand () % 100 < 10) {
				r.random (a);
				v1.push_back (pair <size_t, typename Field::element> (j, a));
			}

			if (rand () % 100 < 10) {
				r.random (a);
				v2.push_back (pair <size_t, typename Field::element> (j, a));
			}
		}

		do r.random (a); while (F.isZero (a));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		commentator.indent (report);
		report << "Element a:  ";
		F.write (report, a);
		report << endl;

		F.inv (ainv, a);
		F.neg (aneg, a);
		VD.axpy (v3, v1, a, v2);
		VD.axpy (v4, v2, ainv, v1);
		VD.axpyin (v3, aneg, v4);

		commentator.indent (report);
		report << "Output vector:  ";
		printVector<Field> (F, report, v3);

		for (k = v3.begin (); k < v3.end (); k++)
			if (!F.isZero ((*k).second))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: (x + a*y) - a*(y + a^-1*x) != 0" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSparseAXPY");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static long n = 100;
	static integer q = 101;
	static int iterations = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to N (default 100)",   TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",   TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Vector domain test suite" << endl << endl;
	cout.flush ();

	commentator.setBriefReportParameters (Commentator::OUTPUT_CONSOLE, false, false, false);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!testDenseDotProduct<Modular<long> >       (F, n, iterations)) pass = false;
	if (!testDenseSparseDotProduct<Modular<long> > (F, n, iterations)) pass = false;
	if (!testDenseAXPY<Modular<long> >             (F, n, iterations)) pass = false;
	if (!testSparseAXPY<Modular<long> >            (F, n, iterations)) pass = false;

	return pass ? 0 : -1;
}
