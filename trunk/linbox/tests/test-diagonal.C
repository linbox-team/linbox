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
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
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
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityApply (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing identity apply", "testIdentityApply", iterations);

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
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		iter_passed = true;

		for (j = 0; j < n; j++)
			r.random (v[j]);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		printVector<Field> (F, report, v);

		D.apply (w, v);

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

/* Test 2: Constant term in minimal polynomial of diagonal map
 *
 * Generates a random diagonal nonsingular matrix and computes its minimal
 * polynomial. Checks that the constant term thereof equals the product of the
 * entries on the diagonal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomMinpoly (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef vector <typename Field::element> Polynomial;
	typedef vector <pair <size_t, typename Field::element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing random minpoly", "testRandomMinpoly", iterations);

	bool ret = true;

	int i, j;
	typename Field::element pi;
	Polynomial m_D;

	Vector d(n);
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator.start (buf);

		F.init (pi, 1);
		for (j = 0; j < n; j++) {
			do r.random (d[j]); while (F.isZero (d[j]));
			F.mulin (pi, d[j]);
		}

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal vector: ";
		printVector<Field> (F, report, d);

		commentator.indent (report);
		report << "Product: ";
		F.write (report, pi);
		report << endl;

		Blackbox D (F, d);
		minpoly<Field, Polynomial, Vector> (m_D, D, F);

		commentator.indent (report);
		report << "Minimal polynomial: ";
		printPolynomial<Field, Polynomial> (F, report, m_D);

		if (!F.areEqual (m_D[0], pi)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: m_D(0) != det(D)" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomMinpoly");

	return ret;
}

/* Test 3: Random transpose
 *
 * Compute a random diagonal matrix and use the transpose test in test-generic.h
 * to check consistency of transpose apply.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomTranspose (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::element> Vector;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing random transpose", "testRandomTranspose", iterations);

	Vector d(n);
	typename Field::RandIter r (F);

	for (int i = 0; i < n; i++)
		r.random (d[i]);

	Blackbox D (F, d);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	report << "Diagonal vector: ";
	printVector<Field> (F, report, d);

	bool ret = testTranpose<Field> (F, D, iterations);

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Diagonal matrix black box test suite" << endl << endl;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	if (!testIdentityApply<Modular<long> >    (F, n, iterations)) pass = false;
	if (!testRandomMinpoly<Modular<long> >    (F, n, iterations)) pass = false;
	if (!testRandomTranspose<Modular<long> >  (F, n, iterations)) pass = false;

	return pass ? 0 : -1;
}
