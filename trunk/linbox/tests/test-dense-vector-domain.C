/* -*- mode: C++; style: linux -*- */

/* tests/test-vector-domain.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by Zhendong Wan<wan@mail.eecis.udel.edu>
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
#include "linbox/field/dense-vector-domain.h"
#include "linbox/util/vector-factory.h"
#include "linbox/util/field-axpy.h"

#include "test-common.h"

using namespace std;
using namespace LinBox;

/* Test 1: Dot product of dense vectors
 *
 * Construct two random dense vectors and compute their doc product
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make vectors
 * factory1 - Factory for first family of vectors
 * factory2 - Factory for second family of vectors
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDenseDotProduct (Field &F, long n,
				 VectorFactory<vector<typename Field::Element> > &factory1,
				 VectorFactory<vector<typename Field::Element> > &factory2) 
{
	typedef vector <typename Field::Element> Vector;

	commentator.start ("Testing dense/dense dot product", "testDenseDotProduct", factory1.m ());

	bool ret = true;

	Vector v1 (n), v2 (n);
	typename Field::Element sigma, rho;
	typename Field::RandIter r (F);

	DenseVectorDomain<Field> VD (F);

	int j;

	while (factory1 && factory2) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", factory1.j ());
		commentator.start (buf);

		F.init (sigma, 0);

		factory1.next (v1);
		factory2.next (v2);

		FieldAXPY<Field> r (F);

		for (j = 0; j < n; j++)
			r.accumulate (v1[j], v2[j]);

		r.get (sigma);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector 1:  ";
		printVector<Field> (F, report, v1);

		commentator.indent (report);
		report << "Input vector 2:  ";
		printVector<Field> (F, report, v2);

		VD.dotproduct (rho, v1, v2);

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
	typedef vector <typename Field::Element> Vector;

	commentator.start ("Testing dense vector axpy", "testDenseAXPY", iterations);

	bool ret = true;
	bool iter_passed;

	Vector v1 (n);
	Vector v2 (n);
	Vector v3 (n);
	Vector v4 (n);
	typename Field::Element a;
	typename Field::Element ainv;
	typename Field::Element aneg;
	typename Field::RandIter r (F);

	DenseVectorDomain<Field> VD (F);

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
		VD.axpy (v3, a, v1, v2);
		VD.axpy (v4, ainv,v2, v1);
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
	typedef vector <pair <size_t, typename Field::Element> > Vector;

	commentator.start ("Testing sparse vector axpy", "testSparseAXPY", iterations);

	bool ret = true;
	bool iter_passed;

	Vector v1;
	Vector v2;
	Vector v3;
	Vector v4;
	typename Field::Element a;
	typename Field::Element ainv;
	typename Field::Element aneg;
	typename Field::RandIter r (F);

	DenseVectorDomain<Field> VD (F);

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
				v1.push_back (pair <size_t, typename Field::Element> (j, a));
			}

			if (rand () % 100 < 10) {
				r.random (a);
				v2.push_back (pair <size_t, typename Field::Element> (j, a));
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

	RandomDenseVectorFactory<Modular<long> > factory1 (F, n, iterations), factory2 (F, n, iterations);

	if (!testDenseDotProduct<Modular<long> >       (F, n, factory1, factory2)) pass = false;

	if (!testDenseAXPY<Modular<long> >             (F, n, iterations)) pass = false;

	return pass ? 0 : -1;
}
