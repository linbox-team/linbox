/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-rank.C
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

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/rank.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Rank of diagonal matrix with n/2 distinct nonzero entries
 *
 * Construct a random diagonal matrix with n/2 distinct nonzero entries and see
 * that its computed determinant equals the product of diagonal entries
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalRank1 (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing diagonal rank (1)", "testDiagonalRank1", iterations);

	bool ret = true;
	bool done;
	int i;
	size_t j, k;

	Vector d(n);
	unsigned long _rank;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		for (j = 0; j < n / 2; j++) {
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
		}

		for (j = n / 2; j < n; j++)
			F.init (d[j], 0);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		printVector<Field> (F, report, d);

		commentator.indent (report);
		report << "True rank: " << n / 2 << endl;

		Blackbox D (F, d);

		rank <Field, Vector> (_rank, D, F);

		commentator.indent (report);
		report << "Computed rank: " << _rank << endl;

		if (_rank != n / 2) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed rank is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalRank1");

	return ret;
}

/* Test 2: Rank of diagonal matrix with n / 2 nondistinct nonzero entries
 *
 * Construct a random diagonal matrix with n / 2 nonzero, nondistinct entries
 * and see that its computed rank equals n / 2
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalRank2 (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing diagonal rank (2)", "testDiagonalRank2", iterations);

	bool ret = true;
	int i;
	size_t j, k;

	Vector d(n);
	unsigned long _rank;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		for (j = 0; j < n / 4; j++)
			do r.random (d[j]); while (F.isZero (d[j]));

		for (j = n / 4; j < n / 2; j++) {
			k = rand () % (n / 4);
			d[j] = d[k];
		}

		for (j = n / 2; j < n; j++)
			F.init (d[j], 0);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		printVector<Field> (F, report, d);

		commentator.indent (report);
		report << "True rank: " << n / 2 << endl;

		Blackbox D (F, d);

		rank <Field, Vector> (_rank, D, F);

		commentator.indent (report);
		report << "Computed rank: " << _rank << endl;

		if (_rank != n / 2) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed rank is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalRank2");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 256;
	static integer q = 2147483647U;
	static int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 256)",       TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Black box rank test suite" << endl << endl;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	if (!testDiagonalRank1<Modular<long> > (F, n, iterations)) pass = false;
	if (!testDiagonalRank2<Modular<long> > (F, n, iterations)) pass = false;

	return pass ? 0 : -1;
}
