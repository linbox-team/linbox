/* -*- mode: c; style: linux -*- */

/* tests/test-trace.C
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

#include "test-common.h"

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/trace.h"

using namespace LinBox;

/* Test 1: Trace of random diagonal matrix
 *
 * Construct a random diagonal matrix and check that its computed trace is the
 * same as the sum of its entries
 *
 * F - Field over which to perform computations
 * factory - Factory that comprises source of diagonal vectors
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalTrace (Field &F, VectorFactory<vector<typename Field::Element> > &factory) 
{
	typedef vector <typename Field::Element> Vector;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing diagonal trace", "testDiagonalTrace", factory.m ());

	bool ret = true;
	bool done;
	int i, j, k;

	Vector d;
	typename Field::Element sigma, res;

	while (factory) {
		commentator.startIteration (factory.j ());

		factory.next (d);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		printVector<Field> (F, report, d);

		F.init (sigma, 0);
		for (i = 0; i < factory.n (); i++)
			F.addin (sigma, VectorWrapper::constRef<Field, Vector> (d, i));

		commentator.indent (report);
		report << "True trace: ";
		F.write (report, sigma);
		report << endl;

		Blackbox D (F, d);

		trace <Field, Vector> (res, D, F);

		commentator.indent (report);
		report << "Computed trace: ";
		F.write (report, res);
		report << endl;

		if (!F.areEqual (sigma, res)) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed trace is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalTrace");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 256;
	static integer q = 101;
	static int iterations = 10;
	static int numVectors = 100;
	static int k = 3;
	static int N = 20;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 256)", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)",  TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",     TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Black box trace test suite" << endl << endl;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	RandomDenseVectorFactory<Modular<long> > factory (F, n, iterations);

	if (!testDiagonalTrace<Modular<long> > (F, factory)) pass = false;

	return pass ? 0 : -1;
}
