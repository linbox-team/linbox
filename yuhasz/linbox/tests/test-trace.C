/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-trace.C
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

#include "test-common.h"

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/trace.h"
#include "linbox/vector/stream.h"

using namespace LinBox;

/* Test 1: Trace of random diagonal matrix
 *
 * Construct a random diagonal matrix and check that its computed trace is the
 * same as the sum of its entries
 *
 * F - Field over which to perform computations
 * stream - Stream that comprises source of diagonal vectors
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalTrace (const Field &F, VectorStream<vector<typename Field::Element> > &stream) 
{
	typedef vector <typename Field::Element> Vector;
	typedef Diagonal <Field> Blackbox;

	commentator.start ("Testing diagonal trace", "testDiagonalTrace", stream.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;
	size_t i;

	Vector d;
	typename Field::Element sigma, res;

	VectorWrapper::ensureDim (d, stream.dim ());

	while (stream) {
		commentator.startIteration (stream.j ());

		stream.next (d);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, d);
		report << endl;

		F.init (sigma, 0);
		for (i = 0; i < stream.n (); i++)
			F.addin (sigma, VectorWrapper::constRef<Field, Vector> (d, i));

		report << "True trace: ";
		F.write (report, sigma);
		report << endl;

		Blackbox D (F, d);

		trace<Vector,Field,Blackbox> (res, D, F);

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

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 256)", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)",  TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",     TYPE_INT,     &iterations },
	};

	typedef Modular<uint32> Field;
	typedef vector<Field::Element> Vector;

	parseArguments (argc, argv, args);
	Field F (q);

	cout << endl << "Black box trace test suite" << endl;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	RandomDenseStream<Field, Vector> stream (F, n, iterations);

	if (!testDiagonalTrace (F, stream)) pass = false;

	return pass ? 0 : -1;
}
