/* -*- mode: c; style: linux -*- */

/* tests/test-sum.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/util/vector-factory.h"
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/field/vector-domain.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/sum.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Application of zero matrix onto random vectors
 *
 * Construct a random diagonal matrix and its opposite, then construct
 * the sum of the two matrices. Apply to random vectors and check that
 * the result is zero.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testZeroApply (Field &F, VectorFactory<Vector> &factory1, VectorFactory<Vector> &factory2) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing zero apply", "testZeroApply", factory1.m ());

	bool ret = true;
	bool iter_passed = true;

	int j;

	Vector d1, d2, v, w, zero;
	VectorDomain<Field> VD (F);
	typename Field::Element neg_one;

	VectorWrapper::ensureDim (zero, factory1.n ());
	VectorWrapper::ensureDim (d2, factory1.n ());
	VectorWrapper::ensureDim (v, factory1.n ());
	F.init (neg_one, -1);

	while (factory1) {
		commentator.startIteration (factory1.j ());
		iter_passed = true;

		factory1.next (d1);
		VD.axpy (d2, zero, neg_one, d1);

		Blackbox D1 (F, d1), D2 (F, d2);
		Sum <Field, Vector> A (F, &D1, &D2);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal matrix:  ";
		printVector<Field> (F, report, d1);

		factory2.reset ();

		while (factory2) {
			factory2.next (w);

			commentator.indent (report);
			report << "Input vector:  ";
			printVector<Field> (F, report, w);

			A.apply (v, w);

			commentator.indent (report);
			report << "Output vector:  ";
			printVector<Field> (F, report, v);

			for (j = 0; j < factory1.n (); j++)
				if (!F.isZero (v[j]))
					ret = iter_passed = false;
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vector is not zero" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testZeroApply");

	return ret;
}

#if 0

/* Test 2: Random transpose
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
	typedef vector <typename Field::Element> Vector;
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

	bool ret = testTranspose<Field> (F, D, iterations);

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}

#endif

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations1 = 100;
	static int iterations2 = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations1 },
		{ 'j', "-j J", "Apply test matrix to J vectors (default 1)",                TYPE_INT,     &iterations2 },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Matrix sum black box test suite" << endl << endl;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	RandomDenseVectorFactory<Modular<long> > factory1 (F, n, iterations1), factory2 (F, n, iterations2);

	if (!testZeroApply<Modular<long> > (F, factory1, factory2)) pass = false;

	return pass ? 0 : -1;
}
