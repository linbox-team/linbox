/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-sum.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/vector/stream.h"
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/field/ntl-lzz_p.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/scalar-matrix.h"
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
static bool testZeroApply (Field &F, VectorStream<Vector> &stream1, VectorStream<Vector> &stream2) 
{
	commentator.start ("Testing zero apply", "testZeroApply", stream1.m ());

	bool ret = true;
	bool iter_passed = true;

	Vector d1, d2, v, w, zero;
	VectorDomain<Field> VD (F);
	typename Field::Element neg_one;

	VectorWrapper::ensureDim (zero, stream1.dim ());
	VectorWrapper::ensureDim (d1, stream1.dim ());
	VectorWrapper::ensureDim (d2, stream1.dim ());
	VectorWrapper::ensureDim (v, stream1.dim ());
	VectorWrapper::ensureDim (w, stream2.dim ());
// 	F.init (neg_one, 1);
// 	F.negin (neg_one);
	F.init (neg_one, -1);

	while (stream1) {
		commentator.startIteration (stream1.j ());
		iter_passed = true;

		stream1.next (d1);
		VD.mul (d2, d1, neg_one);

		Diagonal <Field> D1 (F, d1), D2 (F, d2);

		Sum <Diagonal<Field>,Diagonal <Field> > A (&D1, &D2);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal matrix:  ";
		VD.write (report, d1);
		report << endl;

		report << "Negative diagonal matrix:  ";
		VD.write (report, d2);
		report << endl;

		stream2.reset ();

		while (stream2) {
			stream2.next (w);

			report << "Input vector:  ";
			VD.write (report, w);
			report << endl;

			A.apply (v, w);

			report << "Output vector:  ";
			VD.write (report, v);
			report << endl;

			if (!VD.isZero (v))
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

	commentator.start ("Testing random transpose", "testRandomTranspose", iterations);

	Vector d(n);
	typename Field::RandIter r (F);

	for (int i = 0; i < n; i++)
		r.random (d[i]);

	Diagonal <Field, Vector> D (F, d);

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
	static integer q = 101;
	static int iterations1 = 10;
	static int iterations2 = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations1 },
		{ 'j', "-j J", "Apply test matrix to J vectors.", TYPE_INT,     &iterations2 },
		{ '\0' }
	};

//        typedef UnparametricField<NTL::zz_p> Field;
        typedef NTL_zz_p Field;
	NTL::zz_p::init(q);
	Field F(q);
	typedef vector<Field::Element> Vector;

	parseArguments (argc, argv, args);

	commentator.start("Sum black box test suite", "sum");

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	RandomDenseStream<Field> stream1 (F, n, iterations1), stream2 (F, n, iterations2);

	if (!testZeroApply (F, stream1, stream2)) pass = false;

	n = 10;
	RandomDenseStream<Field> stream3 (F, n, iterations1), stream4 (F, n, iterations2);

	Vector d1(n), d2(n);
	stream3.next (d1);
	stream4.next (d2);

//	Diagonal <Field, Vector> D1 (F, d1), D2 (F, d2);

	Field::Element d; F.init(d, 5);
	ScalarMatrix<Field> D1(F, 10, d), D2(F, 10, d); 
	typedef ScalarMatrix<Field> Blackbox;

	Sum <Blackbox, Blackbox> A (D1, D2);
	pass = pass && testBlackbox(F, A);

	Sum <Blackbox, Blackbox> Aref (&D1, &D2);
	pass = pass && testBlackbox(F, Aref);

	commentator.stop("Sum black box test suite");
	return pass ? 0 : -1;
}
