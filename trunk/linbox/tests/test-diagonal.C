/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-diagonal.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
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
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/util/vector-factory.h"

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

template <class Field, class Vector>
static bool testIdentityApply (Field &F, VectorFactory<Vector> &factory) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing identity apply", "testIdentityApply", factory.m ());

	bool ret = true;
	bool iter_passed = true;

	VectorDomain<Field> VD (F);
	Vector d;

	size_t i;

	VectorWrapper::ensureDim (d, factory.n ());

	for (i = 0; i < factory.n (); i++)
		F.init (VectorWrapper::ref<Field> (d, i), 1);

	Blackbox D (F, d);

	Vector v, w;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, factory.n ());

	while (factory) {
		commentator.startIteration (i);

		iter_passed = true;

		factory.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		D.apply (w, v);

		commentator.indent (report);
		report << "Output vector: ";
		VD.write (report, w);
		report << endl;

		if (!VD.areEqual (w, v))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	factory.reset ();

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

template <class Field, class Vector>
static bool testRandomMinpoly (Field &F, VectorFactory<Vector> &factory) 
{
	typedef vector <typename Field::Element> Polynomial;
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing random minpoly", "testRandomMinpoly", factory.m ());

	bool ret = true;

	size_t j;
	typename Field::Element pi;
	Polynomial m_D;
	VectorDomain<Field> VD (F);

	Vector d;

	VectorWrapper::ensureDim (d, factory.n ());

	while (factory) {
		commentator.startIteration (factory.j ());

		F.init (pi, 1);

		factory.next (d);

		for (j = 0; j < factory.n (); j++)
			F.mulin (pi, VectorWrapper::constRef<Field> (d, j));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal vector: ";
		VD.write (report, d);
		report << endl;

		commentator.indent (report);
		report << "Product: ";
		F.write (report, pi);
		report << endl;

		Blackbox D (F, d);
		minpoly (m_D, D, F);

		commentator.indent (report);
		report << "Minimal polynomial: ";
		printPolynomial (F, report, m_D);

		if (!F.areEqual (m_D[0], pi)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: m_D(0) != det(D)" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	factory.reset ();

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

template <class Field, class Vector>
static bool testRandomTranspose (Field &F,
				 VectorFactory<Vector> &factory1,
				 VectorFactory<Vector> &factory2) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing random transpose", "testRandomTranspose", factory1.m ());

	Vector d;
	typename Field::RandIter r (F);

	VectorWrapper::ensureDim (d, factory1.n ());

	for (size_t i = 0; i < factory1.n (); i++)
		r.random (VectorWrapper::ref<Field> (d, i));

	Blackbox D (F, d);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	report << "Diagonal vector: ";
	printVector (F, report, d);

	bool ret = testTranspose (F, D, factory1, factory2);

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

	typedef Modular<long> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	srand (time (NULL));

	cout << "Diagonal matrix black box test suite" << endl << endl;

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (2);

	RandomDenseVectorFactory<Field> factory1 (F, n, iterations), factory2 (F, n, iterations);
	RandomDenseVectorFactory<Field, NonzeroRandIter<Field> >
		factory3 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, iterations);

	if (!testIdentityApply    (F, factory1)) pass = false;
	if (!testRandomMinpoly    (F, factory3)) pass = false;
	if (!testRandomTranspose  (F, factory1, factory2)) pass = false;

	return pass ? 0 : -1;
}
