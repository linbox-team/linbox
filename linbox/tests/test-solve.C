/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-solve.C
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
#include "linbox/field/modular.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/util/vector-factory.h"
#include "linbox/solutions/solve.h"

#include "linbox/solutions/minpoly.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Solve Ix=b, check that x == b
 *
 * Constructs a black box for the inverse of an n x n identity matrix and checks
 * that that inverse is itself the identity operator by applying it to random
 * vectors.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply identity inverse
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testIdentitySolve (const Field &F, VectorFactory<Vector> &factory) 
{
	typedef ScalarMatrix <Field, Vector> Blackbox;

	commentator.start ("Testing identity solve", "testIdentitySolve", factory.m ());

	bool ret = true;
	bool iter_passed = true;

	VectorDomain<Field> VD (F);

	typename Field::Element s;
	F.init (s, 1);
	Blackbox I (F, factory.n (), s);

	Vector v, w;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, factory.n ());

	while (factory) {
		commentator.startIteration (factory.j ());

		iter_passed = true;

		factory.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		solve (I, w, v, F);

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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentitySolve");

	return ret;
}

#if 0

/* Test 3: Inverse of Vandermonde matrix
 *
 * Computes a random Vandermonde matrix and its inverse. This inverse is a
 * linear operator that interpolates the values given in the input vector to
 * produce a polynomial whose coefficients are the elements of the output
 * vector. We then evaluate the polynomial in Horner fashion at each of the
 * evaluation points generated above to check whether the result is the original
 * input vector.
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 * N - Number of random vectors to which to apply random Vandermonde matrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testVandermondeInverse (const Field           &F,
				    VectorFactory<Vector> &x_factory,
				    VectorFactory<Vector> &v_factory) 
{
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing Vandermonde inverse", "testVandermondeInverse", x_factory.m ());

	bool ret = true;
	bool inner_iter_passed;

	VectorDomain<Field> VD (F);
	size_t j, k;

	Blackbox V (F, x_factory.n (), x_factory.n ());

	Vector x, v, w, z;
	typename Field::Element t;

	VectorWrapper::ensureDim (x, x_factory.n ());
	VectorWrapper::ensureDim (v, x_factory.n ());
	VectorWrapper::ensureDim (w, x_factory.n ());
	VectorWrapper::ensureDim (z, x_factory.n ());

	while (x_factory) {
		commentator.startIteration (x_factory.j ());

		/* Evaluation points */
		x_factory.next (x);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Evaluation points: ";
		VD.write (report, x);
		report << endl;
		report.flush ();

		/* Build the Vandermonde matrix */
		for (j = 0; j < x_factory.n (); j++) {
			F.init (t, 1);

			for (k = 0; k < x_factory.n (); k++) {
				V.setEntry (j, k, t);
				F.mulin (t, VectorWrapper::ref<Field> (x, j));
			}
		}

		Inverse<Field, Vector> VT (F, &V);

		v_factory.reset ();

		while (v_factory) {
			inner_iter_passed = true;

			/* Random vector of evaluation results */
			v_factory.next (v);

			commentator.indent (report);
			report << "Input vector: ";
			VD.write (report, v);
			report << endl;

			/* w should now be the requisite polynomial */
			VT.apply (w, v);

			commentator.indent (report);
			report << "Output vector: ";
			VD.write (report, w);
			report << endl;

			/* Multipoint evaluation to check whether w is correct */
			multiEvalPoly (F, z, w, x);

			commentator.indent (report);
			report << "Evaluation results: ";
			VD.write (report, z);
			report << endl;

			if (!VD.areEqual (z, v))
				ret = inner_iter_passed = false;

			if (!inner_iter_passed)
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Vectors are not equal" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	x_factory.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testVandermondeInverse");

	return ret;
}

#endif

/* Test 3: Solution of diagonal system
 *
 * Constructs a random nonsingular diagonal matrix D and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testDiagonalSolve (const Field &F,
			       VectorFactory<Vector> &factory1, 
			       VectorFactory<Vector> &factory2) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing diagonal solve", "testDiagonalSolve", factory1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;
	bool iter_passed;

	Vector d, b, x, y;

	VectorWrapper::ensureDim (d, factory1.n ());
	VectorWrapper::ensureDim (b, factory1.n ());
	VectorWrapper::ensureDim (x, factory1.n ());
	VectorWrapper::ensureDim (y, factory1.n ());

	while (factory1 && factory2) {
		commentator.startIteration (factory1.j ());

		iter_passed = true;

		factory1.next (d);
		factory2.next (b);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		commentator.indent (report);
		report << "Right-hand side:  ";
		VD.write (report, b);
		report << endl;

		Blackbox D (F, d);

		solve (D, x, b, F);

		commentator.indent (report);
		report << "System solution:  ";
		VD.write (report, x);
		report << endl;

		D.apply (y, x);

		if (!VD.areEqual (y, b))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed inverse does not match expected inverse" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	factory1.reset ();
	factory2.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalSolve");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 100;
	static int N = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations },
		{ 'N', "-N N", "Apply Vandermonde inverse to N vectors (default 1)",        TYPE_INT,     &N },
	};

	typedef Modular<long> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	srand (time (NULL));

	cout << "Solve test suite" << endl << endl;

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);

	RandomDenseVectorFactory<Field > factory1 (F, n, iterations), factory2 (F, n, iterations);

	if (!testIdentitySolve (F, factory1)) pass = false;
	if (!testDiagonalSolve (F, factory1, factory2)) pass = false;

	return pass ? 0 : -1;
}
