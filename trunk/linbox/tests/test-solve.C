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
#include "linbox/vector/stream.h"
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
static bool testIdentitySolve (const Field &F, VectorStream<Vector> &stream) 
{
	typedef ScalarMatrix <Field, Vector> Blackbox;

	commentator.start ("Testing identity solve", "testIdentitySolve", stream.m ());

	bool ret = true;
	bool iter_passed = true;

	VectorDomain<Field> VD (F);

	typename Field::Element s;
	F.init (s, 1);
	Blackbox I (F, stream.n (), s);

	Vector v, w;

	VectorWrapper::ensureDim (v, stream.n ());
	VectorWrapper::ensureDim (w, stream.n ());

	while (stream) {
		commentator.startIteration (stream.j ());

		iter_passed = true;

		stream.next (v);

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

	stream.reset ();

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
				    VectorStream<Vector> &x_stream,
				    VectorStream<Vector> &v_stream) 
{
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing Vandermonde inverse", "testVandermondeInverse", x_stream.m ());

	bool ret = true;
	bool inner_iter_passed;

	VectorDomain<Field> VD (F);
	size_t j, k;

	Blackbox V (F, x_stream.n (), x_stream.n ());

	Vector x, v, w, z;
	typename Field::Element t;

	VectorWrapper::ensureDim (x, x_stream.n ());
	VectorWrapper::ensureDim (v, x_stream.n ());
	VectorWrapper::ensureDim (w, x_stream.n ());
	VectorWrapper::ensureDim (z, x_stream.n ());

	while (x_stream) {
		commentator.startIteration (x_stream.j ());

		/* Evaluation points */
		x_stream.next (x);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Evaluation points: ";
		VD.write (report, x);
		report << endl;
		report.flush ();

		/* Build the Vandermonde matrix */
		for (j = 0; j < x_stream.n (); j++) {
			F.init (t, 1);

			for (k = 0; k < x_stream.n (); k++) {
				V.setEntry (j, k, t);
				F.mulin (t, VectorWrapper::ref<Field> (x, j));
			}
		}

		Inverse<Field, Vector> VT (F, &V);

		v_stream.reset ();

		while (v_stream) {
			inner_iter_passed = true;

			/* Random vector of evaluation results */
			v_stream.next (v);

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

	x_stream.reset ();

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
			       VectorStream<Vector> &stream1, 
			       VectorStream<Vector> &stream2) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing diagonal solve", "testDiagonalSolve", stream1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;
	bool iter_passed;

	Vector d, b, x, y;

	VectorWrapper::ensureDim (d, stream1.n ());
	VectorWrapper::ensureDim (b, stream1.n ());
	VectorWrapper::ensureDim (x, stream1.n ());
	VectorWrapper::ensureDim (y, stream1.n ());

	while (stream1 && stream2) {
		commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (d);
		stream2.next (b);

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

		commentator.indent (report);
		report << "Output:           ";
		VD.write (report, y);
		report << endl;

		if (!VD.areEqual (y, b))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed solution is incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	stream1.reset ();
	stream2.reset ();

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

	typedef Modular<uint32> Field;
	typedef vector<Field::Element> DenseVector;

	parseArguments (argc, argv, args);
	Field F (q);

	srand (time (NULL));

	cout << "Solve test suite" << endl << endl;

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	RandomDenseStream<Field, DenseVector> stream1 (F, n, iterations), stream2 (F, n, iterations);

	if (!testIdentitySolve (F, stream1)) pass = false;
	if (!testDiagonalSolve (F, stream1, stream2)) pass = false;

	return pass ? 0 : -1;
}
