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

/* Test 2: Solution of diagonal system
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
static bool testNonsingularSolve (const Field &F,
				  VectorStream<Vector> &stream1, 
				  VectorStream<Vector> &stream2) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing nonsingular solve", "testNonsingularSolve", stream1.m ());

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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNonsingularSolve");

	return ret;
}

/* Test 3: Solution of a consistent singular diagonal system with nonsingular
 * leading principal minor
 *
 * Constructs a random diagonal matrix D of rank r and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testSingularConsistentSolve (const Field &F,
					 unsigned int n,
					 VectorStream<Vector> &stream1,
					 VectorStream<Vector> &stream2) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing singular consistent solve", "testSingularConsistentSolve", stream1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;
	bool iter_passed;

	Vector d1, b1, d, b, x, y;

	VectorWrapper::ensureDim (d, n);
	VectorWrapper::ensureDim (b, n);
	VectorWrapper::ensureDim (x, n);
	VectorWrapper::ensureDim (y, n);

	SolverTraits traits (SolverTraits::METHOD_WIEDEMANN, false);

	while (stream1 && stream2) {
		commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (d1);
		stream2.next (b1);

		VD.copy (d, d1);
		VD.copy (b, b1);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		commentator.indent (report);
		report << "Right-hand side:  ";
		VD.write (report, b);
		report << endl;

		Blackbox D (F, d);

		solve (D, x, b, F, traits);

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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSingularConsistentSolve");

	return ret;
}

/* Test 4: Solution of an inconsistent singular diagonal system with nonsingular
 * leading principal minor
 *
 * Constructs a random diagonal matrix D of rank r and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testSingularInconsistentSolve (const Field &F,
					   VectorStream<Vector> &stream1,
					   VectorStream<Vector> &stream2) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing singular inconsistent solve", "testSingularInconsistentSolve", stream1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;
	bool iter_passed;

	Vector d1, d, b, x, y;

	VectorWrapper::ensureDim (d, stream2.dim ());
	VectorWrapper::ensureDim (b, stream2.dim ());
	VectorWrapper::ensureDim (x, stream2.dim ());
	VectorWrapper::ensureDim (y, stream2.dim ());

	while (stream1 && stream2) {
		commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (d1);
		stream2.next (b);

		VD.copy (d, d1);

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

		commentator.stop ("done");
		commentator.progress ();
	}

	stream1.reset ();
	stream2.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSingularInconsistentSolve");

	return ret;
}

/* Test 5: Solution of an inconsistent singular diagonal system with singular
 * leading principal minor
 *
 * Constructs a random diagonal matrix D of rank r and a random right-hand
 * side b, and computes the solution to the Dx=b, checking the result
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector, class SparseVector>
static bool testSingularPreconditionedSolve (const Field &F,
					     VectorStream<SparseVector> &stream1,
					     VectorStream<Vector> &stream2) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing singular preconditioned solve", "testSingularPreconditionedSolve", stream1.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;
	bool iter_passed;

	SparseVector d1;
	Vector d, b, x, y;

	VectorWrapper::ensureDim (d, stream2.dim ());
	VectorWrapper::ensureDim (b, stream2.dim ());
	VectorWrapper::ensureDim (x, stream2.dim ());
	VectorWrapper::ensureDim (y, stream2.dim ());

	while (stream1 && stream2) {
		commentator.startIteration (stream1.j ());

		iter_passed = true;

		stream1.next (d1);
		stream2.next (b);

		VD.copy (d, d1);

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

		commentator.stop ("done");
		commentator.progress ();
	}

	stream1.reset ();
	stream2.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testSingularPreconditionedSolve");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 100;
	static size_t r = 20;
	static integer q = 2147483647U;
	static int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 100)",       TYPE_INT,     &n },
		{ 'r', "-r R", "Set singular system rank to R (default 20)",                TYPE_INT,     &r },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
	};

	typedef Modular<uint32> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	cout << "Solve test suite" << endl << endl;

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	RandomDenseStream<Field> stream1 (F, n, iterations), stream2 (F, n, iterations);
	RandomDenseStream<Field> stream3 (F, r, iterations), stream4 (F, r, iterations);
	RandomSparseStream<Field> stream6 (F, n, (double) r / (double) n, iterations);

	if (!testIdentitySolve               (F, stream1)) pass = false;
	if (!testNonsingularSolve            (F, stream1, stream2)) pass = false;
	if (!testSingularConsistentSolve     (F, n, stream3, stream4)) pass = false;
	if (!testSingularInconsistentSolve   (F, stream3, stream2)) pass = false;
	if (!testSingularPreconditionedSolve (F, stream6, stream2)) pass = false;

	return pass ? 0 : -1;
}
