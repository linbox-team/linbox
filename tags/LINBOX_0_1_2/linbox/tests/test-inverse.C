/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-inverse.C
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
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/hilbert.h"
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/inverse.h"
#include "linbox/util/vector-factory.h"

#include "linbox/solutions/minpoly.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Inverse of the identity matrix
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
static bool testIdentityInverse (const Field &F, VectorFactory<Vector> &factory) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing identity inverse", "testIdentityInverse", factory.m ());

	bool ret = true;
	bool iter_passed = true;

	Vector d;
	VectorDomain<Field> VD (F);

	size_t i;

	VectorWrapper::ensureDim (d, factory.n ());

	for (i = 0; i < factory.n (); i++)
		F.init (VectorWrapper::ref<Field> (d, i), 1);

	Blackbox D (F, d);
	Inverse<Field, Vector> DT (F, &D);

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

		DT.apply (w, v);

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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityInverse");

	return ret;
}

/* Test 2: Inverse of Hilbert matrix
 *
 * Constructs an n x n Hilbert matrix and a black box for its inverse. Applies
 * each to random vectors and checks that the results are equal.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testHilbertInverse (const Field &F, VectorFactory<Vector> &factory) 
{
	typedef Hilbert <Field, Vector> Blackbox;

	commentator.start ("Testing Hilbert inverse", "testHilbertInverse", factory.m ());

	bool ret = true;
	bool iter_passed;

	VectorDomain<Field> VD (F);

	Blackbox H (F, factory.n ());
	Inverse<Field, Vector> HT (F, &H);

	Vector v, w, z;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, factory.n ());
	VectorWrapper::ensureDim (z, factory.n ());

	while (factory) {
		commentator.startIteration (factory.j ());

		iter_passed = true;

		factory.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector: ";
		VD.write (report, v);
		report << endl;

		H.apply (z, v);
		HT.apply (w, z);

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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testHilbertInverse");

	return ret;
}

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

/* Test 3: Inverse of diagonal inverse
 *
 * Constructs a random nonsingular diagonal matrix and its inverse, and extracts
 * the values along the diagonal of the inverse. Checks that those values are in
 * fact the inverses of the diagonal elements in the original.
 * 
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testDiagonalInverse (const Field &F, VectorFactory<Vector> &factory) 
{
	typedef Diagonal <Field, Vector> Blackbox;

	commentator.start ("Testing diagonal inverse", "testDiagonalInverse", factory.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;
	bool iter_passed;

	size_t j;

	Vector d, di, dt, e, DTe;

	VectorWrapper::ensureDim (d, factory.n ());
	VectorWrapper::ensureDim (di, factory.n ());
	VectorWrapper::ensureDim (dt, factory.n ());
	VectorWrapper::ensureDim (e, factory.n ());
	VectorWrapper::ensureDim (DTe, factory.n ());

	while (factory) {
		commentator.startIteration (factory.j ());

		iter_passed = true;

		factory.next (d);

		for (j = 0; j < factory.n (); j++)
			F.inv (VectorWrapper::ref<Field> (di, j), VectorWrapper::ref<Field> (d, j));

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		commentator.indent (report);
		report << "Expeted diagonal entries of inverse: ";
		VD.write (report, di);
		report << endl;

		Blackbox D (F, d);
		Inverse <Field, Vector> DT (F, &D);

		for (j = 0; j < factory.n (); j++) {
			F.init (VectorWrapper::ref<Field> (e, j), 1);
			DT.apply (DTe, e);
		}

		VD.copy (dt, DTe);

		commentator.indent (report);
		report << "Diagonal entries of computed inverse: ";
		VD.write (report, dt);
		report << endl;

		if (!VD.areEqual (di, dt))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed inverse does not match expected inverse" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	factory.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalInverse");

	return ret;
}

/* Test 3: Random transpose
 *
 * Compute the inverse of a random dense matrix and apply its transpose to
 * random vectors
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
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing random transpose", "testRandomTranspose", factory1.m ());

	size_t i, j;
	typename Field::Element x;
	typename Field::RandIter r (F);

	Blackbox A (F, factory1.n (), factory2.n ());
	Inverse<Field, Vector> Ainv (F, &A);

	for (i = 0; i < factory1.n (); i++) {
		for (j = 0; j < factory2.n (); j++) {
			r.random (x);
			A.setEntry (i, j, x);
		}
	}

	bool ret = testTranspose (F, Ainv, factory1, factory2);

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

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

	parseArguments (argc, argv, args);
	Field F (q);

	srand (time (NULL));

	cout << "Black box inverse test suite" << endl << endl;

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	RandomDenseVectorFactory<Field > factory1 (F, n, iterations), factory2 (F, n, iterations);
	RandomDenseVectorFactory<Field > factory3 (F, n, N);

	if (!testIdentityInverse    (F, factory1)) pass = false;
	if (!testVandermondeInverse (F, factory1, factory3)) pass = false;
	if (!testDiagonalInverse    (F, factory1)) pass = false;
	if (!testRandomTranspose    (F, factory1, factory2)) pass = false;

	return pass ? 0 : -1;
}
