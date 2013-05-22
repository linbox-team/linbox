
/* tests/test-inverse.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */


/*! @file  tests/test-inverse.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/hilbert.h"
#include "linbox/blackbox/inverse.h"
#include "linbox/vector/stream.h"

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
static bool testIdentityInverse (const Field &F, VectorStream<Vector> &stream)
{
	typedef Diagonal<Field> Blackbox;

	commentator().start ("Testing identity inverse", "testIdentityInverse", stream.m ());

	bool ret = true;
	bool iter_passed = true;

	Vector d(F);
	VectorDomain<Field> VD (F);

	size_t i;

	VectorWrapper::ensureDim (d, stream.n ());

	for (i = 0; i < stream.n (); i++)
		F.init (VectorWrapper::ref<Field> (d, i), 1);

	Blackbox D (d);
	Inverse<Blackbox> DT (&D);

	Vector v(F), w(F);

	VectorWrapper::ensureDim (v, stream.n ());
	VectorWrapper::ensureDim (w, stream.n ());

	while (stream) {
		commentator().startIteration ((unsigned)stream.j ());

		iter_passed = true;

		stream.next (v);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		DT.apply (w, v);

		report << "Output vector: ";
		VD.write (report, w);
		report << endl;

		if (!VD.areEqual (w, v))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator().stop ("done");
		commentator().progress ();
	}

	stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testIdentityInverse");

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
static bool testHilbertInverse (const Field &F, VectorStream<Vector> &stream)
{
	typedef Hilbert <Field> Blackbox;

	commentator().start ("Testing Hilbert inverse", "testHilbertInverse", stream.m ());

	bool ret = true;
	bool iter_passed;

	VectorDomain<Field> VD (F);

	Blackbox H (F, stream.n ());
	Inverse<Blackbox> HT (&H);

	Vector v(F), w(F), z(F);

	VectorWrapper::ensureDim (v, stream.n ());
	VectorWrapper::ensureDim (w, stream.n ());
	VectorWrapper::ensureDim (z, stream.n ());

	while (stream) {
		commentator().startIteration (stream.j ());

		iter_passed = true;

		stream.next (v);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector: ";
		VD.write (report, v);
		report << endl;

		H.apply (z, v);
		HT.apply (w, z);

		report << "Output vector: ";
		VD.write (report, w);
		report << endl;

		if (!VD.areEqual (w, v))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator().stop ("done");
		commentator().progress ();
	}

	stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testHilbertInverse");

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
				    VectorStream<Vector> &x_stream,
				    VectorStream<Vector> &v_stream)
{
	typedef BlasMatrix <Field> Blackbox;

	commentator().start ("Testing Vandermonde inverse", "testVandermondeInverse", x_stream.m ());

	bool ret = true;
	bool inner_iter_passed;

	VectorDomain<Field> VD (F);
	size_t j, k;

	Blackbox V (F, x_stream.n (), x_stream.n ());

	Vector x(F), v(F), w(F), z(F);
	typename Field::Element t;

	VectorWrapper::ensureDim (x, x_stream.n ());
	VectorWrapper::ensureDim (v, x_stream.n ());
	VectorWrapper::ensureDim (w, x_stream.n ());
	VectorWrapper::ensureDim (z, x_stream.n ());

	while (x_stream) {
		commentator().startIteration ((unsigned)x_stream.j ());

		/* Evaluation points */
		x_stream.next (x);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
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

		Inverse<Blackbox> VT (&V);

		v_stream.reset ();

		while (v_stream) {
			inner_iter_passed = true;

			/* Random vector of evaluation results */
			v_stream.next (v);

			report << "Input vector: ";
			VD.write (report, v);
			report << endl;

			/* w should now be the requisite polynomial */
			VT.apply (w, v);

			report << "Output vector: ";
			VD.write (report, w);
			report << endl;

			/* Multipoint evaluation to check whether w is correct */
			multiEvalPoly (F, z, w, x);

			report << "Evaluation results: ";
			VD.write (report, z);
			report << endl;

			if (!VD.areEqual (z, v))
				ret = inner_iter_passed = false;

			if (!inner_iter_passed)
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Vectors are not equal" << endl;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	x_stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testVandermondeInverse");

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
static bool testDiagonalInverse (const Field &F, VectorStream<Vector> &stream)
{
	typedef Diagonal <Field> Blackbox;

	commentator().start ("Testing diagonal inverse", "testDiagonalInverse", stream.m ());

	VectorDomain<Field> VD (F);

	bool ret = true;
	bool iter_passed;

	size_t j;

	Vector d(F), di(F), dt(F), e(F), DTe(F);

	VectorWrapper::ensureDim (d, stream.n ());
	VectorWrapper::ensureDim (di, stream.n ());
	VectorWrapper::ensureDim (dt, stream.n ());
	VectorWrapper::ensureDim (e, stream.n ());
	VectorWrapper::ensureDim (DTe, stream.n ());

	while (stream) {
		commentator().startIteration ((unsigned)stream.j ());

		iter_passed = true;

		stream.next (d);

		for (j = 0; j < stream.n (); j++)
			F.inv (VectorWrapper::ref<Field> (di, j), VectorWrapper::ref<Field> (d, j));

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		VD.write (report, d);
		report << endl;

		report << "Expeted diagonal entries of inverse: ";
		VD.write (report, di);
		report << endl;

		Blackbox D (d);
		Inverse <Blackbox> DT (&D);

		for (j = 0; j < stream.n (); j++) {
			F.init (VectorWrapper::ref<Field> (e, j), 1);
			DT.apply (DTe, e);
		}

		VD.copy (dt, DTe);

		report << "Diagonal entries of computed inverse: ";
		VD.write (report, dt);
		report << endl;

		if (!VD.areEqual (di, dt))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed inverse does not match expected inverse" << endl;

		commentator().stop ("done");
		commentator().progress ();
	}

	stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalInverse");

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
				 VectorStream<Vector> &stream1,
				 VectorStream<Vector> &stream2)
{
	typedef BlasMatrix <Field> Blackbox;

	commentator().start ("Testing random transpose", "testRandomTranspose", stream1.m ());

	size_t i, j;
	typename Field::Element x;
	typename Field::RandIter r (F);

	Blackbox A (F, stream1.n (), stream2.n ());
	Inverse<Blackbox> Ainv (&A);

	for (i = 0; i < stream1.n (); i++) {
		for (j = 0; j < stream2.n (); j++) {
			r.random (x);
			A.setEntry (i, j, x);
		}
	}

	bool ret = testTranspose (F, Ainv, stream1, stream2);

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}


int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static unsigned int iterations = 100;
	static int N = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 'N', "-N N", "Apply Vandermonde inverse to N vectors.", TYPE_INT,     &N },
		END_OF_ARGUMENTS
	};

	typedef Modular<uint32_t> Field; //C.Pernet: avoids confusion with givaro::uint32_t
	typedef BlasVector<Field> Vector;

	parseArguments (argc, argv, args);
	Field F (q);

	srand ((unsigned)time (NULL));

	commentator().start("Inverse black box test suite", "Inverse");

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	RandomDenseStream<Field, Vector> stream1 (F, n, iterations), stream2 (F, n, iterations);
	RandomDenseStream<Field, Vector> stream3 (F, n, (size_t)N);

	if (!testIdentityInverse    (F, stream1)) pass = false;
	if (!testVandermondeInverse (F, stream1, stream3)) pass = false;
	if (!testDiagonalInverse    (F, stream1)) pass = false;
	if (!testRandomTranspose    (F, stream1, stream2)) pass = false;

	commentator().stop("Inverse black box test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
