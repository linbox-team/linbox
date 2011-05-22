
/* tests/test-diagonal.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Time-stamp: <22 Jun 10 15:59:39 Jean-Guillaume.Dumas@imag.fr>
 * --------------------------------------------------------
 *
 * See COPYING for license information
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/blackbox/diagonal.h"
#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/solutions/rank.h"
#include "linbox/vector/stream.h"

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
static bool testIdentityApply (Field &F, VectorStream<Vector> &stream) 
{
        typedef LinBox::Diagonal<Field> Blackbox;

	commentator.start ("Testing identity apply", "testIdentityApply", stream.m ());

	bool ret = true;
	bool iter_passed = true;

	VectorDomain<Field> VD (F);
	Vector d;

	size_t i;

	VectorWrapper::ensureDim (d, stream.n ());

	for (i = 0; i < stream.n (); i++)
		F.init (VectorWrapper::ref<Field> (d, i), 1);

	Blackbox D (F, d);

	Vector v, w;

	VectorWrapper::ensureDim (v, stream.n ());
	VectorWrapper::ensureDim (w, stream.n ());

	while (stream) {
		commentator.startIteration (i);

		iter_passed = true;

		stream.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector:  ";
		VD.write (report, v);
		report << endl;

		D.apply (w, v);

		report << "Output vector: ";
		VD.write (report, w);
		report << endl;

		if (!VD.areEqual (w, v))
			ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop (MSG_STATUS (ret));
		commentator.progress ();
	}

	stream.reset ();

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
static bool testRandomMinpoly (Field &F, VectorStream<Vector> &stream) 
{
	typedef vector <typename Field::Element> Polynomial;
	typedef LinBox::Diagonal <Field> Blackbox;

	commentator.start ("Testing random minpoly", "testRandomMinpoly", stream.m ());

	bool ret = true;

	size_t j;
	typename Field::Element pi;
	Polynomial m_D;
	VectorDomain<Field> VD (F);

	Vector d;

	VectorWrapper::ensureDim (d, stream.n ());

	while (stream) {
		commentator.startIteration (stream.j ());

		F.init (pi, 1);

		stream.next (d);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal vector: ";
		VD.write (report, d);
		report << endl;

		for (j = 0; j < stream.n (); j++)
			F.mulin (pi, VectorWrapper::constRef<Field> (d, j));

		report << "Product: ";
		F.write (report, pi);
		report << endl;

		Blackbox D (F, d);
		minpoly (m_D, D);

		report << "Minimal polynomial: ";
		printPolynomial (F, report, m_D);

		if (!F.areEqual (m_D[0], pi)) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: m_D(0) != det(D)" << endl;
			ret = false;
		}

		commentator.stop (MSG_STATUS(ret));
		commentator.progress ();
	}

	stream.reset ();

	// try it with the random cstor of diagonal
	LinBox::Diagonal <Field> D(F, 10);
	unsigned long r;
	rank(r, D, Method::Wiedemann());
	if (r != 10)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: zeroes in random diagonal" << endl;
	ret = ret && r == 10;

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomMinpoly");

	return ret;
}

/* Test 3: Random linearity
 *
 * Compute a random diagonal matrix and use the linearity test in test-generic.h
 * to ensure that the diagonal black box does indeed correspond to a linear
 * mapping.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply matrix
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector>
static bool testRandomLinearity (Field &F,
				 VectorStream<std::vector<typename Field::Element> > &d_stream,
				 VectorStream<Vector> &stream1,
				 VectorStream<Vector> &stream2) 
{
	typedef LinBox::Diagonal <Field> Blackbox;

	commentator.start ("Testing random transpose", "testRandomLinearity", stream1.m ());

	VectorDomain<Field> VD (F);

	std::vector<typename Field::Element> d;
	VectorWrapper::ensureDim (d, stream1.n ());

	d_stream.next (d);
	Blackbox D (F, d);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	report << "Diagonal vector: ";
	VD.write (report, d) << endl;

	bool ret = testLinearity (F, D, stream1, stream2);

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomLinearity");

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
				 VectorStream<std::vector<typename Field::Element> > &d_stream,
				 VectorStream<Vector> &stream1,
				 VectorStream<Vector> &stream2) 
{
	typedef LinBox::Diagonal <Field> Blackbox;

	commentator.start ("Testing random transpose", "testRandomTranspose", stream1.m ());

	VectorDomain<Field> VD (F);

	std::vector<typename Field::Element> d;
	VectorWrapper::ensureDim (d, stream1.n ());

	d_stream.next (d);
	Blackbox D (F, d);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	report << "Diagonal vector: ";
	VD.write (report, d) << endl;

	bool ret = testTranspose (F, D, stream1, stream2);

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 2; // was 100

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	typedef Modular<LinBox::uint32> Field; //C.Pernet: avoids confusion with givaro::uint32
	typedef vector<Field::Element> Vector;

	parseArguments (argc, argv, args);
	Field F (q);

	srand (time (NULL));

	commentator.start("Diagonal matrix black box test suite", "diagonal");

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	RandomDenseStream<Field, Vector> stream1 (F, n, iterations), stream2 (F, n, iterations), d_stream (F, n, 1);
	RandomDenseStream<Field, Vector, NonzeroRandIter<Field> >
		stream3 (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, iterations);

	if (!testIdentityApply    (F, stream1)) pass = false;
	if (!testRandomMinpoly    (F, stream3)) pass = false;
	if (!testRandomLinearity  (F, d_stream, stream1, stream2)) pass = false;
	if (!testRandomTranspose  (F, d_stream, stream1, stream2)) pass = false;

        Field::RandIter iter(F);
	LinBox::Diagonal<Field> D(F, 10, iter);
	pass = pass && testBlackbox(D);

	commentator.stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
