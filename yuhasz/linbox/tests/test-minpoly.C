/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-minpoly.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 * 2002-04-03: William J. Turner <wjturner@acm.org>
 *
 * changed name of sparse-matrix file.
 *
 * --------------------------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/vector/stream.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Minimal polynomial of the identity matrix
 *
 * Construct the identity matrix and compute its minimal polynomial. Confirm
 * that the result is 1-x
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityMinpoly (Field &F, size_t n) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef pair <vector <size_t>, vector <typename Field::Element> > Row;
	typedef SparseMatrix <Field> Blackbox;

	commentator.start ("Testing identity minpoly", "testIdentityMinpoly");

	bool ret = true;

	typename Field::Element c0, c1;

	StandardBasisStream<Field, Row> stream (F, n);
	Blackbox A (F, stream);

	Polynomial phi;

	minpoly (phi, A, F);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Minimal polynomial is: ";
	printPolynomial<Field, Polynomial> (F, report, phi);

	F.init (c0, -1);
	F.init (c1, 1);

	if (phi.size () != 2 ||
	    !F.areEqual (phi[0], c0) ||
	    !F.areEqual (phi[1], c1)) {
		ret = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Minimal polynomial is incorrect" << endl;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityMinpoly");

	return ret;
}

/* Test 2: Minimal polynomial of a nilpotent matrix
 *
 * Construct an index-n nilpotent matrix and compute its minimal
 * polynomial. Confirm that the result is x^n
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testNilpotentMinpoly (Field &F, size_t n) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef pair <vector <size_t>, vector <typename Field::Element> > Row;
	typedef SparseMatrix <Field> Blackbox;

	commentator.start ("Testing nilpotent minpoly", "testNilpotentMinpoly");

	bool ret = true;
	bool lowerTermsCorrect = true;

	size_t i;

	StandardBasisStream<Field, Row> stream (F, n);
	Row v;
	stream.next (v);
	Blackbox A (F, stream);

	Polynomial phi;

	minpoly (phi, A, F);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << "Minimal polynomial is: ";
	printPolynomial (F, report, phi);

	for (i = 0; i < n - 1; i++)
		if (!F.isZero (phi[i]))
			lowerTermsCorrect = false;

	if (phi.size () != n + 1 || !F.isOne (phi[n]) || !lowerTermsCorrect) {
		ret = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: Minimal polynomial is incorrect" << endl;
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNilpotentMinpoly");

	return ret;
}

/* Test 3: Random minpoly of sparse matrix
 *
 * Generates a random sparse matrix with K nonzero elements per row and computes
 * its minimal polynomial. Then computes random vectors and applies the
 * polynomial to them in Horner style, checking whether the result is 0.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * K - Number of nonzero elements per row
 * numVectors - Number of random vectors to which to apply the minimal polynomial
 *
 * Return true on success and false on failure
 */

template <class Field, class Row, class Vector>
bool testRandomMinpoly (Field                 &F,
			int                    iterations,
			VectorStream<Row>    &A_stream,
			VectorStream<Vector> &v_stream) 
{
	typedef vector <typename Field::Element> Polynomial;
	typedef SparseMatrix <Field> Blackbox;

	commentator.start ("Testing sparse random minpoly", "testRandomMinpoly", iterations);

	bool ret = true;
	bool iter_passed;

	VectorDomain<Field> VD (F);

	Vector v, w;

	VectorWrapper::ensureDim (v, v_stream.n ());
	VectorWrapper::ensureDim (w, v_stream.n ());

	for (int i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		A_stream.reset ();
		Blackbox A (F, A_stream);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Matrix:" << endl;
		A.write (report, Blackbox::FORMAT_PRETTY);

		Polynomial phi;

		minpoly (phi, A, F);

		report << "Minimal polynomial is: ";
		printPolynomial (F, report, phi);

		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "deg minpoly (A) = " << phi.size () - 1 << endl;

		v_stream.reset ();

		while (v_stream) {
			v_stream.next (v);

			report << "Input vector  " << v_stream.j () << ": ";
			VD.write (report, v);
			report << endl;

			applyPoly (F, w, A, phi, v);

			report << "Output vector " << v_stream.j () << ": ";
			VD.write (report, w);
			report << endl;

			if (!VD.isZero (w))
				ret = iter_passed = false;
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector was incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomMinpoly");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 10;
	static int numVectors = 100;
	static int k = 3;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)",                 TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)",          TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",                    TYPE_INT,     &iterations },
		{ 'v', "-v V", "Use V test vectors for the random minpoly tests (default 100)",      TYPE_INT,     &numVectors },
		{ 'k', "-k K", "K nonzero Elements per row in sparse random apply test (default 3)", TYPE_INT,     &k },
	};

	typedef Modular<uint32> Field;
	typedef vector<Field::Element> DenseVector;
	typedef pair<vector<size_t>, vector<Field::Element> > SparseVector;

	parseArguments (argc, argv, args);
	Field F (q);

	srand (time (NULL));

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	cout << endl << "Black box minimal polynomial test suite" << endl;

	RandomDenseStream<Field, DenseVector, NonzeroRandIter<Field> >
		v_stream (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, numVectors);
	RandomSparseStream<Field, SparseVector, NonzeroRandIter<Field> >
		A_stream (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), (double) k / (double) n, n, n);

	if (!testIdentityMinpoly  (F, n)) pass = false;
	if (!testNilpotentMinpoly (F, n)) pass = false;
	if (!testRandomMinpoly    (F, iterations, A_stream, v_stream)) pass = false;

	return pass ? 0 : -1;
}
