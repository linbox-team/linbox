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

#include "linbox/field/modular.h"
#include <linbox/field/gmp-integers.h>
#include "linbox/blackbox/sparse.h"
#include "linbox/util/commentator.h"
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

template <class Field, class Meth>
static bool testIdentityMinpoly (Field &F, size_t n, bool symmetrizing, const Meth& M) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef SparseMatrix<Field> Blackbox;
	typedef typename Blackbox::Row Row;

	commentator.start ("Testing identity minpoly", "testIdentityMinpoly");

	bool ret = true;

	typename Field::Element c0, c1;

	StandardBasisStream<Field, Row> stream (F, n);
	Blackbox A (F, stream);

	Polynomial phi;

	//if (symmetrizing) minpolySymmetric (phi, A);
	//else minpoly (phi, A);
	minpoly (phi, A, M );

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

template <class Field>
static bool testIdentityMinpoly (Field &F, size_t n, bool symmetrizing=false) 
{
    return testIdentityMinpoly(F, n, symmetrizing, Method::Wiedemann());
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

template <class Field, class Meth>
static bool testNilpotentMinpoly (Field &F, size_t n, const Meth& M)
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef SparseMatrix <Field> Blackbox;
	typedef typename Blackbox::Row Row;

	commentator.start ("Testing nilpotent minpoly", "testNilpotentMinpoly");

	bool ret = true;
	bool lowerTermsCorrect = true;

	size_t i;

	StandardBasisStream<Field, Row> stream (F, n);
	Row v;
	stream.next (v);
	Blackbox A (F, stream); // first subdiagonal is 1's.

	Polynomial phi(n+1);

	minpoly (phi, A, M);

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
template <class Field>
static bool testNilpotentMinpoly (Field &F, size_t n)
{
    return testNilpotentMinpoly(F, n, Method::Wiedemann());
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

template <class Field, class BBStream, class VectStream, class Meth>
bool testRandomMinpoly (Field                 &F,
			int                    iterations,
			BBStream    &A_stream,
			VectStream &v_stream,
                        const Meth& M)
{
	typedef std::vector <typename Field::Element> Polynomial;
	typedef SparseMatrix <Field> Blackbox;
        typedef typename VectStream::Vector Vector;

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
		A.write (report, FORMAT_MAPLE);

		Polynomial phi;

		minpoly (phi, A, M );

		report << "Minimal polynomial is: ";
		printPolynomial (F, report, phi);

		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
			<< "deg minpoly (A) = " << phi.size () - 1 << endl;

		v_stream.reset ();

		while (iter_passed && v_stream) {
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
template <class Field, class BBStream, class VectStream>
bool testRandomMinpoly (Field                 &F,
			int                    iterations,
			BBStream    &A_stream,
			VectStream &v_stream)
{
    return testRandomMinpoly(F, iterations, A_stream, v_stream, Method::Wiedemann());
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


	parseArguments (argc, argv, args);

// /////////////// finite field part //////////////////
	typedef Modular<LinBox::uint32> Field;
	typedef vector<Field::Element> DenseVector;
	typedef SparseMatrix<Field>::Row SparseVector;
	//typedef pair<vector<size_t>, vector<Field::Element> > SparseVector;
	Field F (q);
	srand (time (NULL));

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	cout << endl << "Wiedemann minimal polynomial of a matrix over a prime field test suite" << endl;

	RandomDenseStream<Field, DenseVector, NonzeroRandIter<Field> >
		v_stream (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), n, numVectors);
	RandomSparseStream<Field, SparseVector, NonzeroRandIter<Field> >
		A_stream (F, NonzeroRandIter<Field> (F, Field::RandIter (F)), (double) k / (double) n, n, n);

	//no symmetrizing
	if (!testIdentityMinpoly  (F, n)) pass = false;
	if (!testNilpotentMinpoly (F, n)) pass = false;
	if (!testRandomMinpoly    (F, iterations, A_stream, v_stream)) pass = false;

	// symmetrizing
	if (!testIdentityMinpoly  (F, n, true)) pass = false;
	//need other tests...

	cout << endl << "minimal polynomial (basic methods) of a matrix over a prime field test suite" << endl;

std::cout << "Hybrid" << std::endl;
	if (!testIdentityMinpoly  (F, n, false,  Method::Hybrid())) pass = false;
	if (!testNilpotentMinpoly (F, n, Method::Hybrid())) pass = false;
std::cout << "Blackbox" << std::endl;
	if (!testIdentityMinpoly  (F, n, false,  Method::Blackbox())) pass = false;
	if (!testNilpotentMinpoly (F, n, Method::Blackbox())) pass = false;
std::cout << "Elimination" << std::endl;
	if (!testIdentityMinpoly  (F, n, false,  Method::Elimination())) pass = false;
	if (!testNilpotentMinpoly (F, n, Method::Elimination())) pass = false;

// /////////////// integer part //////////////////
	typedef vector<PID_integer::Element> ZDenseVector;
	typedef SparseMatrix<PID_integer>::Row ZSparseVector;
	//typedef pair<vector<size_t>, vector<Field::Element> > SparseVector;
	PID_integer Z;
	srand (time (NULL));

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	cout << endl << "Wiedemann minimal polynomial of an integer matrix test suite" << endl;

	RandomDenseStream<PID_integer, ZDenseVector, NonzeroRandIter<PID_integer> >
		zv_stream (Z, NonzeroRandIter<PID_integer> (Z, PID_integer::RandIter (Z)), n, numVectors);
	RandomSparseStream<PID_integer, ZSparseVector, NonzeroRandIter<PID_integer> >
		zA_stream (Z, NonzeroRandIter<PID_integer> (Z, PID_integer::RandIter (Z)), (double) k / (double) n, n, n);

            // Hybrid
	//no symmetrizing
	if (!testIdentityMinpoly  (Z, n)) pass = false;
	if (!testNilpotentMinpoly (Z, n)) pass = false;

	if (!testRandomMinpoly    (Z, iterations, zA_stream, zv_stream)) pass = false;

	// symmetrizing
	if (!testIdentityMinpoly  (Z, n, true)) pass = false;

	cout << endl << "minimal polynomial (basic methods) of an integer matrix test suite" << endl;

std::cout << "Hybrid" << std::endl;
	if (!testIdentityMinpoly  (Z, n, false,  Method::Hybrid())) pass = false;
	if (!testNilpotentMinpoly (Z, n, Method::Hybrid())) pass = false;
std::cout << "Blackbox" << std::endl;
	if (!testIdentityMinpoly  (Z, n, false,  Method::Blackbox())) pass = false;
	if (!testNilpotentMinpoly (Z, n, Method::Blackbox())) pass = false;
std::cout << "Elimination" << std::endl;
	if (!testIdentityMinpoly  (Z, n, false,  Method::Elimination())) pass = false;
	if (!testNilpotentMinpoly (Z, n, Method::Elimination())) pass = false;

	return pass ? 0 : -1;
}
