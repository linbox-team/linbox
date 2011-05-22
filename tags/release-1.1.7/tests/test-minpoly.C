
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

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/field/modular.h"
#include "linbox/field/givaro-gfq.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/util/commentator.h"
#include "linbox/solutions/minpoly.h"
#include "linbox/vector/stream.h"

#include "test-common.h"

using namespace LinBox;

/* Test 0: Minimal polynomial of the zero matrix
*/
template <class Field, class Meth>
static bool testZeroMinpoly (Field &F, size_t n, bool symmetrizing, const Meth& M) 
{
	commentator.start ("Testing zero minpoly", "testZeroMinpoly");
	typedef vector <typename Field::Element> Polynomial;
	Polynomial phi;
	SparseMatrix<Field> A(F, n, n);
	minpoly(phi, A, M);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        A.write(report, FORMAT_MAPLE);
	report << "Minimal polynomial is: ";

	printPolynomial<Field, Polynomial> (F, report, phi);

	bool ret;
	if (phi.size () == 2 && F.isZero (phi[0]) && F.isOne(phi[1]) ) 
		ret = true;
	else 
	{
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A = 0, should get x, got ";
		printPolynomial<Field, Polynomial> (F, report, phi);
		report << endl;
		ret = false;
	}
	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testZeroMinpoly");
	return ret;
}
template <class Field>
static bool testZeroMinpoly (Field &F, size_t n) 
{
	return testZeroMinpoly(F, n, false, Method::Blackbox());
}

/* Test 1: Minimal polynomial of the identity matrix
 *
 * Construct the identity matrix and compute its minimal polynomial. Confirm
 * that the result is x-1
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

	bool ret;

	F.init (c0, -1);
	F.init (c1, 1);

	if (phi.size () == 2 && F.areEqual (phi[0], c0) && F.areEqual (phi[1], c1)) 
		ret = true;
	else {
		ret = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A = I, should get x-1, got ";
		printPolynomial<Field, Polynomial> (F, report, phi);
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityMinpoly");

	return ret;
}

template <class Field>
static bool testIdentityMinpoly (Field &F, size_t n, bool symmetrizing=false) 
{
    return testIdentityMinpoly(F, n, symmetrizing, Method::Blackbox());
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

	bool ret;
	bool lowerTermsCorrect = true;

	size_t i;

	StandardBasisStream<Field, Row> stream (F, n);
	Row v;
	stream.next (v);
	Blackbox A (F, stream); // first subdiagonal is 1's.

	Polynomial phi(n+1);

	minpoly (phi, A, M);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        A.write (report, FORMAT_MAPLE);
	report << "Minimal polynomial is: ";
	printPolynomial (F, report, phi);

	for (i = 0; i < n - 1; i++)
		if (!F.isZero (phi[i]))
			lowerTermsCorrect = false;

	if (phi.size () == n + 1 && F.isOne (phi[n]) && lowerTermsCorrect)
		ret = true;
	else {
		ret = false;
		commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "ERROR: A^n = 0, should get x^" << n <<", got ";
		printPolynomial (F, report, phi);
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testNilpotentMinpoly");

	return ret;
}
template <class Field>
static bool testNilpotentMinpoly (Field &F, size_t n)
{
    return testNilpotentMinpoly(F, n, Method::Blackbox());
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
				<< "ERROR: A = rand, purported-minpoly(A) is not zero." << endl;
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
    return testRandomMinpoly(F, iterations, A_stream, v_stream, Method::Blackbox());
}

/* Test 4: test gram matrix.
     A test of behaviour with self-orthogonal rows and cols.

	 Gram matrix is 1's offdiagonal and 0's on diagonal.  Self orthogonality when characteristic | n-1.

	 Arg m is ignored. 
	 if p := characteristic is small then size is p+1 and minpoly is x^2 + x
	 (because A^2 = (p-1)A).
	 if p is large, size is 2 and minpoly is x^2 -1.
*/
template <class Field, class Meth>
static bool testGramMinpoly (Field &F, size_t m, bool symmetrizing, const Meth& M) 
{
	commentator.start ("Testing gram minpoly", "testGramMinpoly");
	typedef vector <typename Field::Element> Polynomial;
	integer n;
	F.characteristic(n); n += 1;
	if (n > 30) n = 2;
	Polynomial phi;
	typename Field::Element one, zero, neg1; F.init(one, 1); F.init(zero, 0); F.init(neg1); F.neg(neg1, one);
	DenseMatrix<Field> A(F, n, n);
	for (int i = 0; i < n; ++i) for (int j = 0; j < n; ++j) A.setEntry(i, j, one);
	for (int i = 0; i < n; ++i) A.setEntry(i, i, zero);
	minpoly(phi, A, M);

	ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        A.write (report);
	report << "Minimal polynomial is: ";
	printPolynomial<Field, Polynomial> (F, report, phi);

	bool ret;
	if (n == 2)
		if ( phi.size() == 3 && F.areEqual(phi[0], neg1) && F.isZero(phi[1]) && F.isOne(phi[2]))
			ret = true;
		else
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A = gram, should get x^2 - x, got ";
			printPolynomial<Field, Polynomial> (F, report, phi);
			ret = false;
		}
	else 
		if (phi.size() == 3 && F.isZero(phi[0]) && F.isOne(phi[1]) && F.isOne(phi[2]))
			ret = true;
		else
		{
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: A = gram, should get x^2 + x, got ";
			printPolynomial<Field, Polynomial> (F, report, phi);
			ret = false;
		}
		commentator.stop (MSG_STATUS (ret), (const char *) 0, "testGramMinpoly");
		return ret;
}

template <class Field>
static bool testGramMinpoly (Field &F, size_t n) 
{
	return testGramMinpoly(F, n, false, Method::Blackbox());
}


int main (int argc, char **argv)
{
        commentator.setMaxDetailLevel (-1);

        commentator.setMaxDepth (-1);
	bool pass = true;

	static size_t n = 10;
	//static integer q = 2147483647U;
	static integer q = 1000003; // ok for both Modular<int> and Modular<double>
	static int iterations = 1;
	static int numVectors = 1;
	static int k = 3;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 'v', "-v V", "Use V test vectors for the random minpoly tests.", TYPE_INT,     &numVectors },
		{ 'k', "-k K", "K nonzero Elements per row in sparse random apply test.", TYPE_INT,     &k },
		{ '\0' }
	};


	parseArguments (argc, argv, args);
	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

// /////////////// finite field part //////////////////
	if (q > 5 && q % 2 != 0 && q % 3 != 0 && q % 5 != 0 )
	{
	//typedef Modular<LinBox::uint32> Field;
	//typedef Modular<int> Field;
	typedef Modular<double> Field;
	Field F (q);
	srand (time (NULL));

	commentator.start("Blackbox prime field minpoly test suite", "Wminpoly");

	//no symmetrizing
	if (!testZeroMinpoly  	  (F, n)) pass = false;
	if (!testIdentityMinpoly  (F, n)) pass = false;
	if (!testNilpotentMinpoly (F, n)) pass = false;
	//if (!testRandomMinpoly    (F, n)) pass = false;
	if (!testGramMinpoly      (F, n)) pass = false;

	// symmetrizing
	//if (!testZeroMinpoly  	  (F, n, true)) pass = false;
	if (!testIdentityMinpoly  (F, n, true)) pass = false;
	//if (!testNilpotentMinpoly (F, n, true)) pass = false;
	//if (!testRandomMinpoly    (F, n, true)) pass = false;
	//if (!testGramMinpoly      (F, n, true)) pass = false;
	//need other tests...

	commentator.stop("Blackbox prime field minpoly test suite");
	}else{

	int p;  
	if (q % 2 == 0) p = 2;
	if (q % 3 == 0) p = 3;
	if (q % 5 == 0) p = 5;
	int e = 0;  do {++e; q = q/p; } while (q > 1);
	typedef GivaroGfq Field;
	Field F (p, e);
	srand (time (NULL));

	commentator.start("Blackbox non-prime field minpoly test suite", "Wminpoly");
        ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        F.write(report);

	//no symmetrizing
	if (!testZeroMinpoly  	  (F, n)) pass = false;
	if (!testIdentityMinpoly  (F, n)) pass = false;
	if (!testNilpotentMinpoly (F, n)) pass = false;
	//if (!testRandomMinpoly    (F, n)) pass = false;
	if (!testGramMinpoly      (F, n)) pass = false;

	// symmetrizing
	//if (!testZeroMinpoly  	  (F, n, true)) pass = false;
	if (!testIdentityMinpoly  (F, n, true)) pass = false;
	//if (!testNilpotentMinpoly (F, n, true)) pass = false;
	//if (!testRandomMinpoly    (F, n, true)) pass = false;
	//if (!testGramMinpoly      (F, n, true)) pass = false;
	//need other tests...

	commentator.stop("Blackbox non-prime field minpoly test suite");
	}

#if 0

        Modular<LinBox::uint32> F (q);


	commentator.start("Hybrid prime field minpoly test suite", "Hminpoly");
	if (!testIdentityMinpoly  (F, n, false,  Method::Hybrid())) pass = false;
	if (!testNilpotentMinpoly (F, n, Method::Hybrid())) pass = false;
	commentator.stop("Hybrid prime field minpoly test suite");

	commentator.start("Blackbox prime field minpoly test suite", "Bminpoly");
	if (!testIdentityMinpoly  (F, n, false,  Method::Blackbox())) pass = false;
	if (!testNilpotentMinpoly (F, n, Method::Blackbox())) pass = false;
	commentator.stop("Blackbox prime field minpoly test suite");

	commentator.start("Elimination prime field minpoly test suite", "Eminpoly");
	if (!testIdentityMinpoly  (F, n, false,  Method::Elimination())) pass = false;
	if (!testNilpotentMinpoly (F, n, Method::Elimination())) pass = false;
	commentator.stop("Elimination prime field minpoly test suite");

// /////////////// integer part //////////////////
	typedef vector<PID_integer::Element> ZDenseVector;
	typedef SparseMatrix<PID_integer>::Row ZSparseVector;
	//typedef pair<vector<size_t>, vector<Field::Element> > SparseVector;
	PID_integer Z;
	srand (time (NULL));

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator.start("Blackbox integer minpoly test suite", "WIminpoly");

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

	commentator.stop("Blackbox integer minpoly test suite");

	commentator.start("Hybrid integer minpoly test suite", "HIminpoly");
	if (!testIdentityMinpoly  (Z, n, false,  Method::Hybrid())) pass = false;
	if (!testNilpotentMinpoly (Z, n, Method::Hybrid())) pass = false;
	commentator.stop("Hybrid integer minpoly test suite");

	commentator.start("Blackbox integer minpoly test suite", "BIminpoly");
	if (!testIdentityMinpoly  (Z, n, false,  Method::Blackbox())) pass = false;
	if (!testNilpotentMinpoly (Z, n, Method::Blackbox())) pass = false;
	commentator.stop("Blackbox integer minpoly test suite");

	commentator.start("Elimination integer minpoly test suite", "EIminpoly");
	if (!testIdentityMinpoly  (Z, n, false,  Method::Elimination())) pass = false;
	if (!testNilpotentMinpoly (Z, n, Method::Elimination())) pass = false;
	commentator.stop("Elimination integer minpoly test suite");

#endif
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
