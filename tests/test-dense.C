/* tests/test-dense.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Modified by Zhendong Wan <wan@cis.udel.edu>
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


/*! @file  tests/test-dense.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



// something currently commented out
#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/ring/modular.h"
#include "linbox/matrix/dense-matrix.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Identity matrix in dense representation
 *
 * Construct a dense representation of an n x n identity matrix and check
 * whether the output of its application to a series of random vectors is equal
 * to the input.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random vectors to which to apply identity inverse
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentity (Field &F, size_t n, int iterations = 1)
{
	typedef DenseVector<Field> Vector;
	typedef DenseMatrix<Field> Base;
	typedef DenseMatrix<Field>           Blackbox;

	commentator().start ("Testing identity apply", "testIdentity", (unsigned int)iterations);

	bool ret = true;
	// bool iter_passed = true;

	Blackbox I(F, n, n);
	// Matrix K(I);
	//typename Field::Element x; F.init(x);
	//F.write(std::cout, K.getEntry(x, i, j)) << std::endl;
	//Matrix L(K);


	for (size_t i = 0; i < n; i++)
		I.setEntry (i, i, F.one);

	Vector v(F,n), w(F,n);
	typename Field::RandIter r (F);

	for (int i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator().start (buf);

		bool iter_passed = true;

		for (size_t j = 0; j < n; j++)
			r.random (v[(size_t)j]);

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Input vector: ";
		printVector<Field> (F, report, v);

		I.apply (w, v);
		printVector<Field> (F, report, w);

		Base J (I);
		Blackbox KK( J);
		KK.apply (w, v);
		report << "Output vector: ";
		printVector<Field> (F, report, w);

		for (size_t j = 0; j < (size_t)n; j++)
			if (!F.areEqual (w[(size_t)j], v[(size_t)j]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testIdentity");

	return ret;
}

/* Test 2: Application of Vandermonde matrix in dense representation
 *
 * Computes a random Vandermonde matrix and applies it to a series of random
 * vectors. The random vectors contain the coefficients of polynomials over the
 * ground field. The output of the application is the result of evaluating these
 * polynomials at the points given by the second column of the matrix. This
 * function interpolates (using Lagrange interpolants) the evaluation points to
 * get the original polynomials and checks whether the coefficients match the
 * original vectors.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of random diagonal matrices to construct
 * N - Number of random vectors to which to apply random Vandermonde matrix
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testVandermonde (Field &F, size_t n, int iterations = 1, int N = 1)
{
	typedef DenseVector<Field> Vector;
	typedef DenseMatrix <Field> Blackbox;

	commentator().start ("Testing Vandermonde apply", "testVandermonde", (unsigned int)iterations);

	bool ret = true;
	bool inner_iter_passed;

	int i, j, k;

	Blackbox V(F, n, n); 

	Vector x(F,n), v(F,n), y(F,n), f(F,n);
	typename Field::RandIter r (F);
	typename Field::Element t;

	for (i = 0; i < iterations; i++) {
		char buf[80];
		snprintf (buf, 80, "Iteration %d", i);
		commentator().start (buf);

		/* Evaluation points */
		for (j = 0; j < (int) n; j++) {
			bool flag = true;

			// Make sure points are all distinct
			while (flag) {
				r.random (x[(size_t)j]);
				flag = false;
				for (k = 0; k < j; k++)
					if (F.areEqual (x[(size_t)j], x[(size_t)k]))
						flag = true;
			}
		}

		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Evaluation points: ";
		printVector<Field> (F, report, x);

		/* Build the Vandermonde matrix */
		for (j = 0; j < (int) n; j++) {
			F.assign(t, F.one);

			for (k = 0; k < (int) n; k++) {
				V.setEntry ((size_t)j,(size_t) k, t);
				F.mulin (t, x[(size_t)j]);
			}
		}

		for (j = 0; j < (int) N; j++) {
			inner_iter_passed = true;

			/* Random vector of evaluation results */
			for (k = 0; k < (int)n; k++)
				r.random (v[(size_t)k]);

			report << "Input vector: ";
			printVector<Field> (F, report, v);

			/* w should now be a vector of polynomial evaluations */
			V.apply (y, v);

			report << "Output vector: ";
			printVector<Field> (F, report, y);

			/* Polynomial interpolation to check whether w is correct */
			interpolatePoly (F, f, x, y);

			report << "Interpolation results: ";
			printVector<Field> (F, report, f);

			for (k = 0; k < (int) n; k++)
				if (!F.areEqual (f[(size_t)k], v[(size_t)k]))
					ret = inner_iter_passed = false;

			if (!inner_iter_passed)
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Vectors are not equal" << endl;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testVandermonde");

	return ret;
}

/* Test 3: Random linearity
 *
 * Construct a random dense matrix and a submatrix thereof. Call testLinearity
 * in test-generic.h to test that the submatrix is a linear operator
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrices
 * iterations - Number of iterations to run
 * N - Number of random vectors to which to apply
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomLinearity ( const Field & F,
				 VectorStream<DenseVector<Field> > &A_stream,
				 VectorStream<DenseVector<Field> > &v1_stream,
				 VectorStream<DenseVector<Field> > &v2_stream)
{
	commentator().start ("Testing random linearity", "testRandomLinearity", v1_stream.size ());

	DenseMatrix<Field> A (F, A_stream);

	bool ret = testLinearity (A, v1_stream, v2_stream);

	A_stream.reset ();
	v1_stream.reset ();
	v2_stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandomLinearity");

	return ret;
}

/* Test 4: Random transpose
 *
 * Construct a random dense matrix and a submatrix thereof. Call testLinearity
 * in test-generic.h to test that the submatrix is a linear operator
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrices
 * iterations - Number of iterations to run
 * N - Number of random vectors to which to apply
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomTranspose (const Field                                 &F,
				 VectorStream<DenseVector<Field> > &A_stream,
				 VectorStream<DenseVector<Field> > &v1_stream,
				 VectorStream<DenseVector<Field> > &v2_stream)
{
	commentator().start ("Testing random transpose", "testRandomTranspose", v1_stream.size ());

	DenseMatrix<Field> A (F, A_stream);

	bool ret = testTranspose (F, A, v1_stream, v2_stream);

	A_stream.reset ();
	v1_stream.reset ();
	v2_stream.reset ();

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}


int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;
	static int iterations = 2; // was 100
	//static int N = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",   TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	typedef Givaro::Modular<uint32_t> Field;

	parseArguments (argc, argv, args);
	Field F (q); Field::RandIter gen(F);

	commentator().start("Dense matrix black box test suite", "DenseMatrix");

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	RandomDenseStream<Field> A_stream (F, gen, n, n);
	RandomDenseStream<Field> v1_stream (F, gen, n, (unsigned int)iterations);
	RandomDenseStream<Field> v2_stream (F, gen, n, (unsigned int)iterations);

	if (!testIdentity    (F, n)) pass = false;
	if (!testVandermonde (F, n)) pass = false;
	DenseMatrix<Field> A(F, A_stream);
	if (!testBlackbox(A)) pass = false;
	//if (!testRandomLinearity (F, A_stream, v1_stream, v2_stream)) pass = false;
	//if (!testRandomTranspose (F, A_stream, v1_stream, v2_stream)) pass = false;

	commentator().stop("dense matrix black box test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
