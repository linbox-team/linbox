
/* tests/test-submatrix.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/blackbox/dense.h"
#include "linbox/vector/stream.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Test 1: Application of submatrices onto random vectors
 *
 * Construct 9 dense matrices and place them, as submatrices, into a 3x3 dense
 * grid. Apply submatrices of the result onto random vectors and check equality
 * with the result of applying the submatrices directly.
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrices
 * iterations - Number of iterations to run
 * N - Number of random vectors to which to apply
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomApply (Field                                       &F,
			     unsigned int                                 iterations,
			     size_t                                       n,
			     VectorStream<typename Vector<Field>::Dense> &stream) 
{
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing random apply", "testRandomApply", iterations);

	bool ret = true;
	bool iter_passed;

	typename Vector<Field>::Dense v, w1(n), w2(n);

	size_t i, j, k, l;

	Blackbox *Ai[9];
	Blackbox A (F, n * 3, n * 3);

	for (i = 0; i < 9; i++)
		Ai[i] = new Blackbox (F, n, n);

	typename Field::Element x;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		iter_passed = true;

		for (j = 0; j < 9; j++) {
			for (k = 0; k < n; k++) {
				for (l = 0; l < n; l++) {
					r.random (x);
					Ai[j]->setEntry (k, l, x);
					A.setEntry (k + n * (j / 3), l + n * (j % 3), x);
				}
			}
		}

		stream.reset ();

		while (stream) {
			ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

			stream.next (v);

			report << "Input vector: ";
			printVector<Field> (F, report, v);

			for (k = 0; k < 9; k++) {
				report << "Checking section " << k / 3 + 1 << "x" << k % 3 + 1 << endl;

				Submatrix<Blackbox> B (F, &A, n * (k / 3), n * (k % 3), n, n);
				B.apply (w1, v);

				report << "Output vector (computed): ";
				printVector<Field> (F, report, w1);

				Ai[k]->apply (w2, v);

				report << "Output vector (true):     ";
				printVector<Field> (F, report, w2);

				for (l = 0; l < n; l++)
					if (!F.areEqual (w1[l], w2[l]))
						ret = iter_passed = false;
			}
		}

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Vectors are not equal" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	stream.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply");

	return ret;
}

/* Test 2: Random linearity
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
static bool testRandomLinearity (const Field                                 &F,
				 VectorStream<typename Vector<Field>::Dense> &A_stream,
				 VectorStream<typename Vector<Field>::Dense> &v1_stream,
				 VectorStream<typename Vector<Field>::Dense> &v2_stream) 
{
	commentator.start ("Testing random linearity", "testRandomLinearity", v1_stream.size ());

	DenseMatrix<Field> A (F, A_stream);
	Submatrix<DenseMatrix<Field> > Ap (&A, 0, 0, v1_stream.dim (), v2_stream.dim ());

	bool ret = testLinearity (F, Ap, v1_stream, v2_stream);

	A_stream.reset ();
	v1_stream.reset ();
	v2_stream.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomLinearity");

	return ret;
}

/* Test 3: Random transpose
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
				 VectorStream<typename Vector<Field>::Dense> &A_stream,
				 VectorStream<typename Vector<Field>::Dense> &v1_stream,
				 VectorStream<typename Vector<Field>::Dense> &v2_stream) 
{
	commentator.start ("Testing random transpose", "testRandomTranspose", v1_stream.size ());

	DenseMatrix<Field> A (F, A_stream);
	Submatrix<DenseMatrix<Field> > Ap (&A, 0, 0, v1_stream.dim (), v2_stream.dim ());

	bool ret = testTranspose (F, Ap, v1_stream, v2_stream);

	A_stream.reset ();
	v1_stream.reset ();
	v2_stream.reset ();

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomTranspose");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;
	static int iterations = 100;
	static int N = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 'N', "-N N", "Perform each test on N vectors.", TYPE_INT,     &N },
		{ '\0' }
	};

	typedef Modular<uint32> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	commentator.start("Submatrix black box test suite", "Submatrix");

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	RandomDenseStream<Field> stream (F, n, N);
	RandomDenseStream<Field> A_stream (F, n, n);
	RandomDenseStream<Field> v1_stream (F, n / 2, iterations);
	RandomDenseStream<Field> v2_stream (F, n / 2, iterations);

#if 0
	if (!testRandomApply (F, iterations, n, stream)) pass = false;
#endif
	if (!testRandomLinearity (F, A_stream, v1_stream, v2_stream)) pass = false;
	if (!testRandomTranspose (F, A_stream, v1_stream, v2_stream)) pass = false;

	commentator.stop("Submatrix black box test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
