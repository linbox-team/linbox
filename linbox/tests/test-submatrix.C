/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#include "linbox-config.h"

#include <iostream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/blackbox/dense.h"
#include "linbox/util/vector-factory.h"

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
static bool testRandomApply (Field                                           &F,
			     unsigned int                                     iterations,
			     size_t                                           n,
			     VectorFactory<vector<typename Field::Element> > &factory) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef DenseMatrix <Field> Blackbox;

	commentator.start ("Testing random apply", "testRandomApply", iterations);

	bool ret = true;
	bool iter_passed;

	Vector v, w1(n), w2(n);

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

		factory.reset ();

		while (factory) {
			ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

			factory.next (v);

			report << "Input vector: ";
			printVector<Field> (F, report, v);

			for (k = 0; k < 9; k++) {
				commentator.indent (report);
				report << "Checking section " << k / 3 + 1 << "x" << k % 3 + 1 << endl;

				Submatrix<Vector> B (&A, n * (k / 3), n * (k % 3), n, n);
				B.apply (w1, v);

				commentator.indent (report);
				report << "Output vector (computed): ";
				printVector<Field> (F, report, w1);

				Ai[k]->apply (w2, v);

				commentator.indent (report);
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

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply");

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
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 10)", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",   TYPE_INT,     &iterations },
		{ 'N', "-N N", "Perform each test on N vectors (default 1)",         TYPE_INT,     &N },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "Submatrix matrix black box test suite" << endl << endl;

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

	RandomDenseVectorFactory<Modular<long> > factory (F, n, N);

	if (!testRandomApply<Modular<long> > (F, iterations, n, factory)) pass = false;

	return pass ? 0 : -1;
}
