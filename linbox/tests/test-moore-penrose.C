/* -*- mode: c; style: linux -*- */

/* tests/test-moore-penrose.C
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/blackbox/sparse0.h"
#include "linbox/blackbox/moore-penrose.h"
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
 * n - Row dimension of matrix
 * m - Column dimension of matrix
 * r - Rank of matrix
 * factory - Factory for random vectors to apply
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityApply (Field                                           &F,
			       size_t                                           n,
			       size_t                                           m,
			       size_t                                           r,
			       VectorFactory<vector<typename Field::Element> > &factory) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing random apply", "testIdentityApply", factory.m ());

	bool ret = true;
	bool iter_passed;

	Vector v, w;

	VectorWrapper::ensureDim (w, factory.n ());

	int i, j, k, l;

	Blackbox A (F, n, m);

	typename Field::Element x;

	F.init (x, 1);

	for (i = 0; i < r; i++)
		A.put_value (pair<size_t, size_t> (i, i), x);

	MoorePenrose<Field, Vector> Adagger (F, &A, r);

	while (factory) {
		commentator.startIteration (i);

		iter_passed = true;

		factory.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

		report << "Input vector:  ";
		printVector<Field> (F, report, v);

		Adagger.apply (w, v);

		report << "Output vector: ";
		printVector<Field> (F, report, w);

		for (l = 0; l < r; l++)
			if (!F.areEqual (v[l], w[l]))
				ret = iter_passed = false;

		for (l = r; l < n; l++)
			if (!F.isZero (w[l]))
				ret = iter_passed = false;

		if (!iter_passed)
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Output vector is incorrect" << endl;

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testIdentityApply");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 100;
	static size_t m = 100;
	static size_t r = 10;
	static integer q = 101;
	static int iterations = 100;

	static Argument args[] = {
		{ 'n', "-n N", "Set row dimension of test matrices to N (default 100)",    TYPE_INT,     &n },
		{ 'm', "-m M", "Set column dimension of test matrices to M (default 100)", TYPE_INT,     &m },
		{ 'r', "-r R", "Set rank of test matrices to R (default 10)",              TYPE_INT,     &r },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 101)",       TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",         TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "MoorePenrose black box test suite" << endl << endl;

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

	RandomDenseVectorFactory<Modular<long> > factory (F, m, iterations);

	if (!testIdentityApply<Modular<long> > (F, n, m, r, factory)) pass = false;

	return pass ? 0 : -1;
}
