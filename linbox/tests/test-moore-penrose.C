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
#include "linbox/solutions/rank.h"
#include "linbox/util/vector-factory.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Build a random sparse n x m matrix of rank r with a nonsingular leading principal minor */

template <class Field>
static SparseMatrix0<Field, vector<pair <size_t, typename Field::Element> >, vector<typename Field::Element> >
	*buildRandomSparseMatrix (Field &F, size_t n, size_t m, size_t r, double K, vector<typename Field::Element> &d) 
{
	SparseMatrix0<Field, vector<pair <size_t, typename Field::Element> >, vector<typename Field::Element> >
		*A = new SparseMatrix0<Field, vector<pair <size_t, typename Field::Element> >, vector<typename Field::Element> >
		(F, n, m);
	typename Field::RandIter rnd (F);
	typename Field::Element x, factor, val;
	vector<map<size_t, typename Field::Element> > bottom_left_data (n - r);
	vector<map<size_t, typename Field::Element> > top_right_data (r);
	int i, j, k;
	map<size_t, typename Field::Element>::iterator li1, li2;

	// Build diagonal part
	for (i = 0; i < r; i++) {
		do rnd.random (x); while (F.isZero (x));
		A->put_value (pair<size_t, size_t> (i, i), x);
		d[i] = x;
	}

	// Build top right part
	for (k = 0; k < K * r * (m - r); k++) {
		rnd.random (x);
		i = random () % r;
		j = random () % (m - r) + r;
		A->put_value (pair<size_t, size_t> (i, j), x);
		top_right_data[i][j] = x;
	}

	// Build bottom left part
	for (k = 0; k < K * r * (n - r); k++) {
		rnd.random (x);
		i = random () % (n - r) + r;
		j = random () % r;
		A->put_value (pair<size_t, size_t> (i, j), x);
		bottom_left_data[i - r][j] = x;
	}

	// Fill in bottom right part
	for (i = 0; i < n - r; i++) {
		map<size_t, typename Field::Element> row_data;
		map<size_t, typename Field::Element>::iterator mi;

		for (li1 = bottom_left_data[i].begin (); li1 != bottom_left_data[i].end (); li1++) {
			F.div (factor, (*li1).second, d[(*li1).first]);

			for (li2 = top_right_data[(*li1).first].begin ();
			     li2 != top_right_data[(*li1).first].end ();
			     li2++)
			{
				F.mul (val, factor, (*li2).second);
				F.addin (row_data[(*li2).first], val);
			}
		}

		for (mi = row_data.begin (); mi != row_data.end (); mi++)
			A->put_value (pair<size_t, size_t> (i + r, (*mi).first), (*mi).second);
	}

	return A;
}

/* Test 1: Application of identity onto random vectors
 *
 * Construct a rank-r matrix where the r x r leading principal minor is the
 * identity matrix, and all other entries are 0. Apply the Moore-Penrose inverse
 * of this to random vectors and check the correctness of the results.
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

	commentator.start ("Testing identity apply", "testIdentityApply", factory.m ());

	bool ret = true;
	bool iter_passed;

	Vector v, w;

	VectorWrapper::ensureDim (v, factory.n ());
	VectorWrapper::ensureDim (w, m);

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

		for (l = r; l < m; l++)
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

/* Test 2: Application of random matrix onto random vectors
 *
 * Construct a random sparse rank-r matrix with a nonsingular leading principal
 * minor. Apply the Moore-Penrose inverse of this to random vectors and check
 * the results.
 *
 * F - Field over which to perform computations
 * n - Row dimension of matrix
 * m - Column dimension of matrix
 * r - Rank of matrix
 * iterations - Number of iterations to run
 * K proportion of entries to make nonzero
 * factory - Factory for random vectors to apply
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testRandomApply1 (Field                                           &F,
			      size_t                                           n,
			      size_t                                           m,
			      size_t                                           r,
			      unsigned                                         iterations,
			      double                                           K,
			      VectorFactory<vector<typename Field::Element> > &factory) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix0 <Field, Row, Vector> Blackbox;

	commentator.start ("Testing random apply", "testRandomApply1", iterations);

	bool ret = true;
	bool iter_passed;

	unsigned long rank_A;

	VectorDomainBase<Field, Vector> VD (F);
	Vector w, lambda, mu, ATmu, x_correct, x_computed, d(r);

	VectorWrapper::ensureDim (lambda, factory.n ());
	VectorWrapper::ensureDim (mu, factory.n ());
	VectorWrapper::ensureDim (ATmu, m);
	VectorWrapper::ensureDim (x_correct, m);
	VectorWrapper::ensureDim (x_computed, m);
	VectorWrapper::ensureDim (w, n);

	int i, j, k, l;

	typename Field::Element x;

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);
		iter_passed = true;

		Blackbox *A = buildRandomSparseMatrix (F, n, m, r, K, d);
		MoorePenrose<Field, Vector> Adagger (F, A, r);
		
		{
			ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
			report << "Input matrix" << endl;
			A->prettyPrint (report, 4, 6);
		}

		Submatrix<Vector> Aprime (A, 0, 0, MIN (n, m), MIN (n, m));
		rank (rank_A, Aprime, F);

		if (rank_A == r) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
				<< "Rank is correct. Good." << endl;
		} else
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "Rank is incorrect (" << rank_A << "). Not good." << endl;

		while (factory) {
			factory.next (lambda);
			A->applyTranspose (x_correct, lambda);

			ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

			A->apply (w, x_correct);

			// Get a random element of N(A^T) and add it to w
			factory.next (mu);
			A->applyTranspose (ATmu, mu);

			for (j = 0; j < r; j++) {
				F.div (x, ATmu[j], d[j]);
				F.subin (x, mu[j]);
				F.addin (w[j], x);
			}

			report << "Right hand side: ";
			printVector<Field> (F, report, w);

			commentator.start ("Computing Moore-Penrose inverse");
			Adagger.apply (x_computed, w);
			commentator.stop ("done");

			report << "Correct output:  ";
			printVector<Field> (F, report, x_correct);

			report << "Computed output: ";
			printVector<Field> (F, report, x_computed);

			for (l = 0; l < x_correct.size (); l++)
				if (!F.areEqual (x_correct[l], x_computed[l]))
					ret = iter_passed = false;

			if (!iter_passed)
				commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
					<< "ERROR: Vectors are not equal" << endl;
		}


		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testRandomApply1");

	return ret;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 100;
	static size_t m = 100;
	static size_t r = 10;
	static integer q = 2147483647U;
	static int iterations = 100;
	static int k = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set row dimension of test matrices to N (default 100)",     TYPE_INT,     &n },
		{ 'm', "-m M", "Set column dimension of test matrices to M (default 100)",  TYPE_INT,     &m },
		{ 'r', "-r R", "Set rank of test matrices to R (default 10)",               TYPE_INT,     &r },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 100)",          TYPE_INT,     &iterations },
		{ 'k', "-k K", "Apply random Moore-Penrose to K vectors (default 1)",       TYPE_INT,     &k },
	};

	parseArguments (argc, argv, args);
	Modular<long> F (q);

	srand (time (NULL));

	cout << "MoorePenrose black box test suite" << endl << endl;

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

	RandomDenseVectorFactory<Modular<long> > factory1 (F, n, iterations);
	RandomDenseVectorFactory<Modular<long> > factory2 (F, n, k);

	if (!testIdentityApply<Modular<long> > (F, n, m, r, factory1)) pass = false;
	if (!testRandomApply1<Modular<long> > (F, n, m, r, iterations, 1.0 / (double) r, factory2)) pass = false;
#if 0
	if (!testRandomApply2<Modular<long> > (F, n, m, r, iterations, factory2)) pass = false;
#endif

	return pass ? 0 : -1;
}
