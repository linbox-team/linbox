
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

#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/archetype.h"
#include "linbox/field/modular.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/blackbox/submatrix.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/moore-penrose.h"
#include "linbox/solutions/rank.h"
#include "linbox/vector/stream.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/* Build a random sparse n x m matrix of rank r with a nonsingular leading principal minor */

template <class Vector, class Field, class Row>
static SparseMatrix<Field,  Row>
*buildRandomSparseMatrix (Field                           &F,
			  size_t                           n,
			  size_t                           m,
			  size_t                           r,
			  double                           K,
			  vector<typename Field::Element> &dinv,
			  VectorStream<Row>               &top_right_stream,
			  VectorStream<Row>               &bottom_left_stream) 
{
	typedef SparseMatrix<Field, Row> Blackbox;

	Blackbox *A = new Blackbox (F, n, m);
	typename Field::RandIter rnd_p (F);
	NonzeroRandIter<Field> rnd (F, rnd_p);
	typename Field::Element factor;
	vector<Row> bottom_left_data (n - r);
	vector<Row> top_right_data (r);
	VectorDomain<Field> VD (F);

	// Build diagonal part
	for (size_t i = 0; i < r; i++) {
		typename Field::Element &x = A->refEntry (i, i);
		rnd.random (x);
		F.inv (dinv[i], x);
	}

	// Build top right part
	for (typename vector<Row>::iterator i = top_right_data.begin (); i != top_right_data.end (); i++) {
		top_right_stream.next (*i);
		VD.copy (A->getRow (top_right_stream.j () - 1), *i, r);
	}

	// Build bottom left part
	for (typename vector<Row>::iterator i = bottom_left_data.begin (); i != bottom_left_data.end (); i++) {
		bottom_left_stream.next (*i);
		VD.copy (A->getRow (r + bottom_left_stream.j () - 1), *i);
	}

	// Fill in bottom right part
	for (size_t i = 0; i < n - r; i++) {
		Row bottom_right_data;
		typename Row::first_type::iterator j_idx = bottom_left_data[i].first.begin ();
		typename Row::second_type::iterator j_elt = bottom_left_data[i].second.begin ();

		for (; j_idx != bottom_left_data[i].first.end (); ++j_idx, ++j_elt)
			VD.axpyin (bottom_right_data, F.mul (factor, *j_elt, dinv[*j_idx]),
				   top_right_data[*j_idx]);

		VD.copy (A->getRow (i + r), bottom_right_data, r);
	}

	top_right_stream.reset ();
	bottom_left_stream.reset ();

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
 * stream - Stream for random vectors to apply
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testIdentityApply (Field                                           &F,
			       size_t                                           n,
			       size_t                                           m,
			       size_t                                           r,
			       VectorStream<vector<typename Field::Element> > &stream) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef SparseMatrix <Field, Row> Blackbox;

	commentator.start ("Testing identity apply", "testIdentityApply", stream.m ());

	bool ret = true;
	bool iter_passed;

	Vector v, w;

	VectorWrapper::ensureDim (v, stream.n ());
	VectorWrapper::ensureDim (w, m);

	size_t i, l;

	Blackbox A (F, n, m);

	typename Field::Element x;

	F.init (x, 1);

	for (i = 0; i < r; i++)
		A.setEntry (i, i, x);

	MoorePenrose<Blackbox> Adagger (&A, r);

	while (stream) {
		commentator.startIteration (i);

		iter_passed = true;

		stream.next (v);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

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
 * stream - Stream for random vectors to apply
 *
 * Return true on success and false on failure
 */

template <class Field, class Vector, class Row>
static bool testRandomApply1 (Field                 &F,
			      size_t                 n,
			      size_t                 m,
			      size_t                 r,
			      unsigned               iterations,
			      double                 K,
			      VectorStream<Row>    &M_stream1,
			      VectorStream<Row>    &M_stream2,
			      VectorStream<Vector> &stream) 
{
	typedef SparseMatrix <Field, Row> Blackbox;

	commentator.start ("Testing random apply", "testRandomApply1", iterations);

	bool ret = true;
	bool iter_passed;

	unsigned long rank_A;

	Vector w, lambda, mu, ATmu, x_correct, x_computed;
	vector<typename Field::Element> dinv (r);

	VectorWrapper::ensureDim (lambda, stream.n ());
	VectorWrapper::ensureDim (mu, stream.n ());
	VectorWrapper::ensureDim (ATmu, m);
	VectorWrapper::ensureDim (x_correct, m);
	VectorWrapper::ensureDim (x_computed, m);
	VectorWrapper::ensureDim (w, n);

	size_t i, j, l;

	typename Field::Element x;

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);
		iter_passed = true;

		commentator.start ("Building requisite random sparse matrix");
		Blackbox *A = buildRandomSparseMatrix<Vector> (F, n, m, r, K, dinv, M_stream1, M_stream2);
		commentator.stop ("done");

		commentator.start ("Constructing Moore-Penrose inverse");
		MoorePenrose<Blackbox> Adagger (A, r);
		commentator.stop ("done");
		
		{
			ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
			report << "Input matrix" << endl;
			A->write (report, FORMAT_PRETTY);
		}

		Submatrix<Blackbox> Aprime (A, 0, 0, MIN (n, m), MIN (n, m));
		rank (rank_A, Aprime, Method::Wiedemann());

		if (rank_A == r) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION)
				<< "Rank is correct. Good." << endl;
		} else
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "Rank is incorrect (" << rank_A << "). Not good." << endl;

		while (stream) {
			stream.next (lambda);
			A->applyTranspose (x_correct, lambda);

			ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);

			A->apply (w, x_correct);

			// Get a random element of N(A^T) and add it to w
			stream.next (mu);
			A->applyTranspose (ATmu, mu);

			for (j = 0; j < r; j++) {
				F.mul (x, ATmu[j], dinv[j]);
				F.subin (x, mu[j]);
				F.addin (w[j], x);
			}

			report << "Right hand side: ";
			printVector<Field> (F, report, w);

			commentator.start ("Applying Moore-Penrose inverse");
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

		delete A;

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
		{ 'n', "-n N", "Set row dimension of test matrices to N.", TYPE_INT,     &n },
		{ 'm', "-m M", "Set column dimension of test matrices to M.", TYPE_INT,     &m },
		{ 'r', "-r R", "Set rank of test matrices to R.", TYPE_INT,     &r },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 'k', "-k K", "Apply random Moore-Penrose to K vectors.", TYPE_INT,     &k },
		{ '\0' }
	};

	typedef Modular<LinBox::uint32> Field;  //C.Pernet: avoids confusion with givaro::uint32
	typedef vector<Field::Element> DenseVector;
	typedef pair<vector<size_t>, vector<Field::Element> > SparseVector;

	parseArguments (argc, argv, args);
	Field F (q);

	srand (time (NULL));

	commentator.start("MoorePenrose black box test suite", "MoorePenrose");

	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_IMPORTANT);

	RandomDenseStream<Field, DenseVector> stream1 (F, n, iterations);
	RandomDenseStream<Field, DenseVector> stream2 (F, n, k);

	RandomSparseStream<Field, SparseVector> M_stream1 (F, 0.1, n - r, r);
	RandomSparseStream<Field, SparseVector> M_stream2 (F, 0.1, r, m - r);

	if (!testIdentityApply (F, n, m, r, stream1)) pass = false;
	if (!testRandomApply1 (F, n, m, r, iterations, 1.0 / (double) r, M_stream1, M_stream2, stream2)) pass = false;
#if 0
	if (!testRandomApply2 (F, n, m, r, iterations, stream2)) pass = false;
#endif

	commentator.stop("MoorePenrose black box test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
