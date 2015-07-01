/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-rank.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * -----------------------------------------------------
 *
 * This file is part of LinBox, licensed under the GNU Lesser General
 * Public License. See COPYING for more information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/solutions/rank.h"

#include "test-common.h"

using namespace LinBox;

/* Test 1: Rank of diagonal matrix with n/2 distinct nonzero entries
 *
 * Construct a random diagonal matrix with n/2 distinct nonzero entries and see
 * that its computed determinant equals the product of diagonal entries
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalRank1 (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field> Blackbox;

	commentator.start ("Testing diagonal rank (1)", "testDiagonalRank1", iterations);

	bool ret = true;
	bool done;
	int i;
	size_t j, k;

	Vector d(n);
	unsigned long _rank;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		for (j = 0; j < n / 2; j++) {
			do {
				do r.random (d[j]); while (F.isZero (d[j]));

				done = true;
				for (k = 0; k < j; k++) {
					if (F.areEqual (d[j], d[k])) {
						done = false;
						break;
					}
				}
			} while (!done);
		}

		for (j = n / 2; j < n; j++)
			F.init (d[j], 1);

// 		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		printVector<Field> (F, report, d);

		report << "True rank: " << n / 2 << endl;

		Blackbox D (F, d);

		rank (_rank, D, F, MethodTrait::Wiedemann ());

		report << "Computed rank: " << _rank << endl;

		if (_rank != n) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed rank is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalRank1");

	return ret;
}

/* Test 2: Rank of diagonal matrix with n / 2 nondistinct nonzero entries
 *
 * Construct a random diagonal matrix with n / 2 nonzero, nondistinct entries
 * and see that its computed rank equals n / 2
 *
 * F - Field over which to perform computations
 * n - Dimension to which to make matrix
 * iterations - Number of iterations to run
 *
 * Return true on success and false on failure
 */

template <class Field>
static bool testDiagonalRank2 (Field &F, size_t n, int iterations) 
{
	typedef vector <typename Field::Element> Vector;
	typedef vector <typename Field::Element> Polynomial;
	typedef vector <pair <size_t, typename Field::Element> > Row;
	typedef Diagonal <Field> Blackbox;

	commentator.start ("Testing diagonal rank (2)", "testDiagonalRank2", iterations);

	bool ret = true;
	int i;
	size_t j, k;

	Vector d(n);
	unsigned long _rank;
	typename Field::RandIter r (F);

	for (i = 0; i < iterations; i++) {
		commentator.startIteration (i);

		for (j = 0; j < n / 4; j++)
			do r.random (d[j]); while (F.isZero (d[j]));

		for (j = n / 4; j < n / 2; j++) {
			k = rand () % (n / 4);
			d[j] = d[k];
		}

		for (j = n / 2; j < n; j++)
			F.init (d[j], 0);

		ostream &report = commentator.report (Commentator::LEVEL_UNIMPORTANT, INTERNAL_DESCRIPTION);
		report << "Diagonal entries: ";
		printVector<Field> (F, report, d);

		report << "True rank: " << n / 2 << endl;

		Blackbox D (F, d);

		rank (_rank, D, F, MethodTrait::Wiedemann ());

		report << "Computed rank: " << _rank << endl;

		if (_rank != n / 2) {
			ret = false;
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Computed rank is incorrect" << endl;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testDiagonalRank2");

	return ret;
}

/* Test 3: Rank of a random sparse matrix
 *
 * Constructs a random sparse matrix and computes its rank using both Gaussian
 * elimination and Wiedemann's algorithm. Checks that the results match
 */

template <class Field>
bool testEliminationRank (const Field &F, size_t n, unsigned int iterations) 
{
	typedef SparseMatrix<Field,typename Vector<Field>::SparseSeq> Blackbox;

	commentator.start ("Testing elimination-based rank", "testEliminationRank", iterations);

	bool ret = true;
	unsigned int i;

	unsigned long rank_Wiedemann, rank_elimination;

	typename Field::RandIter ri (F);

	for (i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		RandomSparseStream<Field, typename Vector<Field>::SparseSeq> stream (F, ri, 0.005, n, n);
		Blackbox A (F, stream);

		rank (rank_Wiedemann, A, F, MethodTrait::Wiedemann ());
		rank (rank_elimination, A, F, MethodTrait::Elimination (EliminationTraits::PIVOT_LINEAR));

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Rank computed by Wiedemann: " << rank_Wiedemann << endl
			<< "Rank computed by elimination: " << rank_elimination << endl;

		if (rank_Wiedemann != rank_elimination) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Ranks are not equal" << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testEliminationRank");

	return ret;
}

/* Test 4: Rank of zero and identity matrices
 *
 */

template <class Field>
bool testZeroAndIdentRank (const Field &F, size_t n, unsigned int iterations) 
{
	typedef ScalarMatrix<Field> Blackbox;

	commentator.start ("Testing rank of 0 and I", "testZeroAndIdentRank", iterations);

	bool ret = true;
	unsigned int i;

	unsigned long r; // rank

	for (i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		typename Field::Element zero, one;
		F.init(zero, 0);
		F.init(one, 1);

		Blackbox A (F, n, zero);
		Blackbox I (F, n, one);

		rank (r, A, F, MethodTrait::Wiedemann ());
		if (r != 0) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Rank of 0 is not 0 but " << r << endl;
			ret = false;
		}

		rank (r, I, F, MethodTrait::Wiedemann ());
		if (r != n) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Rank of I is not " << n << " but " << r << endl;
			ret = false;
		}

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testZeroAndIdentRank");

	return ret;
}

int main (int argc, char **argv)
{

//     commentator.setMaxDetailLevel( 100000 );
//     commentator.setMaxDepth( 100000 );
   
	bool pass = true;

	static size_t n = 256;
	static integer q = 2147483647U;
	static int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 256)",       TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<uint32> F (q);

	srand (time (NULL));

	cout << endl << "Black box rank test suite" << endl;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	if (!testDiagonalRank1 (F, n, iterations)) pass = false;
	if (!testDiagonalRank2 (F, n, iterations)) pass = false;
	if (!testEliminationRank (F, n, iterations)) pass = false;
	if (!testZeroAndIdentRank (F, n, 1)) pass = false;

	return pass ? 0 : -1;
}
