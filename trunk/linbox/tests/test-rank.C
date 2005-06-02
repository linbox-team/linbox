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
#include "linbox/blackbox/direct-sum.h"
#include "linbox/solutions/rank.h"

#include "test-common.h"

using namespace LinBox;

// tests 1 and 2 were certain diagonals - now deemed unnecessary.  -bds 2005Mar15 (witnesses: ZW and PG)
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

	unsigned long rank_Wiedemann, rank_elimination, rank_blas_elimination;

	typename Field::RandIter ri (F);

	for (i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		RandomSparseStream<Field, typename Vector<Field>::SparseSeq> stream (F, ri, 0.05, n, n);
		Blackbox A (F, stream);

		A.write( commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)) << endl; 

		rank (rank_Wiedemann, A, Method::Wiedemann ());
		rank (rank_elimination, A, Method::SparseElimination(SparseEliminationTraits::PIVOT_LINEAR));
		rank (rank_blas_elimination, A, Method::BlasElimination ());

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Rank computed by Wiedemann: " << rank_Wiedemann << endl
			<< "Rank computed by elimination: " << rank_elimination << endl
			<< "Rank computed by blas_elimination: " << rank_blas_elimination << endl;

		if (rank_Wiedemann != rank_elimination || rank_elimination != rank_blas_elimination) {
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

	commentator.start ("Testing rank of zero and Identity and half/half matrices", "testZeroAndIdentRank", iterations);

	bool ret = true;
	unsigned int i;

	unsigned long r; // rank

	for (i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		typename Field::Element zero, one;

		F.init(zero, 0);
		Blackbox A (F, n, zero);
		rank (r, A, Method::Wiedemann ());
		if (r != 0) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Rank of 0 is not 0, but is " << r << endl;
			ret = false;
		}

		F.init(one, 1);
		Blackbox I (F, n, one);
		rank (r, I, Method::Wiedemann ());
		if (r != n) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Rank of I is not " << n << ", but is " << r << endl;
			ret = false;
		}

		DirectSum<Blackbox> B(A, I);
		rank (r, B, Method::Wiedemann ());
		if (r != n) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Rank of I+0 is not " << n << ", but is " << r << endl;
			ret = false;
		}
                
                rank (r, B, Method::Wiedemann(Method::Wiedemann::SYMMETRIC));
		if (r != n) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Symmetric Rank of I+0 is not " << n << ", but is " << r << endl;
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

	static size_t n = 80;
	static integer q = 65519U;
	static int iterations = 2;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 256)",       TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 2147483647)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
	};

	parseArguments (argc, argv, args);
	Modular<LinBox::uint32> F (q);
	//Modular<int> F (q);

	srand (time (NULL));

	cout << endl << "Black box rank test suite" << endl;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	if (!testEliminationRank (F, n, iterations)) pass = false;
	if (!testZeroAndIdentRank (F, n, 1)) pass = false;
	//Modular<int> G (3);
	//if (!testEliminationRank (G, n, iterations)) pass = false;
	//if (!testZeroAndIdentRank (G, n, 1)) pass = false;

	return pass ? 0 : -1;
}
