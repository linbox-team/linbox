
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
#include "linbox/field/PID-integer.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/direct-sum.h"
#include "linbox/solutions/rank.h"

#include "test-common.h"

using namespace LinBox;

// tests 1 and 2 were certain diagonals - now deemed unnecessary.  -bds 2005Mar15
/* Test 3: Rank of a random sparse matrix
 *
 * Constructs a random sparse matrix and computes its rank using Gaussian
 * elimination (direct and blas) and Wiedemann's algorithm. Checks that the results match.
 */

template <class Field>
bool testRankMethods(const Field &F, size_t n, unsigned int iterations, double sparsity = 0.05) 
{
	typedef SparseMatrix<Field,typename Vector<Field>::SparseSeq> Blackbox;

	commentator.start ("Testing elimination-based and blackbox rank", "testRankMethods", iterations);

	bool ret = true;
	unsigned int i;

	unsigned long rank_blackbox, rank_elimination, rank_hybrid;
	//unsigned long rank_Wiedemann, rank_elimination, rank_blas_elimination;

	typename Field::RandIter ri (F);

	for (i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		RandomSparseStream<Field, typename Vector<Field>::SparseSeq> stream (F, ri, sparsity, n, n);
		Blackbox A (F, stream);

		F.write( commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)) << endl; 
		A.write( commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION),FORMAT_MAPLE ) << endl; 

		rank (rank_blackbox, A, Method::Blackbox ());
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "blackbox rank " << rank_blackbox << endl;
		rank (rank_elimination, A, Method::Elimination());
		if (rank_blackbox != rank_elimination) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != elimination rank " << rank_elimination << endl;
			ret = false;
		}
		rank (rank_hybrid, A, Method::Hybrid());
		if (rank_blackbox != rank_hybrid) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != hybrid rank " << rank_hybrid << endl;
			ret = false;
		}
		
	/*
		rank (rank_Wiedemann, A, Method::Wiedemann ());
		//rank (rank_elimination, B, Method::SparseElimination());
		rank_elimination = rank_Wiedemann;
		rank (rank_blas_elimination, A, Method::BlasElimination ());

		commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Rank computed by Wiedemann: " << rank_Wiedemann << endl
			<< "Rank computed by sparse elimination: " << rank_elimination << endl
			<< "Rank computed by blas_elimination: " << rank_blas_elimination << endl;

		if (rank_Wiedemann != rank_elimination || rank_elimination != rank_blas_elimination) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Ranks are not equal" << endl;
			ret = false;
		}
		*/

		commentator.stop ("done");
		commentator.progress ();
	}

	commentator.stop (MSG_STATUS (ret), (const char *) 0, "testEliminationRank");

	return ret;
}

/* Test 4: Rank of zero and identity matrices by Wiedemann variants
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
				<< "ERROR: Wiedemann Rank of 0 is not 0, but is " << r << endl;
			ret = false;
		}

		F.init(one, 1);
		Blackbox I (F, n, one);
		rank (r, I, Method::Wiedemann ());
		if (r != n) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Wiedemann Rank of I is " << r << ", should be " << n << endl;
			ret = false;
		}

		DirectSum<Blackbox> B(A, I);
		rank (r, B, Method::Wiedemann ());
		if (r != n) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Wiedemann Rank of I+0 is " << r << ", should be " << n << endl;
			ret = false;
		}
                
                rank (r, B, Method::Wiedemann(Method::Wiedemann::SYMMETRIC));
		if (r != n) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Symmetric Wiedemann Rank of I+0 is " << r << ", should be " << n << endl;
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
	//static integer q = 1000003U;
	static int iterations = 2;
        static double sparsity = 0.05;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN (default 80)",       TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 65519)", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations (default 2)",           TYPE_INT,     &iterations },
                { 's', "-s S", "Sparse matrices with density S (default 0.05)",           TYPE_DOUBLE,     &sparsity },
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	cout << endl << "Black box rank test suite" << endl;
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	<< "over Modular<uint32>" << endl; 
	Modular<LinBox::uint32> F (q);
	if (!testRankMethods (F, n, iterations, sparsity)) pass = false;
	if (!testZeroAndIdentRank (F, n, 1)) pass = false;

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
	<< "over Modular<int>" << endl; 
	Modular<double> G (q);
    if (!testRankMethods (G, n, iterations, sparsity)) pass = false;
	if (!testZeroAndIdentRank (G, n, 1)) pass = false;


	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
	<< "over PID_integer" << endl; 
        PID_integer R;
	if (!testRankMethods (R, n, iterations, sparsity)) pass = false;

	return pass ? 0 : -1;
}
