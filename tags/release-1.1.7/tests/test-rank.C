

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

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include "linbox/field/givaro-zpz.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/direct-sum.h"
#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/gauss-gf2.h"
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

bool testRankMethodsGF2(const GF2& F2, size_t n, unsigned int iterations, double sparsity = 0.05) 
{
	typedef ZeroOne<GF2> Blackbox;
	typedef SparseMatrix<Modular<double>,Vector<Modular<double> >::SparseSeq> MdBlackbox;
	Modular<double> MdF2(2);
	GF2::Element one; Modular<double>::Element mdone;
	F2.init(one,true);
	MdF2.init(mdone,1UL);


	commentator.start ("Testing elimination-based and blackbox rank over GF2", "testRankMethodsGF2", iterations);

	bool ret = true;
	unsigned int i;

	unsigned long rank_blackbox, rank_elimination, rank_sparselimination, rank_sparse;
	//unsigned long rank_Wiedemann, rank_elimination, rank_blas_elimination;

	GF2::RandIter ri (F2);

	for (i = 0; i < iterations; ++i) {
		commentator.startIteration (i);

		Blackbox A(F2,n,n);
		MdBlackbox B(MdF2,n,n);
		for(size_t ii=0; ii<n;++ii) {
			for(size_t jj=0; jj<n; ++jj) {
				if (drand48()<sparsity) {
					A.setEntry(ii,jj,one);
					B.setEntry(ii,jj,mdone);
				}
			}
		}

		F2.write( commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)) << endl; 
		B.write( commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION),FORMAT_GUILLAUME ) << endl; 
		A.write( commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION),FORMAT_GUILLAUME ) << endl; 


		rank (rank_blackbox, A, Method::Blackbox ());
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "blackbox rank " << rank_blackbox << endl;

		rank (rank_elimination, B, Method::BlasElimination());
		if (rank_blackbox != rank_elimination) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != BLAS elimination rank " << rank_elimination << endl;
			ret = false;
		}

		rankin (rank_sparselimination, A, Method::SparseElimination());
		if (rank_blackbox != rank_sparselimination) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != sparse elimination GF2 rank " << rank_elimination << endl;
			ret = false;
		}

		
		rankin (rank_sparse, B, Method::SparseElimination());

		if (rank_sparselimination != rank_sparse) {
			commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rank sparse elimination GF2 != sparse rank " << rank_sparse << endl;
			ret = false;
		}
		
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

	static size_t n = 40;
	static integer q = 65519U;
	//static integer q = 1000003U;
	static integer bigQ("12345678901234567890123456789012345678901234568119");
	static int iterations = 1;
        static double sparsity = 0.05;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
        { 's', "-s S", "Sparse matrices with density S.", TYPE_DOUBLE,     &sparsity },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	commentator.start("rank solution test suite", "rank");
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

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
	<< "over GivaroZpz<Integer>" << endl; 
        GivaroZpz<Integer> Gq(bigQ);
	if (!testRankMethods (Gq, n, iterations, sparsity)) pass = false;
	if (!testZeroAndIdentRank (Gq, n, 1)) pass = false;

	commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) 
	<< "over GF2" << endl; 
        GF2 F2;
	if (!testRankMethodsGF2 (F2, n, iterations, sparsity)) pass = false;


	commentator.stop("rank solution test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
