

/* tests/test-rank.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * -----------------------------------------------------
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
 */


/*! @file  tests/test-rank.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <cstdio>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include "linbox/field/gf2.h"
#include "linbox/field/givaro.h"
#include "linbox/blackbox/diagonal.h"
// #include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/sparse-matrix.h"
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
	 // typedef SparseMatrix<Field,SparseMatrixFormat::SparseSeq > Blackbox;
	 // typedef SparseMatrix<Field,SparseMatrixFormat::SparsePar > Blackbox;
	 // typedef SparseMatrix<Field,SparseMatrixFormat::SparseMap > Blackbox;
	typedef SparseMatrix<Field,SparseMatrixFormat::COO> Blackbox;
	// typedef SparseMatrix<Field,SparseMatrixFormat::CSR> Blackbox; // inf loop
	// typedef SparseMatrix<Field,SparseMatrixFormat::ELL> Blackbox;
	// typedef SparseMatrix<Field,SparseMatrixFormat::ELL_R> Blackbox;
	// typedef SparseMatrix<Field,SparseMatrixFormat::HYB> Blackbox;
	// typedef SparseMatrix<Field,SparseMatrixFormat::TPL> Blackbox;

	commentator().start ("Testing elimination-based and blackbox rank", "testRankMethods", (unsigned int)iterations);

	bool ret = true;
	unsigned int i;

	unsigned long rank_blackbox, rank_elimination, rank_hybrid;
	//unsigned long rank_Wiedemann, rank_elimination, rank_blas_elimination;

	typename Field::RandIter ri (F);

	for (i = 0; i < iterations; ++i) {
		commentator().startIteration (i);

		RandomSparseStream<Field, typename Blackbox::Row> stream (F, ri, sparsity, n, n);
		// RandomSparseStream<Field, typename Vector<Field>::SparseSeq> stream (F, ri, sparsity, n, n);
		Blackbox A (F, stream);
		// std::cout << A.rowdim() << ',' << A.coldim() << std::endl;

		F.write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)) << endl;
		A.write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION),Tag::FileFormat::Maple ) << endl;

		LinBox::rank (rank_blackbox, A, Method::Blackbox ());
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "blackbox rank " << rank_blackbox << endl;
		LinBox::rank (rank_elimination, A, Method::Elimination());
		if (rank_blackbox != rank_elimination) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != elimination rank " << rank_elimination << endl;
			ret = false;
		}
		LinBox::rank (rank_hybrid, A, Method::Hybrid());
		if (rank_blackbox != rank_hybrid) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != hybrid rank " << rank_hybrid << endl;
			ret = false;
		}

		{
			size_t rank_Wiedemann ;
			LinBox::rank (rank_Wiedemann, A, Method::Wiedemann ());
			rank_elimination = rank_Wiedemann;
			size_t rank_blas_elimination ;
			if (F.characteristic() < LinBox::BlasBound && F.characteristic() == F.cardinality()) {
				LinBox::rank (rank_blas_elimination, A, Method::BlasElimination ());
			}
			else {
				LinBox::rank (rank_blas_elimination, A, Method::Elimination ());
			}

			commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< "Rank computed by Wiedemann: " << rank_Wiedemann << endl
			<< "Rank computed by sparse elimination: " << rank_elimination << endl
			<< "Rank computed by blas_elimination: " << rank_blas_elimination << endl;

			if (rank_Wiedemann != rank_elimination || rank_elimination != rank_blas_elimination) {
				commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Ranks are not equal" << endl;
				ret = false;
			}
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testEliminationRank");

	return ret;
}

// this test just doesn't work/compile
#if 0
bool testRankMethodsGF2(const GF2& F2, size_t n, unsigned int iterations, double sparsity = 0.05)
{
	typedef ZeroOne<GF2> Blackbox;
	typedef SparseMatrix<Modular<double>,Vector<Modular<double> >::SparseSeq> MdBlackbox;
	Modular<double> MdF2(2);
	GF2::Element one; Modular<double>::Element mdone;
	MdF2.assign(mdone,MdF2.one);


	commentator().start ("Testing elimination-based and blackbox rank over GF2", "testRankMethodsGF2", (unsigned int)iterations);

	bool ret = true;
	unsigned int i;

	unsigned long rank_blackbox, rank_elimination, rank_sparselimination, rank_sparse;
	//unsigned long rank_Wiedemann, rank_elimination, rank_blas_elimination;

	GF2::RandIter ri (F2);

	for (i = 0; i < iterations; ++i) {
		commentator().startIteration (i);

		Blackbox A(F2,n,n);
		MdBlackbox B(MdF2,n,n);
		for(size_t ii=0; ii<n;++ii) {
			for(size_t jj=0; jj<n; ++jj) {
				if (drand48()<sparsity) {
					A.setEntry(ii,jj,F2.one);
					B.setEntry(ii,jj,mdone);
				}
			}
		}

		F2.write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)) << endl;
		B.write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION),Tag::FileFormat::Guillaume ) << endl;
		A.write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION),Tag::FileFormat::Guillaume ) << endl;


		LinBox::rank (rank_blackbox, A, Method::Blackbox ());
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "blackbox rank " << rank_blackbox << endl;

			LinBox::rank (rank_elimination, B, Method::BlasElimination());
		if (rank_blackbox != rank_elimination) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != BLAS elimination rank " << rank_elimination << endl;
			ret = false;
		}

		rankin (rank_sparselimination, A, Method::SparseElimination());
		if (rank_blackbox != rank_sparselimination) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != sparse elimination GF2 rank " << rank_elimination << endl;
			ret = false;
		}


		rankin (rank_sparse, B, Method::SparseElimination());

		if (rank_sparselimination != rank_sparse) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: rank sparse elimination GF2 != sparse rank " << rank_sparse << endl;
			ret = false;
		}

		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testEliminationRank");

	return ret;
}
#endif

/* Test 4: Rank of zero and identity matrices by Wiedemann variants
 *
 */
template <class Field>
bool testZeroAndIdentRank (const Field &F, size_t n, unsigned int iterations)
{
	typedef ScalarMatrix<Field> Blackbox;

	commentator().start ("Testing rank of zero and Identity and half/half matrices", "testZeroAndIdentRank", (unsigned int)iterations);

	bool ret = true;
	unsigned int i;

	unsigned long r; // rank

	for (i = 0; i < iterations; ++i) {
		commentator().startIteration (i);


		Blackbox A (F, n, n, F.zero);
		LinBox::rank (r, A, Method::Wiedemann ());
		if (r != 0) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Wiedemann Rank of 0 is not 0, but is " << r << endl;
			ret = false;
		}

		Blackbox I (F, n, n, F.one);
		LinBox::rank (r, I, Method::Wiedemann ());
		if (r != n) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Wiedemann Rank of I is " << r << ", should be " << n << endl;
			ret = false;
		}

		DirectSum<Blackbox> B(A, I);
		LinBox::rank (r, B, Method::Wiedemann ());
		if (r != n) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Wiedemann Rank of I+0 is " << r << ", should be " << n << endl;
			ret = false;
		}

		LinBox::rank (r, B, Method::Wiedemann(Method::Wiedemann::SYMMETRIC));
		if (r != n) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Symmetric Wiedemann Rank of I+0 is " << r << ", should be " << n << endl;
			ret = false;
		}
		commentator().stop ("done");
		commentator().progress ();
	}

	commentator().stop (MSG_STATUS (ret), (const char *) 0, "testZeroAndIdentRank");

	return ret;
}

int main (int argc, char **argv)
{

//     commentator().setMaxDetailLevel( 100000 );
//     commentator().setMaxDepth( 100000 );

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
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	srand ((unsigned)time (NULL));
	// srand48 ((unsigned)time (NULL));

	commentator().start("rank solution test suite", "rank");
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	<< "over Modular<uint32_t>" << endl;
	Modular<uint32_t> F (q);
	if (!testRankMethods (F, n, (unsigned int)iterations, sparsity)) pass = false;
	if (!testZeroAndIdentRank (F, n, 1)) pass = false;

	commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	<< "over Modular<int>" << endl;
	Modular<double> G (q);
    	if (!testRankMethods (G, n, (unsigned int)iterations, sparsity)) pass = false;
	if (!testZeroAndIdentRank (G, n, 1)) pass = false;

	commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	<< "over PID_integer" << endl;
	PID_integer R;
	if (!testRankMethods (R, n, (unsigned int)iterations, sparsity)) pass = false;

	commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	<< "over GivaroZpz<Integer>" << endl;
        GivaroZpz<Integer> Gq(bigQ);
	if (!testRankMethods (Gq, n, (unsigned int)iterations, sparsity)) pass = false;
	if (!testZeroAndIdentRank (Gq, n, 1)) pass = false;

#if 0
	commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
	<< "over GF2" << endl;
        GF2 F2;
	if (!testRankMethodsGF2 (F2, n, (unsigned int)iterations, sparsity)) pass = false;
#endif


	commentator().stop("rank solution test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
