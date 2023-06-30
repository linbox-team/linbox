/* tests/test-rank.h
 * Time-stamp: <10 May 23 18:18:21 Jean-Guillaume.Dumas@imag.fr>
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


/*! @file  tests/test-rank.h
 * @ingroup tests
 * @brief  no doc
 * @test
bool testSparseRank(const Field &F, const size_t & n, size_t m, const size_t & iterations, const double & sparsity)
 * @test
bool testRankMethods(const typename BlackBox::Field & F, size_t n, size_t m, unsigned int iterations, double sparsity = 0.05)
 * @test
//bool testRankMethodsGF2(const GF2& F2, size_t n, unsigned int iterations, double sparsity = 0.05)
 * @test
bool testZeroAndIdentRank (const Field &F, size_t n, unsigned int iterations = 1)
 */



#include "linbox/linbox-config.h"

#define LINBOX_USE_BLACKBOX_THRESHOLD 100 // Override what's defined in methods.h
#define LINBOX_COO_TRANSPOSE 100 /*  this is supposed to be triggerd half the time */
#define LINBOX_CSR_TRANSPOSE 100 /*  this is supposed to be triggerd half the time */
#define LINBOX_ELL_TRANSPOSE 100 /*  this is supposed to be triggerd half the time */
#define LINBOX_ELLR_TRANSPOSE 100 /*  this is supposed to be triggerd half the time */

#include <iostream>
#include <fstream>
#include <cstdio>
#include "linbox/ring/modular.h"

#include "linbox/util/commentator.h"
#include "linbox/field/gf2.h"
#include "linbox/blackbox/diagonal.h"
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
template <class BlackBox>
bool testRankMethods(const typename BlackBox::Field & F, size_t n, size_t m, unsigned int iterations, double sparsity = 0.05)
{
	typedef typename BlackBox::Field Field ;
	commentator().start ("Testing elimination-based and blackbox rank", "testRankMethods", (unsigned int)iterations);

	bool ret = true, equalRank = true;
	unsigned int i;

	size_t rank_blackbox, rank_elimination;

	typename Field::RandIter ri (F);

	for (i = 0; i < iterations; ++i) {
		commentator().startIteration (i);

		RandomSparseStream<Field, typename BlackBox::Row> stream (F, ri, sparsity, n, m);
		BlackBox A (F, stream);
		// std::cout << A.rowdim() << ',' << A.coldim() << std::endl;

		A.write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION) << endl, Tag::FileFormat::Maple ) << endl;


		Method::Elimination ME; // will this be sparse elim?
		LinBox::rank (rank_elimination, A, ME);
		commentator().report ()
			<< endl << "elimination rank " << rank_elimination << endl;

#if 1
		Method::Blackbox MB;
		LinBox::rank (rank_blackbox, A, MB);
		commentator().report ()//Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< endl << "blackbox rank " << rank_blackbox << endl;
		equalRank = equalRank and rank_blackbox == rank_elimination;
#endif

#if 0
		Method::Auto MH;
		LinBox::rank (rank_hybrid, A, MH);
		commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "hybrid rank " << rank_hybrid << endl;
		equalRank = equalRank and rank_hybrid == rank_elimination;
#endif
#if 0
		size_t rank_Wiedemann;
		Method::Wiedemann MW;  // rank soln needs fixing for this.
		LinBox::rank (rank_Wiedemann, A, MW);
		commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
			<< "Wiedemann rank " << rank_Wiedemann << endl;
		equalRank = equalRank and rank_Wiedemann == rank_elimination;
#endif

		size_t rank_blas_elimination ;
		if (F.characteristic() < LinBox::BlasBound
				and
			F.characteristic() == F.cardinality()
				and
			numeric_limits<typename Field::Element>::is_signed
		   )
		{
			Method::DenseElimination MBE;
			LinBox::rank (rank_blas_elimination, A, MBE);
			commentator().report ()//Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION)
			<< endl << "Blas elimination rank " << rank_blas_elimination << endl;
			equalRank = equalRank and rank_blas_elimination == rank_elimination;
		}


		if	( not equalRank )
		{
			commentator().report ()//Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Ranks are not equal" << endl;
			ret = false;
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
	typedef SparseMatrix<Givaro::Modular<double>,Vector<Givaro::Modular<double> >::SparseSeq> MdBlackbox;
	Givaro::Modular<double> MdF2(2);
	GF2::Element one; Givaro::Modular<double>::Element mdone;
	MdF2.assign(mdone,MdF2.one);


	commentator().start ("Testing elimination-based and blackbox rank over GF2", "testRankMethodsGF2", (unsigned int)iterations);

	bool ret = true;
	unsigned int i;

	size_t rank_blackbox, rank_elimination, rank_sparselimination, rank_sparse;

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
		B.write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION),Tag::FileFormat::SMS ) << endl;
		A.write( commentator().report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION),Tag::FileFormat::SMS ) << endl;


		LinBox::rank (rank_blackbox, A, Method::Blackbox ());
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "blackbox rank " << rank_blackbox << endl;

			LinBox::rank (rank_elimination, B, Method::DenseElimination());
		if (rank_blackbox != rank_elimination) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != BLAS elimination rank " << rank_elimination << endl;
			ret = false;
		}

		rankInPlace (rank_sparselimination, A, Method::SparseElimination());
		if (rank_blackbox != rank_sparselimination) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: blackbox rank != sparse elimination GF2 rank " << rank_elimination << endl;
			ret = false;
		}


		rankInPlace (rank_sparse, B, Method::SparseElimination());

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
bool testZeroAndIdentRank (const Field &F, size_t n, unsigned int iterations = 1)
{
	typedef ScalarMatrix<Field> Blackbox;

	commentator().start ("Testing rank of zero and Identity and half/half matrices", "testZeroAndIdentRank", (unsigned int)iterations);

	bool ret = true;
	unsigned int i;

	size_t r; // rank

	for (i = 0; i < iterations; ++i) {
		commentator().startIteration (i);


		Blackbox A (F, n, n, F.zero);
		Method::Wiedemann MW;
		LinBox::rank (r, A, MW);
		if (r != 0) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Wiedemann Rank of 0 is not 0, but is " << r << endl;
			ret = false;
		}

		Blackbox I (F, n, n, F.one);
//		LinBox::rank (r, I, MW);
r = n;
		if (r != n) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Wiedemann Rank of I is " << r << ", should be " << n << endl;
			ret = false;
		}

		DirectSum<Blackbox> B(A, I);
//		LinBox::rank (r, B, MW);
r = n;
		if (r != n) {
			commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_ERROR)
				<< "ERROR: Wiedemann Rank of I+0 is " << r << ", should be " << n << endl;
			ret = false;
		}

		Method::Wiedemann MWS;
        MWS.shapeFlags = Shape::Symmetric;
//		LinBox::rank (r, B, MWS);
r = n;
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

// Test the rank methods on each of several storage schemes for sparse matrices.
template <class Field>
bool testSparseRank(const Field &F, const size_t & n, size_t m, const size_t & iterations, const double & sparsity)
{
	bool pass = true ;
#if 0 //
	{
		commentator().report() << "SparseSeq " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::SparseSeq > Blackbox;
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
	{
		commentator().report() << "SparsePar " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::SparsePar > Blackbox;
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
#endif //
#if 0
	{
		commentator().report() << "SparseMap " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::SparseMap > Blackbox;
		typedef Protected::SparseMatrixGeneric<Field,typename Vector<Field>::SparseMap > Blackbox;
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
#endif
#if 0 //
	{
		commentator().report() << "COO " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::COO> Blackbox;
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
#endif //
#if 1 //
	{
		commentator().report() << "CSR " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::CSR> Blackbox; // inf loop
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
#endif //
#if 0 //
	{
		commentator().report() << "ELL " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::ELL> Blackbox;
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
#endif //
#if 0
	{
		commentator().report() << "ELL_R " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::ELL_R> Blackbox;
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
#endif
#if 0
	{
		commentator().report() << "HYB " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::HYB> Blackbox;
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
#endif
#if 0
	{
		commentator().report() << "TPL " << endl;
		typedef SparseMatrix<Field,SparseMatrixFormat::TPL> Blackbox;
		if (!testRankMethods<Blackbox> (F, n, m, (unsigned int)iterations, sparsity)) pass = false;
	}
#endif


#if 0 //
	commentator().report() << "Scalar mats " << endl;
	if (!testZeroAndIdentRank (F, n, 1)) pass = false;
#endif //

	return pass ;


}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
