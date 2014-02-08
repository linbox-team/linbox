/* tests/test-sparse2.C
 * Copyright (C) 2014 the LinBox group
 *
 * Written by bds <saunders@udel.edu>
 *            BB <bbboyer@ncsu.edu>
 *
 * --------------------------------------------------------
 *
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
 *
 */

/*! @file  tests/test-sparse.C
 * @ingroup tests
 *
 * @brief no doc
 *
 * @test no doc.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>


#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/sparse-matrix.h"

#include "test-blackbox.h"

using namespace LinBox;

template <class SM>
void buildBySetEntry(SM & A, size_t nnz)
{
	typename SM::Field::RandIter r(A.field(),0);
	srand(0);

	//size_t i, j, k;
	typename SM::Field::Element x;

	for (size_t k = 0; k < nnz; ++k)
	{
		size_t i = rand() % A.rowdim();
		size_t j = rand() % A.coldim();
		r.nonzerorandom(x);
		// r.random(x); // I want to see what happens when one reads in zero. If we don't want it, we just stop permitting setting zero... (hence the clearEntry function)
		// std::cout << "A[ " << i+1 << ',' << j+1 << "]:=" << x << ';' << std::endl;
		// if (A.field().isZero(x)) std::cout << "#is zero" << std::endl;
		// std::cout << i << ',' << j << ',' << x << std::endl;
		A.setEntry(i,j,x);
	}
	// std::cout << "---" << std::endl;
	A.finalize();
	// std::cout << "B :=" << std::endl;
	// A.write(std::cout);
	// std::cout << ';' << std::endl;
}

template <class SM>
bool buildBySetGetEntry(SM & A, const SM &B)
{

	A.resize(B.rowdim(),B.coldim());

	for (size_t i = 0; i < B.rowdim(); ++i)
		for (size_t j = 0; j < B.coldim(); ++j)
		{
			A.setEntry(i,j,B.getEntry(i,j));
		}
	A.finalize();

	MatrixDomain<typename SM::Field> MD(B.field()) ;

	return MD.areEqual(A,B);
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static size_t m = 10;
	static integer q = 11;
	static size_t N = m+n;

	static Argument args[] = {
		{ 'n', "-n N", "Set column dimension of test matrices to N.", TYPE_INT,     &n },
		{ 'm', "-m M", "Set row dimension of test matrices to M.", TYPE_INT,     &m },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'N', "-N N", "N nonzero Elements in sparse matrices.", TYPE_INT,     &N },
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	//typedef	Modular<uint32_t> Field;
	typedef	Modular<double> Field;
	// typedef Field::Element  Element;

	Field F (q);

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator().start("Sparse matrix black box test suite", "Sparse");
	MatrixDomain<Field> MD(F) ;

	 /*  default */
	commentator().start("SparseMatrix<Field>", "Field");
	SparseMatrix<Field> S1(F, m, n);
	buildBySetEntry(S1, N);
	// it is assumed that any other matrix is built with the same random number generator  with the same seed
	if ( testBlackbox(S1,false))
		commentator().stop("SparseMatrix<Field> pass");
	else {
		commentator().stop("SparseMatrix<Field> FAIL");
		pass = false;
	}

	{ /*  COO */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::COO>", "COO");
		SparseMatrix<Field, SparseMatrixFormat::COO> S2(F, m, n);
		SparseMatrix<Field, SparseMatrixFormat::COO> S3(F, m, n);
		buildBySetEntry(S2, N);
		if ( ! testBlackbox(S2,true) ) {
			commentator().stop("Format COO FAIL (testBlackbox)");
			pass = false;
		}
		else if ( ! MD.areEqual(S1,S2) ) {
			commentator().stop("Format COO FAIL (areEqual)");
			pass = false;
		}
		else if (! buildBySetGetEntry(S3,S2) ) {
			commentator().stop("Format COO FAIL (getEntry)");
			pass = false;
		}
		else {
			commentator().stop("Format COO pass");
		}
	}

	{ /*  CSR */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::CSR>", "CSR");
		SparseMatrix<Field, SparseMatrixFormat::CSR> S2(F, m, n);
		SparseMatrix<Field, SparseMatrixFormat::CSR> S3(F, m, n);
		buildBySetEntry(S2, N);
		if ( ! testBlackbox(S2,true) ) {
			commentator().stop("Format CSR FAIL (testBlackbox)");
			pass = false;
		}
		else if ( ! MD.areEqual(S1,S2) ) {
			commentator().stop("Format CSR FAIL (areEqual)");
			pass = false;
		}
		else if (! buildBySetGetEntry(S3,S2) ) {
			commentator().stop("Format CSR FAIL (getEntry)");
			pass = false;
		}
		else {
			commentator().stop("Format CSR pass");
		}
	}

	{ /*  ELL */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::ELL>", "ELL");
		SparseMatrix<Field, SparseMatrixFormat::ELL> S2(F, m, n);
		SparseMatrix<Field, SparseMatrixFormat::ELL> S3(F, m, n);
		buildBySetEntry(S2, N);
		if ( ! testBlackbox(S2,true) ) {
			commentator().stop("Format ELL FAIL (testBlackbox)");
			pass = false;
		}
		else if ( ! MD.areEqual(S1,S2) ) {
			commentator().stop("Format ELL FAIL (areEqual)");
			pass = false;
		}
		else if (! buildBySetGetEntry(S3,S2) ) {
			commentator().stop("Format ELL FAIL (getEntry)");
			pass = false;
		}
		else {
			commentator().stop("Format ELL pass");
		}
	}

	{ /*  ELL_R */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::ELL_R>", "ELL_R");
		SparseMatrix<Field, SparseMatrixFormat::ELL_R> S2(F, m, n);
		SparseMatrix<Field, SparseMatrixFormat::ELL_R> S3(F, m, n);
		buildBySetEntry(S2, N);
		if ( ! testBlackbox(S2,true) ) {
			commentator().stop("Format ELL_R FAIL (testBlackbox)");
			pass = false;
		}
		else if ( ! MD.areEqual(S1,S2) ) {
			commentator().stop("Format ELL_R FAIL (areEqual)");
			pass = false;
		}
		else if (! buildBySetGetEntry(S3,S2) ) {
			commentator().stop("Format ELL_R FAIL (getEntry)");
			pass = false;
		}
		else {
			commentator().stop("Format ELL_R pass");
		}
	}

#if 0 // doesn't compile
	commentator().start("SparseMatrix<Field, SparseMatrixFormat::HYB>", "HYB");
	SparseMatrix<Field, SparseMatrixFormat::HYB> S6(F, m, n);
	buildBySetEntry(S6, N);
	S6.optimise();
	if ( testBlackbox(S6,false)  && MD.areEqual(S1,S6))
		commentator().stop("Format HYB pass");
	else {
		commentator().stop("Format HYB FAIL");
		pass = false;
	}
#endif

	{ /*  SparseSeq */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::SparseSeq>", "SparseSeq");
		SparseMatrix<Field, SparseMatrixFormat::SparseSeq> S2(F, m, n);
		SparseMatrix<Field, SparseMatrixFormat::SparseSeq> S3(F, m, n);
		buildBySetEntry(S2, N);
		if ( ! testBlackbox(S2,true) ) {
			commentator().stop("Format SparseSeq FAIL (testBlackbox)");
			pass = false;
		}
		else if ( ! MD.areEqual(S1,S2) ) {
			commentator().stop("Format SparseSeq FAIL (areEqual)");
			pass = false;
		}
		else if (! buildBySetGetEntry(S3,S2) ) {
			commentator().stop("Format SparseSeq FAIL (getEntry)");
			pass = false;
		}
		else {
			commentator().stop("Format SparseSeq pass");
		}
	}

	{ /*  SparsePar */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::SparsePar>", "SparsePar");
		SparseMatrix<Field, SparseMatrixFormat::SparsePar> S2(F, m, n);
		SparseMatrix<Field, SparseMatrixFormat::SparsePar> S3(F, m, n);
		buildBySetEntry(S2, N);
		if ( ! testBlackbox(S2,true) ) {
			commentator().stop("Format SparsePar FAIL (testBlackbox)");
			pass = false;
		}
		else if ( ! MD.areEqual(S1,S2) ) {
			commentator().stop("Format SparsePar FAIL (areEqual)");
			pass = false;
		}
		else if (! buildBySetGetEntry(S3,S2) ) {
			commentator().stop("Format SparsePar FAIL (getEntry)");
			pass = false;
		}
		else {
			commentator().stop("Format SparsePar pass");
		}
	}

	{ /*  SparseMap */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::SparseMap>", "SparseMap");
		SparseMatrix<Field, SparseMatrixFormat::SparseMap> S2(F, m, n);
		SparseMatrix<Field, SparseMatrixFormat::SparseMap> S3(F, m, n);
		buildBySetEntry(S2, N);
		if ( ! testBlackbox(S2,true) ) {
			commentator().stop("Format SparseMap FAIL (testBlackbox)");
			pass = false;
		}
		else if ( ! MD.areEqual(S1,S2) ) {
			commentator().stop("Format SparseMap FAIL (areEqual)");
			pass = false;
		}
		else if (! buildBySetGetEntry(S3,S2) ) {
			commentator().stop("Format SparseMap FAIL (getEntry)");
			pass = false;
		}
		else {
			commentator().stop("Format SparseMap pass");
		}
	}

	{ /*  TPL */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::TPL>", "TPL");
		SparseMatrix<Field, SparseMatrixFormat::TPL> S2(F, m, n);
		SparseMatrix<Field, SparseMatrixFormat::TPL> S3(F, m, n);
		buildBySetEntry(S2, N);
		if ( ! testBlackbox(S2,true) ) {
			commentator().stop("Format TPL FAIL (testBlackbox)");
			pass = false;
		}
		else if ( ! MD.areEqual(S1,S2) ) {
			commentator().stop("Format TPL FAIL (areEqual)");
			pass = false;
		}
		else if (! buildBySetGetEntry(S3,S2) ) {
			commentator().stop("Format TPL FAIL (getEntry)");
			pass = false;
		}
		else {
			commentator().stop("Format TPL pass");
		}
	}

	{ /*  Default OLD */
		commentator().start("SparseMatrix<Field>", "Field");
		Protected::SparseMatrixGeneric<Field> S11(F, m, n);
		buildBySetEntry(S11, N);
		if ( testBlackbox(S1,true)  && MD.areEqual(S1,S11))
			commentator().stop("SparseMatrix<Field> pass");
		else {
			commentator().stop("SparseMatrix<Field> FAIL");
			pass = false;
		}
	}

	commentator().stop( MSG_STATUS(pass), "Sparse matrix black box test suite pass");


	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

