/* tests/test-sparse.C
 * Copyright (C) 2014 the LinBox group
 *
 * Written by bds <saunders@udel.edu>
 *            Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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
#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"


#include "test-blackbox.h"

using namespace LinBox;

template <class SM, class SM2>
bool buildBySetGetEntry(SM & A, const SM2 &B);

template <class Field, class SMF>
bool testSparseFormat(string format, const SparseMatrix<Field> & S1)
{
	typedef SparseMatrix<Field, SMF> SM;
	bool pass=true;
	string msg = "SparseMatrix<Field, SparseMatrixFormat::" + format + ">";
	commentator().start(msg.c_str(), format.c_str());
	const typename SM::Field& F = S1.field();
	MatrixDomain<typename SM::Field> MD(F);
	SM S2(F,S1.rowdim(),S1.coldim());
	SM S3(F,S1.rowdim(),S1.coldim());
	buildBySetGetEntry(S2, S1);
	buildBySetGetEntry(S3, S2);
	// behaves as a blackbox should
	if ( ! testBlackbox(S2,true) ) {
		msg = format + " FAIL(testBlackbox)";
		commentator().stop(msg.c_str());
		pass = false;
	}
	// equal behaviour when representing same mat
	else if ( ! MD.areEqual(S1,S2) ) { 
		msg = format + " FAIL(areEqual)";
		commentator().stop(msg.c_str());
		pass = false;
	}
	// just a check of S2 getEntry
	else if ( ! MD.areEqual(S3,S2) ) { 
		msg = format + " FAIL(getEntry)";
		commentator().stop(msg.c_str());
		pass = false;
	}
	else {
		msg = format + " pass";
		commentator().stop(msg.c_str());
	}
	return pass;
}

template <class SM, class SM2>
bool buildBySetGetEntry(SM & A, const SM2 &B)
{

	A.resize(B.rowdim(),B.coldim());

	for (size_t i = 0; i < B.rowdim(); ++i)
		for (size_t j = 0; j < B.coldim(); ++j)
		{
			typename SM::Field::Element x = B.getEntry(i,j);
			if (not A.field().isZero(x))
			A.setEntry(i,j,x);
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
		{ 'm', "-m m", "Set row dimension of test matrices to m.", TYPE_INT,     &m },
		{ 'n', "-n n", "Set column dimension of test matrices to n.", TYPE_INT,     &n },
		{ 'N', "-N N", "Set number of nonzeros in test matrices.", TYPE_INT,     &N },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	//typedef	Givaro::Modular<uint32_t> Field;
	typedef	Givaro::Modular<double> Field;

	Field F (q);

	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator().start("Sparse matrix black box test suite", "Sparse");
	MatrixDomain<Field> MD(F) ;
	typename Field::RandIter r(F,1);
	srand(0);

	 /*  default case */
	commentator().start("SparseMatrix<Field>", "Field");

	SparseMatrix<Field> S1(F, m, n);
	typename Field::Element x;
	for (size_t k = 0; k < N; ++k)
	{
		size_t i = rand() % m;
		size_t j = rand() % n;
		while (S1.field().isZero(r.random(x)));
		S1.setEntry(i,j,x);
	}
	S1.finalize();

	//if ( testBlackbox(S1, true, N != 0))
	if ( testBlackbox(S1, false, N != 0))
		commentator().stop("SparseMatrix<Field> pass");
	else {
		commentator().stop("SparseMatrix<Field> FAIL");
		pass = false;
	}

	/* other formats */
	pass = pass and 
		testSparseFormat<Field, SparseMatrixFormat::COO>("COO",S1);
	pass = pass and 
		testSparseFormat<Field, SparseMatrixFormat::CSR>("CSR",S1);
	pass = pass and 
		testSparseFormat<Field, SparseMatrixFormat::ELL>("ELL",S1);
	pass = pass and 
		testSparseFormat<Field, SparseMatrixFormat::ELL_R>("ELL_R",S1);
	pass = pass and 
		testSparseFormat<Field, SparseMatrixFormat::TPL>("TPL",S1);
	pass = pass and 
		testSparseFormat<Field, SparseMatrixFormat::SparseSeq>("SparseSeq",S1);
	pass = pass and 
		testSparseFormat<Field, SparseMatrixFormat::SparsePar>("SparsePar",S1);
	pass = pass and 
		testSparseFormat<Field, SparseMatrixFormat::SparseMap>("SparseMap",S1);
#if 0 // doesn't compile
	commentator().start("SparseMatrix<Field, SparseMatrixFormat::HYB>", "HYB");
	SparseMatrix<Field, SparseMatrixFormat::HYB> S6(F, m, n);
	buildBySetGetEntry(S6, S1);
	S6.optimise();
	if ( testBlackbox(S6,false)  && MD.areEqual(S1,S6))
		commentator().stop("Format HYB pass");
	else {
		commentator().stop("Format HYB FAIL");
		pass = false;
	}
#endif

	{ /*  Default OLD */
		commentator().start("SparseMatrix<Field>", "Field");
		Protected::SparseMatrixGeneric<Field> S11(F, m, n);
		buildBySetGetEntry(S11, S1);
		if ( testBlackbox(S11,true)  && MD.areEqual(S1,S11))
			commentator().stop("SparseMatrix<Field> pass");
		else {
			commentator().stop("SparseMatrix<Field> FAIL");
			pass = false;
		}
	}

	commentator().stop( MSG_STATUS(pass), "Sparse matrix black box test suite pass");


	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
