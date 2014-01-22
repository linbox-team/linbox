
/* tests/test-sparse.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
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

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>


#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/matrix/sparse.h"
#include "linbox/matrix/sparse-matrix.h"

#include "test-blackbox.h"

using namespace LinBox;

template <class SM>
void buildBySetEntry(SM & A, size_t nnz)
{
	typename SM::Field::RandIter r(A.field());
	//size_t i, j, k;
	typename SM::Field::Element x;

	for (size_t k = 0; k < nnz; ++k)
	{
		size_t i = random() % A.rowdim();
		size_t j = random() % A.coldim();
		r.nonzerorandom(x);
		A.setEntry(i,j,x);
	}
	A.finalize();
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static size_t m = 10;
	static integer q = 101;
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

	commentator().start("SparseMatrix<Field>", "Field");
	SparseMatrix<Field> S1(F, m, n);
	buildBySetEntry(S1, N);
	if ( testBlackboxNoRW(S1) )
		commentator().stop("SparseMatrix<Field> pass");
	else {
		commentator().stop("SparseMatrix<Field> FAIL");
		pass = false;
	}

#if 0 // using CSR, segfaults
	commentator().start("SparseMatrix2<Field>", "Field");
	SparseMatrix2<Field> S1a(F, m, n);
	buildBySetEntry(S1a, N);
	if ( testBlackboxNoRW(S1) )
		commentator().stop("SparseMatrix2<Field> pass");
	else {
		commentator().stop("SparseMatrix2<Field> FAIL");
		pass = false;
	}
#endif

#if 0 // acts like zero
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::COO>", "COO");
	SparseMatrix2<Field, SparseMatrixFormat::COO> S2(F, m, n);
	buildBySetEntry(S2, N);
	if ( testBlackboxNoRW(S2) )
		commentator().stop("Format COO pass");
	else {
		commentator().stop("Format COO FAIL");
		pass = false;
	}
#endif

#if 0 // segfaults
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::CSR>", "CSR");
	SparseMatrix2<Field, SparseMatrixFormat::CSR> S3(F, m, n);
	buildBySetEntry(S3, N);
	if ( testBlackboxNoRW(S3) )
		commentator().stop("Format CSR pass");
	else {
		commentator().stop("Format CSR FAIL");
		pass = false;
	}
#endif

#if 0 // acts like zero
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::ELL>", "ELL");
	commentator().report() << "SparseMatrix2<Field, SparseMatrixFormat::ELL>" << std::endl;
	SparseMatrix2<Field, SparseMatrixFormat::ELL> S4(F, m, n);
	buildBySetEntry(S4, N);
	if ( testBlackboxNoRW(S4) )
		commentator().stop("Format ELL pass");
	else {
		commentator().stop("Format ELL FAIL");
		pass = false;
	}
#endif

#if 0 // looks like zero
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::ELL_R>", "ELL_R");
	SparseMatrix2<Field, SparseMatrixFormat::ELL_R> S5(F, m, n);
	buildBySetEntry(S5, N);
	if ( testBlackboxNoRW(S5) )
		commentator().stop("Format ELL_R pass");
	else {
		commentator().stop("Format ELL_R FAIL");
		pass = false;
	}
#endif

#if 0 // doesn't compile
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::HYB>", "HYB");
	SparseMatrix2<Field, SparseMatrixFormat::HYB> S6(F, m, n);
	buildBySetEntry(S6, N);
	S6.optimise();
	if ( testBlackboxNoRW(S6) )
		commentator().stop("Format HYB pass");
	else {
		commentator().stop("Format HYB FAIL");
		pass = false;
	}
#endif

#if 0 /* not ported yet */
// Vector of row, row is Vector of index/value Pair
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::VVP>", "VVP");
	SparseMatrix2<Field, SparseMatrixFormat::VVP> Svvp(F, m, n);
	buildBySetEntry(Svvp, N);
	if ( testBlackboxNoRW(Svvp) )
		commentator().stop("Format VVP pass");
	else {
		commentator().stop("Format VVP FAIL");
		pass = false;
	}
#endif

#if 0 /* not ported yet */
// Vector of row, row is Pair of Vectors, one of indices, one of values.
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::VPV>", "VPV");
	SparseMatrix2<Field, SparseMatrixFormat::VPV> Svpv(F, m, n);
	buildBySetEntry(Svpv, N);
	if ( testBlackboxNoRW(Svpv) )
		commentator().stop("Format VPV pass");
	else {
		commentator().stop("Format VPV FAIL");
		pass = false;
	}
#endif

#if 0 /* not ported yet */
// Vector of row, row is Map, index to value.
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::VMap>", "VMap");
	SparseMatrix2<Field, SparseMatrixFormat::VMap> Svm(F, m, n);
	buildBySetEntry(Svm, N);
	if ( testBlackboxNoRW(Svm) )
		commentator().stop("Format VMap pass");
	else {
		commentator().stop("Format VMap FAIL");
		pass = false;
	}
#endif

#if 1 /* not ported yet */
// Vector of i,j,val triples
	commentator().start("SparseMatrix2<Field, SparseMatrixFormat::TPL>", "TPL");
	SparseMatrix2<Field, SparseMatrixFormat::TPL> Striple(F, m, n);
	buildBySetEntry(Striple, N);
	if ( testBlackboxNoRW(Striple) )
		commentator().stop("Format TPL pass");
	else {
		commentator().stop("Format TPL FAIL");
		pass = false;
	}

	if (pass)
		commentator().stop("Sparse matrix black box test suite pass");
	else
		commentator().stop("Sparse matrix black box test suite FAIL");
	return pass ? 0 : -1;
}
#endif

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

