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
		// r.random(x); // I want to see what happens when one reads in zero. If we don't want it, we just stop permitting setting zero... (hence the clearEntry function)
		// std::cout << "A[ " << i+1 << ',' << j+1 << "]:=" << x << ';' << std::endl;
		// if (A.field().isZero(x)) std::cout << "#is zero" << std::endl;
		A.setEntry(i,j,x);
	}
	A.finalize();
	// std::cout << "B :=" << std::endl;
	// A.write(std::cout);
	// std::cout << ';' << std::endl;
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
	MatrixDomain<Field> MD ;

	 /*  default */
	commentator().start("SparseMatrix<Field>", "Field");
	SparseMatrix<Field> S1(F, m, n);
	buildBySetEntry(S1, N);
	if ( testBlackbox(S1,false))
		commentator().stop("SparseMatrix<Field> pass");
	else {
		commentator().stop("SparseMatrix<Field> FAIL");
		pass = false;
	}

	{ /*  COO */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::COO>", "COO");
		SparseMatrix<Field, SparseMatrixFormat::COO> S2(F, m, n);
		buildBySetEntry(S2, N);
		if ( testBlackbox(S2,false)  && MD.areEqual(S1,S2) )
			commentator().stop("Format COO pass");
		else {
			commentator().stop("Format COO FAIL");
			pass = false;
		}
	}

	{ /*  CSR  */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::CSR>", "CSR");
		SparseMatrix<Field, SparseMatrixFormat::CSR> S3(F, m, n);
		buildBySetEntry(S3, N);
		if ( testBlackbox(S3,false)  && MD.areEqual(S1,S3))
			commentator().stop("Format CSR pass");
		else {
			commentator().stop("Format CSR FAIL");
			pass = false;
		}
	}

	{ /*  ELL  */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::ELL>", "ELL");
		commentator().report() << "SparseMatrix<Field, SparseMatrixFormat::ELL>" << std::endl;
		SparseMatrix<Field, SparseMatrixFormat::ELL> S4(F, m, n);
		buildBySetEntry(S4, N);
		if ( testBlackbox(S4,false)  && MD.areEqual(S1,S4))
			commentator().stop("Format ELL pass");
		else {
			commentator().stop("Format ELL FAIL");
			pass = false;
		}
	}

	{ /*  ELL_R  */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::ELL_R>", "ELL_R");
		SparseMatrix<Field, SparseMatrixFormat::ELL_R> S5(F, m, n);
		buildBySetEntry(S5, N);
		if ( testBlackbox(S5,false)  && MD.areEqual(S1,S5))
			commentator().stop("Format ELL_R pass");
		else {
			commentator().stop("Format ELL_R FAIL");
			pass = false;
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

	{ /* Vector of row, row is Vector of index/value Pair */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::SparseSeq>", "VVP");
		SparseMatrix<Field, SparseMatrixFormat::SparseSeq> S7(F, m, n);
		buildBySetEntry(S7, N);
		if ( testBlackbox(S7,true)  && MD.areEqual(S1,S7))
			commentator().stop("Format VVP pass");
		else {
			commentator().stop("Format VVP FAIL");
			pass = false;
		}
	}

	{ /*  SparsePar */
		// Vector of row, row is Pair of Vectors, one of indices, one of values.
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::SparsePar>", "VPV");
		SparseMatrix<Field, SparseMatrixFormat::SparsePar> S8(F, m, n);
		buildBySetEntry(S8, N);
		if ( testBlackbox(S8,true)  && MD.areEqual(S1,S8))
			commentator().stop("Format VPV pass");
		else {
			commentator().stop("Format VPV FAIL");
			pass = false;
		}
	}

	{ /*  SparseMap */
		// Vector of row, row is Map, index to value.
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::SparseMap>", "VMap");
		SparseMatrix<Field, SparseMatrixFormat::SparseMap> S9(F, m, n);
		buildBySetEntry(S9, N);
		if ( testBlackbox(S9,true)  && MD.areEqual(S1,S9))
			commentator().stop("Format VMap pass");
		else {
			commentator().stop("Format VMap FAIL");
			pass = false;
		}
	}


	{ /* Vector of i,j,val triples */
		commentator().start("SparseMatrix<Field, SparseMatrixFormat::TPL>", "TPL");
		SparseMatrix<Field, SparseMatrixFormat::TPL> S10(F, m, n);
		buildBySetEntry(S10, N);
		if ( testBlackbox(S10,true)  && MD.areEqual(S1,S10))
			commentator().stop("Format TPL pass");
		else {
			commentator().stop("Format TPL FAIL");
			pass = false;
		}

		if (pass)
			commentator().stop("Sparse matrix black box test suite pass");
		else
			commentator().stop("Sparse matrix black box test suite FAIL");
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

	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

