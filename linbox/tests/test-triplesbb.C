/* Copyright (C) LinBox
 *
 * using generic testBlackbox  -bds
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file  tests/test-triplesbb.C
 * @ingroup tests
 *
 * @brief no doc
 *
 * @test no doc.
 */


#include "linbox-config.h"

#include <iostream>
#include <fstream>


#include "linbox/field/modular.h"
// #include "linbox/blackbox/triplesbb.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/blas-domain.h"
#include "test-blackbox.h"
#include "test-common.h"

using namespace LinBox;

template<class BB>
void randBuild(BB & A, size_t nnz){
	for(size_t i = 0; i < nnz; ++i)
	{	typename BB::Field::Element d; A.field().init(d, rand());
		size_t ii = rand()%A.rowdim();
		size_t jj = rand()%A.coldim();
		A.setEntry(ii,jj, d);
	}
}

template<class BB>
void bidiag(BB & A){
	typename BB::Field::Element d; A.field().init(d, 1);
	A.setEntry(0,0, d);
	size_t n = A.coldim() > A.rowdim() ? A.rowdim() : A.coldim();
	for(size_t i = 1; i < n; ++i)
	{	A.setEntry(i,i, d);
		A.setEntry(i,i-1, d);
	}
}

int main (int argc, char **argv)
{
	// ofstream report;

	bool pass = true;

	static size_t m = 4;
	static size_t n = 20;
	static size_t nnz = 0;
	static integer q = 2147483647U;
	q = 101;

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of test matrix to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set col dimension of test matrix to N.", TYPE_INT,     &n },
		{ 'z', "-n NNZ", "Set number of nonzero entries in test matrix.", TYPE_INT,     &nnz },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	if (nnz == 0) nnz = m*n/10;

	srand ((unsigned)time (NULL));

	commentator().start("TriplesBB black box test suite", "triplesbb");

	//Field
	typedef Modular<double> Field;
	typedef Field::Element Element;
	Field F (q);
	MatrixDomain<Field> MD(F);
	typedef MatrixDomain<Field>::OwnMatrix OwnMatrix;
	//Vectors
	VectorDomain<Field> VD(F);
	std::vector<Element> x(n), y(m), z(m);
	for (size_t i = 0; i < n; ++i) F.init(x[i], i+1);

	// TriplesBB<Field> A(F, m, n);
	SparseMatrix<Field,SparseMatrixFormat::TPL> A(F, m, n);
	randBuild(A, nnz);
	pass = pass && testBlackbox(A);
	// just see if it compiles
	size_t b = 10;
	OwnMatrix X(F,n,b), Y(F,m,b), Z(F,b,m), W(F,b,n);
	X.random(); Z.random();
	A.applyLeft(Y, X);
	A.applyRight(W, Z);

	// standard constructor
	// TriplesBB<Field> B(F, m, n);
	SparseMatrix<Field,SparseMatrixFormat::TPL> B(F, m, n);
	bidiag(B);
	pass = pass && testBlackbox(B);
	B.apply(y, x);
	// default cstor plus init
	// TriplesBB<Field> C(F);
	SparseMatrix<Field,SparseMatrixFormat::TPL> C(F);
      	C.init(F, m, n);
	bidiag(C);
	pass = pass && testBlackbox(C);
	// check B == C
	C.apply(z, x);
	if (not VD.areEqual(y, z)) {
		pass = false;
		LinBox::commentator().report() << "fail: cstor and init disagree" << std::endl;
	}
	// copy construction
	// TriplesBB<Field> D(B);
	SparseMatrix<Field,SparseMatrixFormat::TPL> D(B);
	pass = pass && testBlackbox(D);
	// check B == D
	D.apply(z, x);
	if (not VD.areEqual(y, z)) {
		pass = false;
		LinBox::commentator().report() << "copy cstor failure" << std::endl;
	}

	// check that it's deep copy
	Element a;
	B.getEntry(a,0,0);
	F.addin(a, F.one);
	D.setEntry(0,0,a);
	D.apply(z, x);
	if (VD.areEqual(y, z)) {
		pass = false;
		LinBox::commentator().report() << "fail changed copy cstor and original agree" << std::endl;
		VD.write(commentator().report() << "y ", y) << std::endl;
		VD.write(commentator().report() << "z ", z) << std::endl;
	}

	commentator().stop("TriplesBB black box test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

