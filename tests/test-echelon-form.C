/* Copyright (C) 2012 bds
 * tests/test-echelon-form.C
 *
 * Adapted by bds from test-blas-domain written by Pascal Giorgi 
 *
 * ---------------------------------------------------------
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
 *.
 *
 *
 */


/*! @file  tests/test-echelon-form.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */


#include "linbox/linbox-config.h"
#include "test-common.h"
#include <iostream>
#include <string>


#include "linbox/ring/modular.h"
//#include "linbox/field/givaro.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/util/commentator.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/solutions/echelon.h"
#include "givaro/givranditer.h"

using namespace LinBox;

const int maxpretty = 35;

string blank;

const char* pretty(string a)
{

	blank = "     " + a;
	int msgsize= maxpretty - (int)blank.size();
	string dot(".");
	for (int i=0;i<msgsize ;++i)
		 blank+=dot;
	 return blank.c_str();
}
#define mycommentator commentator

/*
 *  Testing the rank of dense matrices using BlasDomain
 *  construct a n*n matrices of rank r and compute the rank
 */
template <class Field>
static bool testRank (const Field& F, size_t m, size_t n, int iterations = 1)
{

	typedef typename Field::Element Element;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing rank"),"testRank",(size_t)iterations);

	typename Field::RandIter G(F);
	typename Field::NonZeroRandIter Gn(G);
	Element tmp;
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
        //EchelonFormDomain<Field> EFD(F);

	for (int k=0;k<iterations; ++k) {

	size_t r;
		mycommentator().progress(k);
		DenseMatrix<Field> A(F,m,n),S(F,m,n), L(F,m,m);

		int mn = (m < n) ? (int)m :(int) n;
		r = (size_t)(random() % mn);
		// create S as an upper triangular matrix with r nonzero rows
		for (size_t i=0;i<r;++i){
			S.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<n;++j)
				S.setEntry(i,j,G.random(tmp));
		}
                BMD.write(commentator().report(), S) << " = S" << std::endl;


		// create L as a lower triangular matrix with nonzero elements on the diagonal
		for (size_t i=0;i<m;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,Gn.random(tmp));
		}
                BMD.write(commentator().report(), L) << " = L" << std::endl;

		//  compute A=LS
		BMD.mul(A,L,S);
                BMD.write(commentator().report(), A) << " = A" << std::endl;

		// compute the rank of A
		DenseMatrix<Field> E1(F,m,n), E2(F,m,n), E3(F,m,n), E4(F,m,n);
		size_t rank1 = rowEchelon(E1, A);
        BMD.write(commentator().report(), E1) << " = rowEchelon(E1, A)" << std::endl;

        size_t rank2 = reducedRowEchelon(E2, A);
        BMD.write(commentator().report(), E2) << " = rowReducedEchelon(E2, A)" << std::endl;
		size_t rank3= colEchelon(E3, A);
                BMD.write(commentator().report(), E3) << " = columnEchelon(E3, A)" << std::endl;

		size_t rank4 = reducedColEchelon(E4, A);
                BMD.write(commentator().report(), E4) << " = columnReducedEchelon(E4, A)" << std::endl;

		size_t rank5 = colEchelonize(A);
                BMD.write(commentator().report(), A) << " = columnEchelon(A)" << std::endl;

		size_t rank6 = reducedColEchelonize(E1);
                BMD.write(commentator().report(), E1) << " = columnReducedEchelon(E1)" << std::endl;

		commentator().report() << "Ranks " << rank1 << " " << rank2 << " " << rank3 << " " << rank4 << " " << rank5 << " " << rank6 << " should be " << r << std::endl;

		if (rank1!=r or rank2 !=r  or rank3 !=r  or rank4 !=r  or rank5 !=r  or rank6 !=r)
			ret=false;
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testRank");

	return ret;
}

/*
 * Test of the PLUQMatrix class
 */
template <class Field>
static bool testPLUQ (const Field& F, size_t m, size_t n, int iterations = 1)
{

	typedef typename Field::Element                  Element;
	typedef DenseMatrix<Field>                       Matrix;
	typedef typename Field::RandIter                RandIter;
	typedef typename Field::NonZeroRandIter         NonZeroRandIter;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing PLUQ factorization"),"testPLUQ", iterations);

	RandIter G(F);
	NonZeroRandIter Gn(G);

	bool ret = true;
	MatrixDomain<Field> MD(F);
	BlasMatrixDomain<Field> BMD(F);
	Element tmp; F.init(tmp);

	for (int k=0;k<iterations;++k) {

		mycommentator().progress(k);

		Matrix A(F, m,n), Abis(F, m,n), B(F, m,m), C(F, m,n);


		// Create B a random matrix of rank n/2
		for (size_t j=0;j<m;++j)
			if ( j % 2 )
				for (size_t i=0;i<m;++i)
				  B.setEntry(i,j,G.random(tmp));
			else
			  for (size_t i=0;i<m;++i)
			    B.setEntry(i,j,F.zero);
		// Create C a random matrix of rank n/2
		for (size_t i=0;i<m;++i)
			if ( i % 2 )
				for (size_t j=0;j<n;++j)
					C.setEntry(i,j,G.random(tmp));
			else
				for (size_t j=0;j<n;++j)
					C.setEntry(i,j,F.zero);

		// A = B*C
		BMD.mul(A, B, C);

		Abis = A;

		BlasPermutation<size_t>  P(A.coldim()),Q(A.rowdim());
		PLUQMatrix<Field> X(A,P,Q);

        Matrix Ldata(F,m,m), Udata(F,m,n);
		TriangularBlasMatrix<Matrix> L(Ldata,Tag::Shape::Lower,Tag::Diag::Unit);
		TriangularBlasMatrix<Matrix> U(Udata,Tag::Shape::Upper,Tag::Diag::NonUnit);
		X.getL(L);
		X.getU(U);
		P=X.getP();

		Q=X.getQ();

		// C = U*Q
		BMD.mul( C, U, Q);
		// A = L*C
		BMD.mul( A, L, C);
		// C = Q*C
		BMD.mulin_right( P, C);

		if (!MD.areEqual(A,Abis))
            ret=false;
        
        // Second pass
		// A = B*C
		BMD.mul(A, B, C);

		Abis = A;

		PLUQMatrix<Field> Y(A,P,Q);

        Matrix L2data(F,m,m), U2data(F,m,n);
		TriangularBlasMatrix<Matrix> L2(L2data,Tag::Shape::Lower,Tag::Diag::Unit);
		TriangularBlasMatrix<Matrix> U2(U2data,Tag::Shape::Upper,Tag::Diag::NonUnit);

	
		Y.getL(L2);
		Y.getU(U2);
		P=Y.getP();

		Q=Y.getQ();

		// C = Q*U2
		BMD.mul( C,L2,U2);
		// C = C*Q
		BMD.mulin_left(  C,Q);
		// A = L*C
		BMD.mul( A, P, C);

		if (!MD.areEqual(A,Abis))        
            ret=false;        
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testPLUQ");

	return ret;
}

// returns true if ok, false if not.
template<class Field>
int launch_tests(Field & F, size_t m, size_t n, int iterations = 1)
{
	bool pass = true ;
	size_t mn = (m < n) ? n : m;
 	if (!testRank (F, mn, mn, iterations))     pass=false;
	if (m != n) {
 		if (!testRank (F, n, m, iterations))   pass=false;
 		if (!testRank (F, m, n, iterations))   pass=false;
	}
    if (!testPLUQ (F,n,n,iterations))                     pass=false;
	return pass ;

}

int main(int argc, char **argv)
{

	static size_t n = 20;
	static size_t m = 10;
	static integer q = 1000003U;
	static int iterations = 1;

    static Argument args[] = {
        { 'm', "-n M", "Set dimension of test matrices to MxN", TYPE_INT,     &m },
        { 'n', "-n N", "Set dimension of test matrices to MxN", TYPE_INT,     &n },
        { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1]",  TYPE_INTEGER, &q },
        { 'i', "-i I", "Perform each test for I iterations",    TYPE_INT,     &iterations },
	END_OF_ARGUMENTS
    };

	parseArguments (argc, argv, args);

	typedef Givaro::Modular<double> Field;
	//typedef Givaro::Modular<int> Field;
	//typedef Givaro::Modular<float> Field;

	Field F1 (q);
	bool pass = true;

	srand ((unsigned)time (NULL));


	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	commentator().start("BlasMatrixDomain test suite", "BlasMatrixDomain");

	pass &= launch_tests(F1,m, n, iterations);
	commentator().stop(MSG_STATUS (pass), (const char *) 0,"BlasMatrixDomain test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
