
/* tests/test-blas-domain.C
 * Copyright (C) 2004 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
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


/*! @file  tests/test-blas-domain.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */


 // where is this used?
#define __MINP_CONSTRUCT
#include "linbox/linbox-config.h"
#include <iostream>
#include <string>
#include "linbox/integer.h"
#include "linbox/field/gf2.h"
#include "linbox/field/modular.h"
#include "linbox/field/modular-balanced.h"
#include "linbox/field/givaro.h"
#ifdef __LINBOX_HAVE_NTL
#include "linbox/field/ntl.h"
#endif
#include "linbox/matrix/blas-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/randiter/nonzero.h"
#include "linbox/util/commentator.h"
#include "linbox/algorithms/blas-domain.h"
#include "linbox/field/PID-integer.h"
// #include "linbox/algorithms/matrix-hom.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/blackbox/scalar-matrix.h"



#include "test-common.h"

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


template <class Field>
static bool testMulAdd (const Field& F, size_t n, int iterations)
{

	typedef typename Field::Element     Element;
	typedef typename Field::RandIter   RandIter;
	typedef BlasMatrix<Field>          Matrix;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing muladd"),"testMulAdd",(unsigned int)iterations);

	RandIter G(F);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);
	MatrixDomain<Field>      MD(F);
	VectorDomain<Field>      VD(F);

	Matrix A;
	for (int k=0;k<iterations; ++k) {

		mycommentator().progress(k);
		A.init(F, n, n);
		Matrix /*A(F, n,n),*/B(F, n,n),C(F, n,n),D(F, n,n),T(F, n,n),R(F, n,n);
		std::vector<Element> x(n),y(n),z(n),t(n);

		Element alpha, beta,malpha,tmp;


		// Create 3 random n*n matrices
		A.random();
		B.random();
		C.random();

		// Create 2 random vectors
		for (size_t i=0;i<n;++i) {
			G.random(x[i]);
			G.random(y[i]);
		}

		// create 2 random element
		G.random(alpha);
		G.random(beta);

		F.neg(malpha,alpha);

		// compute D = -alpha.(A*C+B*C) + alpha.(A+B)*C

		BMD.mul(D,A,C);
		BMD.mul(T,B,C);
		MD.addin(D,T);

		MD.add(T,A,B);
		BMD.muladd(R,malpha,D,alpha,T,C);

		if (!MD.isZero(R))
			ret=false;

		// compute z = beta.y + alpha.A*x

		BMD.muladd(z,beta,y,alpha,A,x);

		MD.vectorMul(t,A,x);
		for (size_t i=0;i<n;++i){
		  F.mulin(t[i],alpha);
		  F.axpyin(t[i],beta,y[i]);
		}

		if (!VD.areEqual(t,z))
			ret=false;
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testMulAdd");

	return ret;
}

// computes D = alpha A B + beta C on integers and check the result is ok mod p.
// actually we check the mod p muladd here...
bool CheckMulAdd( const Integer & alpha ,
		  const BlasMatrix<PID_integer> & A ,
		  const BlasMatrix<PID_integer> & B ,
		  const Integer & beta ,
		  const BlasMatrix<PID_integer> & C)
{

	size_t M = C.rowdim();
	size_t N = C.coldim();

	typedef Modular<double>       Field ;
	typedef Field::Element      Element ;

	PID_integer ZZ ;
	MatrixDomain<PID_integer> ZMD(ZZ);

	BlasMatrix<PID_integer> D(ZZ,M,N);

	Integer p = Integer::random_between(10,12) ;
	nextprime(p,p); //!@bug si p n'est pas premier, fgemm fait n'importe quoi (division par alpha)
	Field Zp (p);

	BlasMatrixDomain<Field> BMD (Zp);
	MatrixDomain<Field>      MD (Zp);

	// Ep = b C + a A B
	ZMD.muladd(D,beta,C,alpha,A,B);
	BlasMatrix<Field> Dp(D,Zp); // D mod p

	BlasMatrix<Field> Ap(A,Zp);
	BlasMatrix<Field> Bp(B,Zp);
	BlasMatrix<Field> Cp(C,Zp);
	// BlasMatrix<Field> Ap(A.rowdim(),A.coldim());
	// BlasMatrix<Field> Bp(B.rowdim(),B.coldim());
	// BlasMatrix<Field> Cp(C.rowdim(),C.coldim());
	// MatrixHom::map(Ap,A,Zp);
	// MatrixHom::map(Bp,B,Zp);
	// MatrixHom::map(Cp,C,Zp);
	BlasMatrix<Field> Ep(Zp,M,N);  // D mod p

	Element ap, bp ;
	Zp.init(ap,alpha);
	Zp.init(bp,beta);

	// Ep = bp Cp + ap Ap Bp mod p
	BMD.muladd(Ep,bp,Cp,ap,Ap,Bp);

	bool pass = MD.areEqual(Ep,Dp);
	if (!pass) {
#if 0 /*  maple check on stdout */
		std::cout << "#########################################" << std::endl;
		std::cout << "p := " << p << ';' << std::endl;
		std::cout << "ap,bp := " << ap << ',' << bp << ';' << std::endl;
		Ap.write(std::cout << "Ap :=", LinBoxTag::FormatMaple) << ";" << std::endl;
		Bp.write(std::cout << "Bp :=", LinBoxTag::FormatMaple) << ";" << std::endl;
		Cp.write(std::cout << "Cp :=", LinBoxTag::FormatMaple) << ";" << std::endl;
		Dp.write(std::cout << "Dp :=", LinBoxTag::FormatMaple) << ";" << std::endl;
		Ep.write(std::cout << "Ep :=", LinBoxTag::FormatMaple) << ";" << std::endl;
		std::cout << "alpha,beta := " << alpha << ',' << beta << ';' << std::endl;
		A.write(std::cout << "A :=",LinBoxTag::FormatMaple) << ';' << std::endl;
		B.write(std::cout << "B :=",LinBoxTag::FormatMaple) << ';' << std::endl;
		C.write(std::cout << "C :=",LinBoxTag::FormatMaple) << ';' << std::endl;
		D.write(std::cout << "E :=",LinBoxTag::FormatMaple) << ';' << std::endl;
		std::cout << "evalm(E-alpha*A.B-beta*C);" << std::endl;
		std::cout << "#########################################" << std::endl;
#endif
		mycommentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") *** " << std::endl;
	}

	// Ep = bp Cp + ap Ap Bp mod p
	MD.muladd(Ep,bp,Cp,ap,Ap,Bp);
	bool all = MD.areEqual(Ep,Dp);
	if (!all) {
		mycommentator().report() << " *** MD ERROR *** " << std::endl;
	}

	return pass&all ;

}

// tests MulAdd for various parameters alpha and beta.
template <class Field>
static bool testMulAddAgain (const Field& , size_t n, int iterations)
{

	PID_integer ZZ ;
	typedef BlasMatrix<PID_integer>         IMatrix;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing muladd again"),"testMulAddAgain",(unsigned int)iterations);

	bool ret = true;

	size_t ll = 17 ;
	Integer::seeding();
#if 0
	size_t lA = 15 ;
	size_t lB = 18 ;
	size_t lC = 19 ;

	PID_integer ZZ ;

	typedef RandomIntegerIter<>                       IntRandIter ;
	typedef RandomDenseMatrix<IntRandIter, PID_integer > IntRand_t;

	IntRandIter RA(lA);
	IntRandIter RB(lB);
	IntRandIter RC(lC);

	IntRand_t Arand (ZZ,RA);
	IntRand_t Brand (ZZ,RB);
	IntRand_t Crand (ZZ,RC);
#endif

	for (int k=0;k<iterations; ++k) {

		mycommentator().progress(k);
		IMatrix A(ZZ, n,n),B(ZZ, n,n),C(ZZ, n,n),D(ZZ, n,n),E(ZZ, n,n);

		Integer a , b;
#if 0
		Arand.random<IMatrix>(A);
		Brand.random<IMatrix>(B);
		Crand.random<IMatrix>(C);
#endif
		for (size_t i=0;i<n;++i)
			for (size_t j=0;j<n;++j){
				A.setEntry(i,j,Integer::random());
				B.setEntry(i,j,Integer::random());
				C.setEntry(i,j,Integer::random());
			}

		a = 1 ; b = 1 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = 1 ; b = -1 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = -1 ; b = 1 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = -1 ; b = -1 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = 0 ; b = 1 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = 1 ; b = 0 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = 0 ; b = -1 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = -1 ; b = 0 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = Integer::random<false>(ll) ; b = 1 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a =Integer::random<false>(ll) ; b = -1 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a =  1 ; b = Integer::random<false>(ll) ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = -1 ; b = Integer::random<false>(ll) ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = 0 ; b = Integer::random<false>(ll) ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = Integer::random<false>(ll) ; b = 0 ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;
		a = Integer::random<false>(ll) ; b = Integer::random<false>(ll) ;
		if (!CheckMulAdd(a,A,B,b,C)) ret = false ;

	}


	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testMulAddAgain");

	return ret;
}

// tests MulAdd for various shapes and values of transposition.
template <class Field>
static bool testMulAddShapeTrans (const Field &F, size_t m, size_t n, size_t k, int iterations)
{
	bool ret = true ;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing muladd for shapes and transposition"),"testMulAddShapeTrans",(unsigned long)iterations);


	typedef typename Field::Element Element;
	typedef BlasMatrix<Field> Matrix ;
	typedef TransposedBlasMatrix<Matrix> TransposedMatrix ;
	typedef typename Field::RandIter Randiter ;
	Randiter R(F) ;
	RandomDenseMatrix<Randiter,Field> RandMat(F,R);

	BlasMatrixDomain<Field> BMD (F);
	MatrixDomain<Field>      MD (F);

	// input matrix
	Matrix A(F, m,k);
	Matrix B(F, k,n);
	Matrix C(F, m,n);
	// result matrix
	Matrix D(F, m,n);
	Matrix E(F, m,n);

	// random A,B
	RandMat.random(A);
	RandMat.random(B);
	RandMat.random(C);

	// hard tranpose A,B
	Matrix A1 (F, k,m) ;
	A.transpose(A1) ;
	Matrix B1 (F, n,k) ;
	B.transpose(B1) ;
	TransposedMatrix tA(A1); // t(tA)=A
	TransposedMatrix tB(B1); // t(tB)=B

	// random alpha, beta
	Element alpha ;
	do {R.random(alpha);} while (F.isZero(alpha)); // nonzerorandom

	Element beta ;
	R.random(beta);

	// témoin.
	MD.muladd(D,beta,C,alpha,A,B);

	// A,B
	BMD.muladd(E,beta,C,alpha,A,B);
	if (!MD.areEqual(E,D)) {
		ret = false ;
		mycommentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") (noTrans, noTrans) *** " << std::endl;
	}

	BMD.muladd(E,beta,C,alpha,A,tB);
	if (!MD.areEqual(E,D))  {
		ret = false ;
		mycommentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") (noTrans, Trans) *** " << std::endl;
	}

	BMD.muladd(E,beta,C,alpha,tA,B);
	if (!MD.areEqual(E,D)) {
		ret = false ;
		mycommentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") (Trans, noTrans) *** " << std::endl;
	}

	BMD.muladd(E,beta,C,alpha,tA,tB);
	if (!MD.areEqual(E,D)) {
		ret = false ;
		mycommentator().report() << " *** BMD ERROR (" << alpha << ',' << beta << ") (Trans, Trans) *** " << std::endl;
	}


	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testMulAddShapeTrans");
	return ret ;
}

// tests MulAdd for various shapes and values of transposition.
template<class Field, bool LeftSide, bool UnitDiag>
static bool testTriangMulShapeTrans (const Field &F, size_t m, size_t n, int iterations)
{
	bool ret = true ;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing triangular matmul for shapes and transposition"),"testTriangMulShapeTrans",(unsigned int)iterations);


	typedef typename Field::Element                               Element;
	typedef BlasMatrix<Field>                                   Matrix ;
	typedef TriangularBlasMatrix<Field>               TriangularMatrix ;
	typedef TransposedBlasMatrix<Matrix>                TransposedMatrix ;
	typedef TransposedBlasMatrix<TriangularMatrix > TransposedTriangular ;
	typedef typename Field::RandIter                            Randiter ;
	Randiter R(F) ;
	RandomDenseMatrix<Randiter,Field> RandMat(F,R);

	Element one ;
	F.init(one,1);

	BlasMatrixDomain<Field> BMD (F);
	MatrixDomain<Field>      MD (F);

	int k =(int) (LeftSide?m:n) ;
	// input matrix
	Matrix A(F, k,k); // A = L+U-I. Either L or U is unit.
	Matrix B(F, m,n);
	// result matrix
	Matrix D(F, m,n);
	Matrix E(F, m,n);

	// random A,B
	RandMat.random(A);
	RandMat.random(B);

	// hard tranpose A,B
	Matrix A1 (F, k,k) ;
	A.transpose(A1) ;



	/*  test (L+U-I) B+B = LB+UB */
	if (LeftSide)
		BMD.muladd(D,one,B,one,A,B);
	else
		BMD.muladd(D,one,B,one,B,A);

	/****  DIRECT ****/
	{
		/*  L */
		TriangularMatrix L (A, LinBoxTag::Lower,
				    (UnitDiag?LinBoxTag::Unit:LinBoxTag::NonUnit));

		/*  U */
		TriangularMatrix U (A, LinBoxTag::Upper,
				    (!UnitDiag?LinBoxTag::Unit:LinBoxTag::NonUnit));

		/*  make product */
		E = B ;
		// Matrix G(m,n);
		// G = E ;
		Matrix G((const Matrix&)E); //!@warning on n'oublie pas l'esperluette !!!
		if(LeftSide) {
			BMD.mulin_right(L,E) ; // B <- AB
			BMD.mulin_right(U,G) ;
		}
		else {
			BMD.mulin_left(G,L) ;  // B <- BA
			BMD.mulin_left(E,U) ;
		}
		BMD.addin(E,G);

		/*  check equality */
		if (!MD.areEqual(E,D)) {
			ret = false ;
			mycommentator().report() << " *** BMD ERROR (" << (LeftSide?"left":"right") << ',' << (UnitDiag?" L":" U") << " is unit) *** " << std::endl;
		}
		else {
			mycommentator().report() << " direct triangular multiplication ok." << std::endl;
		}
	}
	/****  Transpose ****/
	{
		/*  L */
		TriangularMatrix L1 (A1, LinBoxTag::Lower,
				    (UnitDiag?LinBoxTag::Unit:LinBoxTag::NonUnit));

		/*  U */
		TriangularMatrix U1 (A1, LinBoxTag::Upper,
				    (!UnitDiag?LinBoxTag::Unit:LinBoxTag::NonUnit));

		TransposedTriangular L(L1);
		TransposedTriangular U(U1);
		/*  make product */
		E = B ;
		// Matrix G(m,n);
		// G = E ;
		Matrix G((const Matrix&)E); //!@warning on n'oublie pas l'esperluette !!!
		if(LeftSide) {
			BMD.mulin_right(L,E) ; // B <- AB
			BMD.mulin_right(U,G) ;
		}
		else {
			BMD.mulin_left(G,L) ;  // B <- BA
			BMD.mulin_left(E,U) ;
		}
		BMD.addin(E,G);

		/*  check equality */
		if (!MD.areEqual(E,D)) {
			ret = false ;
			mycommentator().report() << " *** BMD ERROR Transpose (" << (LeftSide?"left":"right") << ',' << (UnitDiag?" L":" U") << " is unit) *** " << std::endl;
		}
		else {
			mycommentator().report() << " transposed triangular multiplication ok." << std::endl;
		}
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testMulAddShapeTrans");
	return ret ;
}


/*
 *  Testing the rank of dense matrices using BlasDomain
 *  construct a n*n matrices of rank r and compute the rank
 */
template <class Field>
static bool testRank (const Field& F,size_t n, int iterations)
{

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing rank"),"testRank",(unsigned int)iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	Element tmp;
	unsigned int r;
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations; ++k) {

		mycommentator().progress(k);
		BlasMatrix<Field> A(F,n,n),S(F,n,n), L(F,n,n);

		r = (unsigned int)((size_t)random() % n);
		// create S as an upper triangular matrix with r nonzero rows
		for (size_t i=0;i<r;++i){
			S.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<n;++j)
				S.setEntry(i,j,G.random(tmp));
		}
                BMD.write(commentator().report(), S) << std::endl;


		// create L as a lower triangular matrix with nonzero elements on the diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,Gn.random(tmp));
		}
                BMD.write(commentator().report(), L) << std::endl;

		//  compute A=LS
		BMD.mul(A,L,S);
                BMD.write(commentator().report(), A) << std::endl;

		// compute the rank of A
		unsigned int rank= BMD.rankin(A);
		commentator().report() << "Rank " << rank << " should be " << r << std::endl;

		if (rank!=r)
			ret=false;
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testRank");

	return ret;
}


/*
 *  Testing the determinant of dense matrices using BlasDomain
 *  construct a n*n matrices of determinant d and compute the determinant
 */
template <class Field>
static bool testDet (const Field& F,size_t n, int iterations)
{

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing determinant"),"testDet",(unsigned int)iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	Element tmp,One,d;
	F.init(One,1UL);

	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {

		mycommentator().progress(k);

		G.random(d);

		BlasMatrix<Field> A(F,n,n),S(F,n,n), L(F,n,n);

		// create S as an upper triangular matrix of full rank
		// with diagonal's element equal to 1 except the first entry wich equals to d
		for (size_t i=0;i<n;++i){
			S.setEntry(i,i,One);
			for (size_t j=i+1;j<n;++j)
				S.setEntry(i,j,G.random(tmp));
		}
		S.setEntry(0,0,d);

		// create L as a lower triangular matrix with only 1's on diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,One);
		}


		//  compute A=LS
		BMD.mul(A,L,S);

		// compute the determinant of A
		Element det= BMD.detin(A);

		if (!F.areEqual(det,d))
			ret=false;
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testDet");

	return ret;
}

/*
 *  Testing the inverse of dense matrices using BlasDomain
 *  construct a non-singular n*n matrices
 */
template <class Field>
static bool testInv (const Field& F,size_t n, int iterations)
{

	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef  BlasMatrix<Field> Matrix;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing inverse"),"testInv",(unsigned int)iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	Element One,tmp;
	F.init(One,1UL);

	bool ret = true;
	MatrixDomain<Field> MD(F);
	BlasMatrixDomain<Field> BMD(F);

	Matrix Id(F, n,n);
	for (size_t i=0;i<n;++i)
		Id.setEntry(i,i,One);

	for (int k=0;k<iterations;++k) {

		mycommentator().progress(k);


		Matrix A(F, n,n),S(F, n,n), L(F, n,n), invA(F, n,n);

		// create S as an upper triangular matrix of full rank
		// with nonzero random diagonal's element
		for (size_t i=0;i<n;++i){
			S.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<n;++j)
				S.setEntry(i,j,G.random(tmp));
		}

		// create L as a lower triangular matrix
		// with only 1's on diagonal
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,One);
		}

		//  compute A=LS
		BMD.mul(A,L,S);

		// compute the inverse of A
		BMD.inv(invA,A);

		// compute Ainv*A and A*Ainv
		BMD.mul(L,invA,A);
		BMD.mul(S,A,invA);

		if (!MD.areEqual(L,Id) || !MD.areEqual(S,Id))
			ret=false;
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testInv");

	return ret;
}


/*
 * Test resolution of linear system with a triangular matrix
 */
template <class Field>
static bool testTriangularSolve (const Field& F, size_t m, size_t n, int iterations)
{

	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Field>                       Matrix;
	typedef TriangularBlasMatrix<Field>   TriangularMatrix;
	typedef typename Field::RandIter                RandIter;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing triangular solver"),"testTriangularSolve",(unsigned int)iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	Element One,tmp;
	F.init(One,1UL);

	bool ret = true;
	MatrixDomain<Field> MD(F);
	VectorDomain<Field>  VD(F);
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {

		mycommentator().progress(k);

		Matrix Al(F, m,m),Au(F, m,m);
		Matrix X(F, m,n), B(F, m,n), C(F, m,n);

		std::vector<Element> b(m),x(m),c(m);

		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		//create random vector b
		for( size_t i=0;i<m;++i)
			F.init(b[i],G.random(tmp));

		// Create Au a random full rank upper triangular matrix
		for (size_t i=0;i<m;++i){
			Au.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<m;++j)
				Au.setEntry(i,j,G.random(tmp));
		}

		// Create Al a random full rank lower triangular matrix
		for (size_t i=0;i<m;++i){
			for (size_t j=0;j<i;++j)
				Al.setEntry(i,j,G.random(tmp));
			Al.setEntry(i,i,Gn.random(tmp));
		}

		// Create 2 trinagular matrix as view of matrix
		TriangularMatrix TAl(Al,LinBoxTag::Lower,LinBoxTag::NonUnit), TAu(Au,LinBoxTag::Upper,LinBoxTag::NonUnit);

		// testing solver with matrix right hand side
		BMD.left_solve(X,TAl,B);
		BMD.mul(C,Al,X);
		if (!MD.areEqual(C,B))
			ret=false;

		BMD.left_solve(X,TAu,B);
		BMD.mul(C,Au,X);
		if (!MD.areEqual(C,B))
			ret=false;

		// testing solver with matrix left hand side
		BMD.right_solve(X,TAl,B);
		BMD.mul(C,X,Al);
		if (!MD.areEqual(C,B))
			ret=false;

		BMD.right_solve(X,TAu,B);
		BMD.mul(C,X,Au);
		if (!MD.areEqual(C,B))
			ret=false;


		// testing solver with vector right hand side
		BMD.left_solve(x,TAl,b);
		BMD.mul(c,Al,x);
		if (!VD.areEqual(c,b))
			ret=false;

		BMD.left_solve(x,TAu,b);
		BMD.mul(c,Au,x);
		if (!VD.areEqual(c,b))
			ret=false;

		// testing solver with vector left hand side
		BMD.right_solve(x,TAl,b);
		BMD.mul(c,x,Al);
		if (!VD.areEqual(c,b))
			ret=false;

		BMD.right_solve(x,TAu,b);
		BMD.mul(c,x,Au);
		if (!VD.areEqual(c,b))
			ret=false;

	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testTriangularSolve");

	return ret;
}

/*
 * Test resolution of linear system with a matrix
 */
template <class Field>
static bool testSolve (const Field& F, size_t m, size_t n, int iterations)
{

	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Field>                       Matrix;
	typedef typename Field::RandIter                RandIter;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing solver"),"testTriangularSolve",(unsigned int)iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	Element One,tmp;
	F.init(One,1UL);

	bool ret = true;
	MatrixDomain<Field> MD(F);
	VectorDomain<Field>  VD(F);
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {

		mycommentator().progress(k);

		Matrix A(F, m,m),L(F, m,m),S(F, m,m);
		Matrix X(F, m,n), B(F, m,n), C(F, m,n);

		std::vector<Element> b(m),x(m),c(m);

		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<n;++j)
				B.setEntry(i,j,G.random(tmp));

		//create a random vector b
		for( size_t i=0;i<m;++i)
			F.init(b[i],G.random(tmp));


		// create S as an upper triangular matrix of full rank
		// with nonzero random diagonal's element
		for (size_t i=0;i<m;++i){
			S.setEntry(i,i,Gn.random(tmp));
			for (size_t j=i+1;j<n;++j)
				S.setEntry(i,j,G.random(tmp));
		}

		// create L as a lower triangular matrix
		// with only 1's on diagonal
		for (size_t i=0;i<m;++i){
			for (size_t j=0;j<i;++j)
				L.setEntry(i,j,G.random(tmp));
			L.setEntry(i,i,One);
		}

		//  compute A=LS
		BMD.mul(A,L,S);


		// testing solver with matrix right hand side
		BMD.left_solve(X,A,B);
		BMD.mul(C,A,X);
		if (!MD.areEqual(C,B))
			ret=false;

		BMD.left_solve(X,A,B);
		BMD.mul(C,A,X);
		if (!MD.areEqual(C,B))
			ret=false;

		// testing solver with matrix left hand side
		BMD.right_solve(X,A,B);
		BMD.mul(C,X,A);
		if (!MD.areEqual(C,B))
			ret=false;

		BMD.right_solve(X,A,B);
		BMD.mul(C,X,A);
		if (!MD.areEqual(C,B))
			ret=false;


		// testing solver with vector right hand side
		BMD.left_solve(x,A,b);
		BMD.mul(c,A,x);
		if (!VD.areEqual(c,b))
			ret=false;

		BMD.left_solve(x,A,b);
		BMD.mul(c,A,x);
		if (!VD.areEqual(c,b))
			ret=false;

		// testing solver with vector left hand side
		BMD.right_solve(x,A,b);
		BMD.mul(c,x,A);
		if (!VD.areEqual(c,b))
			ret=false;

		BMD.right_solve(x,A,b);
		BMD.mul(c,x,A);
		if (!VD.areEqual(c,b))
			ret=false;

	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testTriangularSolve");

	return ret;
}

/*
 * Test of the BlasPermutations
 */
template <class Field>
static bool testPermutation (const Field& F, size_t m, int iterations)
{

	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Field>                       Matrix;
	typedef typename Field::RandIter                RandIter;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing permutations"),"testPermutation",(unsigned int)iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	Element One,zero,tmp;
	F.init(One,1UL);
	F.init(zero,0UL);

	bool ret = true;
	MatrixDomain<Field> MD(F);
	VectorDomain<Field>  VD(F);
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {

		mycommentator().progress(k);

		std::vector<size_t> P(m);

		// Field Z2(2);
		// RandIter G2(Z2);

		// for (size_t i=0; i<m; ++i){
			// G.random(tmp);
			// if ( Z2.isZero(G2.random(tmp2) ) )
				// P[i] = i + ( (size_t) random() % (m-i) );
			// else
				// P[i] = i;
		// }

		//std::cerr<<P<<std::endl;
		Matrix A(F, m,m), Abis(F, m,m), B(F, m,m), C(F, m,m), D(F, m,m);
		std::vector<Element> a(m),abis(m),b(m),c(m), d(m);
		BlasPermutation<size_t>  Perm(m);
		RandomBlasPermutation(Perm);

		// Create A a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				A.setEntry(i,j,Gn.random(tmp));
		// Create a a random vector
		for (size_t i=0;i<m;++i)
			F.assign(a[i],Gn.random(tmp));

		/*
		 * Test A.P.P^t == A
		 */

		// B = A.P
		BMD.mul( B, A, Perm);
		// C = B.P^t
		BMD.mul( C, B, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm) );
		// Test C==A
		if (!MD.areEqual(A,C))
			ret=false;
		/*
		 * Test A.P^t.P == A
		 */

		// B = A.P^t
		BMD.mul( B, A, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm));
		// C = B.P
		BMD.mul( C, B, Perm );
		// Test C==A
		if (!MD.areEqual(A,C))
			ret=false;
		/*
		 * Test P.P^t.A == A
		 */

		// B = P.A
		BMD.mul( B, Perm, A);
		// C = P^t.B
		BMD.mul( C, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm) , B);
		// Test C==A
		if (!MD.areEqual(A,C))
			ret=false;
		/*
		 * Test P^t.P.A == A
		 */

		// B = P^t.A
		BMD.mul( B, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm), A);
		// C = P.B
		BMD.mul( C, Perm, B);
		// Test C==A
		if (!MD.areEqual(A,C))
			ret=false;

		/*
		 * Test a.P.P^t == a
		 */

		// b = a.P
		BMD.mul( b, a, Perm);
		// c = b.P^t
		BMD.mul( c, b, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm) );
		// Test c==a
		if (!VD.areEqual(a,c))
			ret=false;

		/*
		 * Test a.P^t.P == a
		 */

		// b = a.P^t
		BMD.mul( b, a, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm));
		// c = B.P
		BMD.mul( c, b, Perm );
		// Test c==a
		if (!VD.areEqual(a,c))
			ret=false;
		/*
		 * Test P.P^t.a == a
		 */

		// b = P.a
		BMD.mul( b, Perm, a);
		// c = P^t.b
		BMD.mul( c, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm) , b);
		// Test c==a
		if (!VD.areEqual(a,c))
			ret=false;

		/*
		 * Test P^t.P.a == a
		 */

		// b = P^t.a
		BMD.mul( b, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm), a);
		// c = P.b
		BMD.mul( c, Perm, b);
		// Test c==a
		if (!VD.areEqual(a,c))
			ret=false;


		/*
		 * Test P^t.A.(P.A)^-1.B == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		// Abis = P.A
		BMD.mul( Abis, Perm, A);
		// C = (P.A)^-1.B
		BMD.left_solve( C, Abis, B);
		// D = A.C (= P^-1.B)
		BMD.mul(D, A, C);
		// D = P.D
		BMD.mulin_right( Perm,D);
		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test A.P^t.(A.P)^-1.B == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		// Abis = A.P
		BMD.mul( Abis, A, Perm);
		// C = (A.P)^-1.B
		BMD.left_solve( C, Abis, B);
		// C = P.C
		BMD.mulin_right( Perm,C);
		// D = A.C (= P^-1.B)
		BMD.mul(D, A, C);

		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test B.P^t.A.(P.A)^-1 == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		// Abis = P.A
		BMD.mul( Abis, Perm, A);
		// C = B.(P.A)^-1
		BMD.right_solve( C, Abis, B);
		// C = C.P
		BMD.mulin_left( C,Perm);
		// D = C.A (=B)
		BMD.mul(D, C, A);
		if (!MD.areEqual(D,B))
		  ret=false;

		/*
		 * Test B.A.P^t.(A.P)^-1 == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		// Abis = A.P
		BMD.mul( Abis, A, Perm);
		// C = B.(A.P)^-1
		BMD.right_solve( C, Abis, B);
		// D = C.A (= B.P^t)
		BMD.mul(D, C, A);
		// C = C.P
		BMD.mulin_left( D, Perm);

		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test P.A.(P^t.A)^-1.B == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		// Abis = P^t.A
		BMD.mul( Abis, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm), A);
		// C = (P^t.A)^-1.B
		BMD.left_solve( C, Abis, B);
		// D = A.C (= P.B)
		BMD.mul(D, A, C);
		// D = P^t.D
		BMD.mulin_right( TransposedBlasMatrix<BlasPermutation<size_t> >(Perm),D);
		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test A.P.(A.P^t)^-1.B == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		// Abis = A.P^t
		BMD.mul( Abis, A, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm));
		// C = (A.P^t)^-1.B
		BMD.left_solve( C, Abis, B);
		// C = P^t.C
		BMD.mulin_right( TransposedBlasMatrix<BlasPermutation<size_t> >(Perm),C);
		// D = A.C (= P.B)
		BMD.mul(D, A, C);

		if (!MD.areEqual(D,B))
			ret=false;
		/*
		 * Test B.P.A.(P^t.A)^-1 == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		// Abis = P^t.A
		BMD.mul( Abis, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm), A);
		// C = B.(P^t.A)^-1
		BMD.right_solve( C, Abis, B);
		// C = C.P^t
		BMD.mulin_left( C,TransposedBlasMatrix<BlasPermutation<size_t> >(Perm));
		// D = C.A (=B)
		BMD.mul(D, C, A);
		if (!MD.areEqual(D,B))
		  ret=false;

		/*
		 * Test B.A.P.(A.P^t)^-1 == B
		 */
		// Create B a random matrix
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<m;++j)
				B.setEntry(i,j,G.random(tmp));

		// Abis = A.P^t
		BMD.mul( Abis, A, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm));
		// C = B.(A.P^t)^-1
		BMD.right_solve( C, Abis, B);
		// D = C.A (= B.P)
		BMD.mul(D, C, A);
		// C = C.P^t
		BMD.mulin_left( D, TransposedBlasMatrix<BlasPermutation<size_t> >(Perm));

		if (!MD.areEqual(D,B))
			ret=false;
	}
	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testLQUP");

	return ret;
}

/*
 * Test of the LQUPMatrix class
 */
template <class Field>
static bool testLQUP (const Field& F, size_t m, size_t n, int iterations)
{

	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Field>                       Matrix;
	typedef typename Field::RandIter                RandIter;

	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing LQUP factorization"),"testLQUP",(unsigned int)iterations);

	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	Element One,zero,tmp;
	F.init(One,1UL);
	F.init(zero,0UL);

	bool ret = true;
	MatrixDomain<Field> MD(F);
	BlasMatrixDomain<Field> BMD(F);

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
			    B.setEntry(i,j,zero);
		// Create C a random matrix of rank n/2
		for (size_t i=0;i<m;++i)
			if ( i % 2 )
				for (size_t j=0;j<n;++j)
					C.setEntry(i,j,G.random(tmp));
			else
				for (size_t j=0;j<n;++j)
					C.setEntry(i,j,zero);

		// A = B*C
		BMD.mul(A, B, C);

		Abis = A;

		BlasPermutation<size_t>  P(A.coldim()),Q(A.rowdim());
		LQUPMatrix<Field> X(A,P,Q);

		TriangularBlasMatrix<Field> L(F,m,m,LinBoxTag::Lower,LinBoxTag::Unit);
		TriangularBlasMatrix<Field> U(F,m,n,LinBoxTag::Upper,LinBoxTag::NonUnit);
		X.getL(L);
		X.getU(U);
		P=X.getP();

		Q=X.getQ();

		// C = U*P
		BMD.mul( C, U, P);
		// C = Q*C
		BMD.mulin_right( Q, C);
		// A = L*C
		BMD.mul( A, L, C);

		if (!MD.areEqual(A,Abis))
			ret=false;

		// Second pass
		// A = B*C
		BMD.mul(A, B, C);

		Abis = A;

		LQUPMatrix<Field> Y(A,P,Q);

		TriangularBlasMatrix<Field> L2(F,m,m,LinBoxTag::Lower,LinBoxTag::Unit);
		TriangularBlasMatrix<Field> U2(F,m,n,LinBoxTag::Upper,LinBoxTag::NonUnit);
		Y.getL(L2);
		Y.getU(U2);
		P=Y.getP();

		Q=Y.getQ();

		// C = Q*U2
		BMD.mul( C,Q,U2);
		// C = Q*C
		BMD.mulin_left(  C,P);
		// A = L*C
		BMD.mul( A, L2, C);

		if (!MD.areEqual(A,Abis))
			ret=false;
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testLQUP");

	return ret;
}

template <class Field>
static bool testMinPoly (const Field& F, size_t n, int iterations)
{
	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Field>                       Matrix;
	typedef typename Field::RandIter                RandIter;
	typedef vector<Element>                       Polynomial;
	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing minpoly"),"testMinPoly",(unsigned int)iterations);
	Element tmp, one, zero,mOne;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	F.init(one, 1UL);
	F.init(zero, 0UL);
	F.neg(mOne, one);
	//F.neg( mOne, one);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {

		mycommentator().progress(k);

		Matrix A(F,n,n);
		Polynomial P;
		// Test MinPoly(In) = X-1
		for (size_t i=0;i<n;++i)
			A.setEntry(i,i,one);

		BMD.minpoly( P, A );

		if ( P.size() !=2 )
			ret = false;
		if ( !F.areEqual(P[0], mOne) )
			ret = false;
		if ( !F.areEqual(P[1], one) )
			ret = false;

		// Test MinPoly(a*In) = X-a
		G.random(tmp);

		for (size_t i=0;i<n;++i)
			A.setEntry(i,i,tmp);
		F.negin(tmp);
		BMD.minpoly( P, A );
		if ( P.size() !=2 )
			ret = false;
		if ( !F.areEqual(P[0], tmp) )
			ret = false;
		if ( !F.areEqual(P[1], one) )
			ret = false;

		for (size_t i=0;i<n-1;++i){
			for (size_t j=0; j<i+1; ++j)
				A.setEntry(i,j,zero);
			A.setEntry(i,i+1,one);
			if (i<n-2)
				for (size_t j=i+2; j<n; ++j)
					A.setEntry(i,j,zero);
		}
		for (size_t j=0;j<n;++j)
			A.setEntry(n-1,j,zero);

		BMD.minpoly( P, A );
		if ( P.size() !=n+1 )
			ret = false;
		for (size_t i=0; i<n;++i)
			if ( !F.areEqual(P[i], zero) )
				ret = false;
		if ( !F.areEqual(P[n], one) )
			ret = false;

	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testMinPoly");

	return ret;
}

template <class Field>
static bool testCharPoly (const Field& F, size_t n, int iterations)
{
	typedef typename Field::Element                  Element;
	typedef BlasMatrix<Field>                       Matrix;
	typedef typename Field::RandIter                RandIter;
	typedef vector<Element>                       Polynomial;
	//Commentator mycommentator;
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	mycommentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);
	mycommentator().start (pretty("Testing charpoly"),"testCharPoly",(unsigned int)iterations);
	Element tmp, one, zero,mOne;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	F.init(one, 1UL);
	F.init(zero, 0UL);
	F.neg(mOne, one);
	//F.neg( mOne, one);
	bool ret = true;
	BlasMatrixDomain<Field> BMD(F);

	for (int k=0;k<iterations;++k) {

		mycommentator().progress(k);

		Matrix A(F,n,n);
		list<Polynomial> P;
		// Test CharPoly(In) = (X-1)^n
		for (size_t i=0;i<n;++i){
			for (size_t j=0;j<i;++j)
				A.setEntry(i,j,zero);
			A.setEntry(i,i,one);
			for (size_t j=i+1;j<n;++j)
				A.setEntry(i,j,zero);
		}
		P.clear();
		BMD.charpoly( P, A );

		typename list<Polynomial>::const_iterator P_it = P.begin();
		while (P_it != P.end()){
			if ( P_it->size() !=2 )
				ret = false;
			if ( !F.areEqual(P_it->operator[](0), mOne) )
				ret = false;
			if ( !F.areEqual(P_it->operator[](1), one) )
				ret = false;

			P_it++;
		}

		// Test CharPoly(a*In) = X-a
		G.random(tmp);

		for (size_t i=0;i<n;++i)
			A.setEntry(i,i,tmp);
		F.negin(tmp);
		P.clear();
		BMD.charpoly( P, A );
		P_it = P.begin();
		while (P_it != P.end()){
			if ( P_it->size() !=2 )
				ret = false;
			if ( !F.areEqual(P_it->operator[](0), tmp) )
				ret = false;
			if ( !F.areEqual(P_it->operator[](1), one) )
			ret = false;
			P_it++;
		}
	}

	mycommentator().stop(MSG_STATUS (ret), (const char *) 0, "testCharPoly");

	return ret;
}

template<class Field>
static bool testBlasMatrixConstructors(const Field& Fld, size_t m, size_t n) {
	bool pass = true;
	typedef typename Field::Element Element;
	BlasMatrixDomain<Field> BMD(Fld);
	MatrixDomain<Field>      MD (Fld);
	//BlasMatrix<Field> A; // nowhere to go

	BlasMatrix<Field> B(Fld);
	B.resize(m, n, Fld.zero);
	// don't understand the variations on integer types for the indices...

	BlasMatrix<Field> C(Fld,m,n);
	pass = pass and MD.areEqual(B, C);

//	MatrixStream<Field> ms; ...
//	BlasMatrix<Field> D(ms);
//	pass = pass and MD.areEqual(B, D);

	ScalarMatrix<Field> Eo(Fld, n, Fld.zero);
	BlasMatrix<Field> E(Eo); // copy a bb
	pass = pass and MD.areEqual(B, E);

#if 0
	ScalarMatrix<Field> Fo(Fld, 2*m, n, Fld.zero);
        BlasMatrix F(Fo, m, 0, m, n) ; // copy subm of a bb
	pass = pass and MD.areEqual(B, F);

	BlasMatrix<Field> G(Eo, Fld); // other field?
	pass = pass and MD.areEqual(B, G);

	BlasMatrix<Field> H(B); // copy cstor
	pass = pass and MD.areEqual(B, H);

	std::vector<Element> v(m*n, Fld.zero);
	BlasMatrix<Field> I(Fld,v,m,n);
	pass = pass and MD.areEqual(B, I);

	Element *p = &v[0];
	BlasMatrix<Field> J(Fld,p,m,n);
	pass = pass and MD.areEqual(B, J);
#endif

	return pass;
}

// returns true if ok, false if not.
template<class Field>
int launch_tests(Field & F, size_t n, int iterations)
{
	bool pass = true ;
	//std::cout << "no blas tests for now" << std::endl;
	// no slow test while I work on io
	if (!testBlasMatrixConstructors(F, n, n))             pass=false;
	if (!testMulAdd (F,n,iterations))                     pass=false;
	if (!testMulAddAgain (F,n,iterations))                pass=false;
	size_t m = n+n/2 ; size_t k = 2*n+1 ;
	if (!testMulAddShapeTrans (F,n,m,k,iterations))       pass=false;
	if (!testMulAddShapeTrans (F,n,k,m,iterations))       pass=false;
	if (!testMulAddShapeTrans (F,m,n,k,iterations))       pass=false;
	if (!testMulAddShapeTrans (F,m,k,n,iterations))       pass=false;
	if (!testMulAddShapeTrans (F,k,n,m,iterations))       pass=false;
	if (!testMulAddShapeTrans (F,k,m,n,iterations))       pass=false;
	if (!testTriangMulShapeTrans<Field,true,true>   (F,m,n,iterations))     pass=false;
	if (!testTriangMulShapeTrans<Field,true,true>   (F,n,m,iterations))     pass=false;
	if (!testTriangMulShapeTrans<Field,false,true>  (F,m,n,iterations))     pass=false;
	if (!testTriangMulShapeTrans<Field,false,true>  (F,n,m,iterations))     pass=false;
	if (!testTriangMulShapeTrans<Field,true,false>  (F,m,n,iterations))     pass=false;
	if (!testTriangMulShapeTrans<Field,true,false>  (F,n,m,iterations))     pass=false;
	if (!testTriangMulShapeTrans<Field,false,false> (F,m,n,iterations))     pass=false;
	if (!testTriangMulShapeTrans<Field,false,false> (F,n,m,iterations))     pass=false;
 	if (!testRank (F, n, iterations))                     pass=false;
 	if (!testDet  (F, n, iterations))                     pass=false;
 	if (!testInv  (F, n, iterations))                     pass=false;
 	if (!testTriangularSolve (F,n,n,iterations))          pass=false;
 	if (!testSolve (F,n,n,iterations))                    pass=false;
 	if (!testPermutation (F,n,iterations))                pass=false;
 	if (!testLQUP (F,n,n,iterations))                     pass=false;
 	if (!testMinPoly (F,n,iterations))                    pass=false;
	if (!testCharPoly (F,n,iterations))                   pass=false;
	//
	//
	return pass ;

}

int main(int argc, char **argv)
{

	static size_t n = 40;
	static integer q = 1000003U;
	static int iterations = 2;

    static Argument args[] = {
        { 'n', "-n N", "Set dimension of test matrices to NxN", TYPE_INT,     &n },
        { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1]",  TYPE_INTEGER, &q },
        { 'i', "-i I", "Perform each test for I iterations",    TYPE_INT,     &iterations },
	END_OF_ARGUMENTS
    };

	parseArguments (argc, argv, args);

	typedef Modular<double> Field;
	//typedef Modular<int> Field;
	//typedef Modular<float> Field;

	Field F1 (q);
	ModularBalanced<double> F2(q);
	Modular<float> F3(2011);
	GF2 F4 ;
	GivaroZpz<Givaro::Unsigned32> F5(2001);
	Modular<bool> F6 ;
	// BB. (Blas)MatrixDomain are not very generic...
#ifdef __LINBOX_HAVE_NTL
	NTL::ZZ_p::init(NTL::to_ZZ((size_t)q));
	NTL_ZZ_p F7;
#endif


	bool pass = true;

	srand ((unsigned)time (NULL));


	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	commentator().start("BlasMatrixDomain test suite", "BlasMatrixDomain");

	pass &= launch_tests(F1,n,iterations);
	pass &= launch_tests(F2,n,iterations);
	pass &= launch_tests(F3,n,iterations);
#pragma message "#warning GF2 is not working at all -> working on m4ri"
	// pass &= launch_tests(F4,n,iterations);
	// pass &= launch_tests(F6,n,iterations);
#pragma message "#warning GivaroZpz is not working at all"
	// pass &= launch_tests(F5,n,iterations);
#ifdef __LINBOX_HAVE_NTL
#pragma message "#warning NTL_ZZp is not working at all"
	// pass &= launch_tests(F7,n,iterations);
#endif

	commentator().stop(MSG_STATUS (pass), (const char *) 0,"BlasMatrixDomain test suite");
	return pass ? 0 : -1;
}



// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

