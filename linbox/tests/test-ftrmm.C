/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/test-ftrmm.inl
 * Copyright (C) 2010 LinBox
 *
 * Written by Brice Boyer <Brice.Boyer@imag.fr>
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
 */


/*! @file  tests/test-ftrmm.C
 * @ingroup tests
 * @brief  tests all variants of applyP, ftrmm, ftrsm and fgemm for varions m,k,n and leading dimensions combinations.
 * @bug ftrmm fails on non-double transposed versions.
 * @bug is fgemm correct on null A and B ?
 * @test FFLAS::ftrmm, FFPACK::applyP, FFLAS::ftrsm, FFLAS::fgemm
 */



//#define DEBUG

#include <cassert>
#include "linbox/field/modular.h"
#include "linbox/linbox-config.h"
#include "fflas-ffpack/fflas/fflas.h"
#include "fflas-ffpack/ffpack/ffpack.h"
#include "linbox/integer.h" // for fflas on integers...

//!@todo !
//#define __LINBOX_HAVE_INT64


#include "linbox/field/modular-balanced.h"
#include "linbox/field/modular.h"
#include "fflas-ffpack/utils/Matio.h"
#include "test-common.h"

//#define _LB_TIME
#define _LB_MAX_SZ 50
#define _LB_ITERS 3
//#define _LB_DEBUG


using namespace LinBox ;

template<class Field, FFLAS::FFLAS_SIDE Side,FFLAS::FFLAS_UPLO UpLo, FFLAS::FFLAS_TRANSPOSE Trans, FFLAS::FFLAS_DIAG Diag >
int test_ftrmm(std::ostream & report, const Field & F)
{
	linbox_check(_LB_MAX_SZ>3);
	size_t M    = random()%(_LB_MAX_SZ/2) ;
	size_t N    = random()%(_LB_MAX_SZ/2) ; // B is £MxN£ in a £ldb x rows£ table
	size_t ldb  = random()%_LB_MAX_SZ ;
	size_t lda  = random()%_LB_MAX_SZ ;
	size_t K    ;                           // A is £KxK£
	while (ldb<N) ldb = random()%_LB_MAX_SZ; // £ldb >= N£
	if (Side == FFLAS::FflasRight){
		K = N ;
	}
	else {
		K = M ;
	}
	while (lda<K) lda = random()%_LB_MAX_SZ; // £ldba>= N£

	assert(N    <= ldb);
	assert(K    <= lda);

	size_t rows = std::max(N,M)+3; // number of rows in A and B as a big table. (it is still M or N as a matrix)

#ifdef _LB_DEBUG
	report << "#M x N :"<< M << 'x' << N << std::endl;
	report << "#lda x ldb :"<< lda << 'x' << ldb << std::endl;
#endif
	assert(rows >= M);

	typedef typename Field::Element Element;

	Element * A = new Element[rows*lda];
	assert(A);
	Element * B = new Element[rows*ldb];
	assert(B);
	Element * C = new Element[M*N] ; // le résultat, le vrai, le bon.
	assert(C);
	Element * D = new Element[rows*ldb] ; // backup de B
	assert(D);

	typedef typename Field::RandIter RandIter;
	report << '#' ;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);

	/* init A,B,C and more...*/
	for (size_t i = 0 ; i < rows*lda ; ++i) G.random( *(A+i) ) ;
#ifdef _LB_DEBUG
	Element * E = new Element[rows*lda]; // copy of A
	assert(E);
	FFLAS::fcopy(F,rows*lda,E,1,A,1);
	// for (size_t i = 0 ; i < rows*lda ; ++i) *(E+i) = *(A+i);
#endif

	for (size_t i = 0 ; i < rows*ldb  ; ++i) G.random( *(B+i) ) ;

	Element zero ; F.init(zero,0UL);
	for (size_t i = 0 ; i < M*N ; ++i) *(C+i) = zero;
	// for (size_t i = 0 ; i < rows*ldb ; ++i) *(D+i) = *(B+i);
	FFLAS::fcopy(F,rows*ldb,D,1,B,1);
	Element alpha ;
	//! @todo F.isInvertible()
	//! @todo InvertibleRandomIter
	Element invalpha ;
	F.init(invalpha,0UL);
	Gn.random(alpha);
	// F.init(alpha,1UL);
	F.inv(invalpha,alpha);

	/* *************** */
	/*  testing ftrmm  */
	/* *************** */
	report << "# testing ftrmm A = (" << alpha << ") " << ((Side == FFLAS::FflasLeft)?(""):("A")) ;
	report << ((UpLo==FFLAS::FflasUpper)?(" U"):(" L")) ;
	report << ((Trans==FFLAS::FflasTrans)?(""):("^t")) ;
	report << ((Side == FFLAS::FflasLeft)?("A"):("")) ;
	report << " on field \""  ;
	F.write(report) ;
	report << "\" "<< ((Diag==FFLAS::FflasUnit)?("and the diagonal is unit"):("")) << std::endl;


	/*  compute C pedestrially :) */
	if (Side == FFLAS::FflasRight) // right
	{
		if (UpLo == FFLAS::FflasLower)
			if (Trans == FFLAS::FflasNoTrans) // C = M * L
				for (size_t i = 0 ; i < M ; ++i)
					for (size_t j = 0 ; j < N ; ++j)
						for (size_t k = j+1 ; k < N ; ++k)
							F.axpyin(*(C+i*N+j),*(B+i*ldb+k),*(A+k*lda+j));
			else // C = M * L^t
				for (size_t i = 0 ; i < M ; ++i)
					for (size_t j = 0 ; j < N ; ++j)
						for (size_t k = 0 ;  k < j ; ++k)
							F.axpyin(*(C+i*N+j),*(B+i*ldb+k),*(A+j*lda+k));
		else
			if (Trans == FFLAS::FflasNoTrans) // C = M * U
				for (size_t i = 0 ; i < M ; ++i)
					for (size_t j = 0 ; j < N ; ++j)
						for (size_t k = 0 ; k < j ; ++k)
							F.axpyin(*(C+i*N+j),*(B+i*ldb+k),*(A+k*lda+j));
			else // C = M * U^t
				for (size_t i = 0 ; i < M ; ++i)
					for (size_t j = 0 ; j < N ; ++j)
						for (size_t k = j+1 ; k < N ; ++k)
							F.axpyin(*(C+i*N+j),*(B+i*ldb+k),*(A+j*lda+k));
	}
	else // left
	{
		if (UpLo == FFLAS::FflasLower)
			if (Trans == FFLAS::FflasNoTrans) // C = L*M
				for (size_t i = 0 ; i < M ; ++i)
					for (size_t j = 0 ; j < N ; ++j)
						for (size_t k = 0 ; k < i ; ++k)
							F.axpyin(*(C+i*N+j),*(A+i*lda+k),*(B+k*ldb+j));
			else // C = L^t * M
				for (size_t i = 0 ; i < M ; ++i)
					for (size_t j = 0 ; j < N ; ++j)
						for (size_t k = i+1 ; k < M ; ++k)
							F.axpyin(*(C+i*N+j),*(A+k*lda+i),*(B+k*ldb+j));
		else
			if (Trans == FFLAS::FflasNoTrans) // C = U * M
				for (size_t i = 0 ; i < M ; ++i)
					for (size_t j = 0 ; j < N ; ++j)
						for (size_t k = i+1 ; k < M ; ++k)
							F.axpyin(*(C+i*N+j),*(A+i*lda+k),*(B+k*ldb+j));
			else // C = U^t M
				for (size_t i = 0 ; i < M ; ++i)
					for (size_t j = 0 ; j < N ; ++j)
						for (size_t k = 0 ; k < i ; ++k)
							F.axpyin(*(C+i*N+j),*(A+k*lda+i),*(B+k*ldb+j));
	}

	if (Diag == FFLAS::FflasUnit) // unit
		if (Side == FFLAS::FflasRight) // right
		{
			if (UpLo == FFLAS::FflasLower)
				if (Trans == FFLAS::FflasNoTrans) // C = M * L
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.addin(*(C+i*N+j),*(B+i*ldb+j));
				else // C = M * L^t
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.addin(*(C+i*N+j),*(B+i*ldb+j));
			else
				if (Trans == FFLAS::FflasNoTrans) // C = M * U
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.addin(*(C+i*N+j),*(B+i*ldb+j));
				else // C = M * U^t
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.addin(*(C+i*N+j),*(B+i*ldb+j));
		}
		else // left
		{
			if (UpLo == FFLAS::FflasLower)
				if (Trans == FFLAS::FflasNoTrans) // C = L*M
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.addin(*(C+i*N+j),*(B+i*ldb+j));
				else // C = L^t * M
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.addin(*(C+i*N+j),*(B+i*ldb+j));
			else
				if (Trans == FFLAS::FflasNoTrans) // C = U * M
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.addin(*(C+i*N+j),*(B+i*ldb+j));
				else // C = U^T M
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.addin(*(C+i*N+j),*(B+i*ldb+j));
		}
	else // non unit
		if (Side == FFLAS::FflasRight) // right
		{
			if (UpLo == FFLAS::FflasLower)
				if (Trans == FFLAS::FflasNoTrans) // C = M * L
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.axpyin(*(C+i*N+j),*(B+i*ldb+j),*(A+j*lda+j));
				else // C = M * L^t
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.axpyin(*(C+i*N+j),*(B+i*ldb+j),*(A+j*lda+j));
			else
				if (Trans == FFLAS::FflasNoTrans) // C = M * U
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.axpyin(*(C+i*N+j),*(B+i*ldb+j),*(A+j*lda+j));
				else // C = M * U^t
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.axpyin(*(C+i*N+j),*(B+i*ldb+j),*(A+j*lda+j));
		}
		else // left
		{
			if (UpLo == FFLAS::FflasLower)
				if (Trans == FFLAS::FflasNoTrans) // C = L*M
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.axpyin(*(C+i*N+j),*(A+i*lda+i),*(B+i*ldb+j));
				else // C = L^t * M
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.axpyin(*(C+i*N+j),*(A+i*lda+i),*(B+i*ldb+j));
			else
				if (Trans == FFLAS::FflasNoTrans) // C = U * M
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.axpyin(*(C+i*N+j),*(A+i*lda+i),*(B+i*ldb+j));
				else // C = M * U^t
					for (size_t i = 0 ; i < M ; ++i)
						for (size_t j = 0 ; j < N ; ++j)
							F.axpyin(*(C+i*N+j),*(A+i*lda+i),*(B+i*ldb+j));
		}

	for (size_t i = 0 ; i < M*N ; ++i)
		F.mulin(*(C+i),alpha);

	/*  compute B with ftrmm */
	FFLAS::ftrmm(F,Side,UpLo,Trans,Diag,M,N,alpha,A,lda,B,ldb);

	/*  checking C==B */
	int err = 0 ;
	for (size_t i = 0 ; i < M && !err ; ++i)
		for (size_t j = 0 ; j < N && !err ; ++j)
			if (!F.areEqual(*(C+i*N+j),*(B+i*ldb+j)))
				err = -1  ;
	if (err)
		report << "# *** error *** ftrmm fails" << std::endl;
	else {
	// checking B has nothing written outside B...
	for (size_t i = M ; i < rows && !err ; ++i)
		for (size_t j = N ; j < ldb  && !err ; ++j)
			if (!F.areEqual(*(D+i*ldb+j),*(B+i*ldb+j)))
				err = -1  ;
	if (err)
		report << "# *** error *** ftrmm writes outside : failure." << std::endl;

	}
#ifdef _LB_DEBUG
	{
		if (err)
		{
			report<<"#-------------T-------------" <<std::endl;
			report << "T :=" ;
			write_field(F,report,UpLo,Diag,E,rows,lda,lda,true);
			report << ':' << std::endl;
			report<<"#------------TT--------------" <<std::endl;
			report << "TT := " ; write_field(F,report,E,rows,lda,lda,true);
			report << ':' << std::endl;

			report<<"#-------------M-------------" <<std::endl;
			report << "M :=" ;
			write_field(F,report,D,M,N,ldb,true);
			report << ':' << std::endl;
			report<<"#------------a--------------" <<std::endl;
			report << "alpha := " << alpha << ':' << std::endl;
			report<<"#------------C--------------" <<std::endl;
			report << "C := " ; write_field(F,report,C,M,N,N,true);
			report << ':' << std::endl;
			report << "#------------B--------------" << std::endl;
			report << "B := " ; write_field(F,report,B,M,N,ldb,true);
			report << ':' << std::endl;
			report << "N := alpha * " ;
			if (Side == FFLAS::FflasRight)
				report << "M." ;
			else
				if (Trans == FFLAS::FflasNoTrans)
					report << "T." ;
				else
					report << "LinearAlgebra:-Transpose(T)." ;
			if (Side == FFLAS::FflasLeft)
				report << "M" ;
			else
				if (Trans == FFLAS::FflasNoTrans)
					report << "T" ;
				else
					report << "LinearAlgebra:-Transpose(T)" ;
			report << "  mod " << F.characteristic() << ':' << std::endl;
			report << "linalg:-iszero(C - N  mod " << F.characteristic() << "),";
			report << "linalg:-iszero(B - N  mod " << F.characteristic() << ");"  <<  std::endl;

		}

		}
#endif
	if (err)
		report << "# \033[1;31m>\033[0m ftrmm fail" ;
	else
		report << "# \033[1;32m>\033[0m ftrmm success";
	report << std::endl;


	/* *************** */
	/*  testing ftrsm  */
	/* *************** */
	report << "# testing ftrsm A = (" << alpha << ") " << ((Side == FFLAS::FflasLeft)?(""):("A")) ;
	report << ((UpLo==FFLAS::FflasUpper)?(" U^(-1)"):(" L^(-1)")) ;
	report << ((Trans==FFLAS::FflasTrans)?(""):("^t")) ;
	report << ((Side == FFLAS::FflasLeft)?("A"):("")) ;
	report << " on field \""  ;
	F.write(report) ;
	report << "\" "<< ((Diag==FFLAS::FflasUnit)?("and the diagonal is unit"):("")) << std::endl;

	int eur = 0 ;

	// for (size_t i = 0 ; i < rows*ldb ; ++i) *(B+i) = *(D+i);
	FFLAS::fcopy(F,rows*ldb,B,1,D,1);
	if (Diag == FFLAS::FflasNonUnit)
		for (size_t i = 0 ; i < K ; ++i) Gn.random(*(A+i*(lda+1))) ; // invertible diag !
	FFLAS::ftrmm(F, Side, UpLo, Trans, Diag, M, N, alpha,    A, lda, B, ldb);
	/* revert with ftrsm  */
	FFLAS::ftrsm(F, Side, UpLo, Trans, Diag, M, N, invalpha, A, lda, B, ldb);
	//! @todo check ftrsm fails nicely with non invertible A !

	for (size_t i = 0 ; i < rows && !eur ; ++i)
		for (size_t j = 0 ; j < ldb && !eur ; ++j)
			if (!F.areEqual(*(D+i*ldb+j),*(B+i*ldb+j)))
				eur = -1  ;

#ifdef _LB_DEBUG
	{
		delete[] E ;
	}
#endif
	if (eur)
		report << "# \033[1;31m>\033[0m ftrsm fail" ;
	else
		report << "# \033[1;32m>\033[0m ftrsm success";
	report << std::endl;


	delete[] A;
	delete[] B;
	delete[] C;
	delete[] D;

	return err+eur ;
}

//!@todo  test \c NULL permutation
template<class Field>
int test_applyP(std::ostream & report, const Field & F)
{
	size_t M    = random()%_LB_MAX_SZ+1 ;
	size_t N    = random()%_LB_MAX_SZ+1 ;
	size_t lda  = random()%_LB_MAX_SZ+1 ;
	if (lda<N) std::swap(lda,N);

#ifdef _LB_DEBUG
	report << "#M x N :"<< M << 'x' << N << std::endl;
	report << "#lda :"<< lda << std::endl;
#endif
	typedef typename Field::Element Element;

	Element * A = new Element[M*lda];
	assert(A);
	Element * B = new Element[M*lda];
	assert(B);

	typedef typename Field::RandIter RandIter;
	report << '#' ;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);

	/* init A,B,C and more...*/
	for (size_t i = 0 ; i < M*lda ; ++i) G.random( *(A+i) ) ;
	for (size_t i = 0 ; i < M*lda ; ++i)  *(B+i) = *(A+i) ;

	report << "# testing row applyP" << std::endl;
	size_t * P = new size_t[M] ;

	size_t r = random()%M ;
	for (size_t i = 0 ; i < r ; ++i){
		P[i] = i + (size_t)(M-i)*( (std::max(0.,drand48()-1e-6)) ); // histoire de pas tomber sur 1...
		//                if (P[i] == M) P[i]-- ;
		assert(P[i] < M);
	}

	FFPACK::applyP(F,FFLAS::FflasLeft,FFLAS::FflasNoTrans,
		       N, 0,r,
		       A,lda,P);

	FFPACK::applyP(F,FFLAS::FflasLeft,FFLAS::FflasTrans,
		       N, 0,r,
		       A,lda,P);
	delete[] P ;
	int err = 0 ;
	for (size_t i = 0 ; i < M && !err ; ++i)
		for (size_t j = 0 ; j < N && !err ; ++j)
			if (!F.areEqual(*(A+i*lda+j),*(B+i*lda+j)))
				err = -1  ;
	if (err)
		report << "# \033[1;31m>\033[0m row applyP failed" << std::endl;
	else
		report << "# \033[1;32m>\033[0m row applyP success" << std::endl;

	report << "# testing col applyP" << std::endl;
	size_t * Q = new size_t[N] ;

	r = random()%N ;
	for (size_t i = 0 ; i < r ; ++i){
		Q[i] = i + (size_t)(N-i)*( std::max(0.,drand48()-1e-6)) ;
		//                if (Q[i] == N) P[i]-- ;
		assert(Q[i] < N);
	}

	FFPACK::applyP(F,FFLAS::FflasRight,FFLAS::FflasNoTrans,
		       M, 0,r,
		       A,lda,Q);

	FFPACK::applyP(F,FFLAS::FflasRight,FFLAS::FflasTrans,
		       M, 0,r,
		       A,lda,Q);
	delete[] Q ;

	int eur = 0 ;
	for (size_t i = 0 ; i < M && !eur ; ++i)
		for (size_t j = 0 ; j < N && !eur ; ++j)
			if (!F.areEqual(*(A+i*lda+j),*(B+i*lda+j)))
				err = -1  ;
	if (eur)
		report << "# \033[1;31m>\033[0mcol applyP failed" << std::endl;
	else
		report << "# \033[1;32m>\033[0mcol applyP success" << std::endl;

	delete[] A;
	delete[] B ;
	return (err+eur);

}

template<class Field, FFLAS::FFLAS_TRANSPOSE ATRANS, FFLAS::FFLAS_TRANSPOSE BTRANS>
int test_fgemm(std::ostream & report, const Field & F)
{
	size_t M    = random()%_LB_MAX_SZ+1 ; // rows of (t)A and C
	size_t N    = random()%_LB_MAX_SZ+1 ; // cols of (t)B and C
	size_t K    = random()%_LB_MAX_SZ+1 ; // cols of (t)A, rows of (t)B

	//!@bug if beta != 0 but A and B are 0 (or don't exist) then C != beta C

	size_t lda  = random()%(_LB_MAX_SZ/2) ;
	if (ATRANS == FFLAS::FflasTrans) lda += M ; else lda += K ;
	size_t ldb  = random()%(_LB_MAX_SZ/2) ;
	if (BTRANS == FFLAS::FflasTrans) ldb += K ; else ldb += N ;
	size_t ldc  = N+random()%(_LB_MAX_SZ/2) ;

	size_t rowA = 5;
	if (ATRANS == FFLAS::FflasTrans) rowA += K ; else rowA += M ;
	size_t rowB = 5;
	if (BTRANS == FFLAS::FflasTrans) rowB += N ; else rowB += K ;
	size_t rowC = M+5;


	// A is M x K
#ifdef _LB_DEBUG
	report <<"# A : " << M << 'x' << K << " ("<< rowA << 'x' << lda << ')' << std::endl;
	report <<"# B : " << K << 'x' << N << " ("<< rowB << 'x' << ldb << ')' << std::endl;
	report <<"# C : " << M << 'x' << N << " ("<< rowC << 'x' << ldc << ')' << std::endl;
#endif

	typedef typename Field::RandIter RandIter;
	report << "#" ;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G);
	typedef typename Field::Element Element;

	Element alpha, beta ;
	// G.random(alpha);
	// G.random(beta);
	alpha = Integer::random<false>(2);
	if (abs(alpha) > 1.5 ) G.random(alpha);
	//!@bug needs p prime.
	beta  = Integer::random<false>(2);
	if (abs(beta) >1.5 ) G.random(beta);
	report << "#testing fgemm C = (" << alpha << ") A" << ((ATRANS==FFLAS::FflasTrans)?("^T"):("")) << ".B"  << ((BTRANS==FFLAS::FflasTrans)?("^T"):("")) << "+ (" << beta << ") C on " ; F.write(report) << std::endl;

	F.init(alpha,alpha);
	F.init(beta,beta);


	Element * A = new Element[rowA*lda];
	assert(A);
	Element * B = new Element[rowB*ldb];
	assert(B);
	Element * C = new Element[rowC*ldc] ; // le résultat, le vrai, le bon.
	assert(C);
	Element * D = new Element[rowC*ldc] ; // backup de C
	assert(D);

	for (size_t i = 0 ; i < rowA*lda ; ++i) G.random( *(A+i) ) ;
	for (size_t i = 0 ; i < rowB*ldb ; ++i) G.random( *(B+i) ) ;
	for (size_t i = 0 ; i < rowC*ldc ; ++i) G.random( *(C+i) ) ;
	for (size_t i = 0 ; i < rowC*ldc ; ++i) *(D+i) = *(C+i)  ;

	if (ATRANS == FFLAS::FflasNoTrans)
		if (BTRANS == FFLAS::FflasNoTrans) // A.B + C
			for (size_t i = 0 ; i < M ; ++i)
				for (size_t j = 0 ; j < N ; ++j) {
					Element temp;
					F.init(temp,0UL);
					for (size_t k = 0 ; k < K ; ++k)
						F.axpyin(temp,A[i*lda+k],B[k*ldb+j]);
					F.mulin(C[i*ldc+j],beta);
					F.axpyin(C[i*ldc+j],alpha,temp);
				}
		else // A.tB + C
			for (size_t i = 0 ; i < M ; ++i)
				for (size_t j = 0 ; j < N ; ++j) {
					Element temp;
					F.init(temp,0UL);
					for (size_t k = 0 ; k < K ; ++k)
						F.axpyin(temp,A[i*lda+k],B[j*ldb+k]);
					F.mulin(C[i*ldc+j],beta);
					F.axpyin(C[i*ldc+j],alpha,temp);
				}
	else
		if (BTRANS == FFLAS::FflasNoTrans) // tA.B + C
			for (size_t i = 0 ; i < M ; ++i)
				for (size_t j = 0 ; j < N ; ++j) {
					Element temp;
					F.init(temp,0UL);
					for (size_t k = 0 ; k < K ; ++k)
						F.axpyin(temp,A[k*lda+i],B[k*ldb+j]);
					F.mulin(C[i*ldc+j],beta);
					F.axpyin(C[i*ldc+j],alpha,temp);
				}
		else // tA.tB + C
			for (size_t i = 0 ; i < M ; ++i)
				for (size_t j = 0 ; j < N ; ++j) {
					Element temp;
					F.init(temp,0UL);
					for (size_t k = 0 ; k < K ; ++k)
						F.axpyin(temp,A[k*lda+i],B[j*ldb+k]);
					F.mulin(C[i*ldc+j],beta);
					F.axpyin(C[i*ldc+j],alpha,temp);
				}

	FFLAS::fgemm(F,ATRANS,BTRANS,M,N,K,alpha,A,lda,B,ldb,beta,D,ldc);

	int err = 0 ;
	for (size_t i = 0 ; i < rowC && !err ; ++i)
		for (size_t j = 0 ; j < ldc && !err ; ++j)
			if (!F.areEqual(*(C+i*ldc+j),*(D+i*ldc+j)))
				err = -1  ;

	delete[] A ;
	delete[] B ;
	delete[] C ;
	delete[] D ;

	if (err) {
		report << "# \033[1;31m>\033[0mfgemm failed  " << std::endl;
	}
	else{
		report << "# \033[1;32m>\033[0mfgemm success " << std::endl;
	}

	return err ;
}

int main(int ac, char ** av)
{
	//typedef ModularBalanced<float>  FieldF;
	typedef Modular<float>          FieldF;
	//typedef ModularBalanced<double> FieldD;
	typedef Modular<double>         FieldD;
	//typedef ModularBalanced<int32_t>  FieldI;
	typedef Modular<int32_t>          FieldI;
	//!@bug : this one completely fails :
	//typedef Modular<Integer>          FieldI;

       	static Argument as[] = {
		END_OF_ARGUMENTS
        };

	parseArguments (ac, av, as);

	srand(time(NULL));
	bool fail = false ;

	LinBox::commentator.start("ftrmm and consorts full test suite", "ftrmm et al");
	std::ostream &report = LinBox::commentator.report (LinBox::Commentator::LEVEL_IMPORTANT,
							   INTERNAL_DESCRIPTION);


	// need to be prime for ftrsm
	FieldD FD(13);
	FieldD FD2(65563);
	FieldF FF(13);
	FieldF FF2(1069);
	FieldI FI(13);
	FieldI FI2(106739);
	int tot = 6;
#ifdef __LINBOX_HAVE_INT64
	//        typedef ModularBalanced<int64_t>  FieldU;
	typedef Modular<int64_t>          FieldU;
	FieldU FU(13);
	FieldU FU2(13132153);
	tot += 2 ;
#endif
	tot *=2*16*_LB_ITERS ; // ftr(s/m)m 16 combinations, repeated _LB_ITERS times
	int ret = tot;
	//-----------//
	/* FTR(S/M)M */
	//-----------//
	for (size_t r = 0 ; r < _LB_ITERS ; ++r)
	{
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>(report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>(report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FD2);

		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FD2);

		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FD2);

		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FF);

		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FF);

		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FF);

		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FF);

		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FF2);

		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FF2);

		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FF2);

		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FF2);


		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FI);

		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FI);

		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FI);

		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FI);

		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FI2);

		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FI2);

		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FI2);

		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FI2);

#ifdef __LINBOX_HAVE_INT64
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FU);

		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FU);

		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FU);

		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FU);

		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FU2);

		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FU2);

		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FU2);

		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasUnit>   (report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,  FFLAS::FflasNonUnit>(report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>   (report,FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(report,FU2);
#endif
	}


#ifdef DEBUG
	report << "# \033[1;33m>\033[0m ftr(s/m)m  passed " << ret << "/" << tot << "tests" <<std::endl;
#endif
	if (ret != tot) fail=true;
#ifdef DEBUG
	if (fail)
		report << "# \033[1;31m>\033[0m ftr(s/m)m failed" << std::endl;
#endif
	int our = tot = 6*_LB_ITERS*2 ;
#ifdef __LINBOX_HAVE_INT64
	our = tot = tot+2*_LB_ITERS*2 ;
#endif

	//-------//
	/* APPLY */
	//-------//
	for (size_t r = 0 ; r < _LB_ITERS ; ++r)
	{
		our+= test_applyP(report,FD);
		our+= test_applyP(report,FD2);
		our+= test_applyP(report,FI);
		our+= test_applyP(report,FI2);
		our+= test_applyP(report,FF);
		our+= test_applyP(report,FF2);
#ifdef __LINBOX_HAVE_INT64
		our+= test_applyP(report,FU);
		our+= test_applyP(report,FU2);
#endif

	}
#ifdef DEBUG
	report << "# \033[1;33m>\033[0m applyP  passed " << our << "/" << tot << "tests" <<std::endl;
#endif
	if (our != tot) fail = true;
#ifdef DEBUG
	if (our != tot)
		report << "# \033[1;31m>\033[0m applyP failed" << std::endl;
#endif


	our = tot = 4*_LB_ITERS*6 ;
#ifdef __LINBOX_HAVE_INT64
	our = tot = tot+4*_LB_ITERS*2 ;
#endif

	//-------//
	/* FGEMM */
	//-------//
	for (size_t r = 0 ; r < _LB_ITERS ; ++r)
	{
		our+= test_fgemm<FieldD,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans>(report,FD);
		our+= test_fgemm<FieldD,FFLAS::FflasNoTrans, FFLAS::FflasTrans>  (report,FD);
		our+= test_fgemm<FieldD,FFLAS::FflasTrans,   FFLAS::FflasNoTrans>(report,FD);
		our+= test_fgemm<FieldD,FFLAS::FflasTrans,   FFLAS::FflasTrans>  (report,FD);

		our+= test_fgemm<FieldD,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans>(report,FD2);
		our+= test_fgemm<FieldD,FFLAS::FflasNoTrans, FFLAS::FflasTrans>  (report,FD2);
		our+= test_fgemm<FieldD,FFLAS::FflasTrans,   FFLAS::FflasNoTrans>(report,FD2);
		our+= test_fgemm<FieldD,FFLAS::FflasTrans,   FFLAS::FflasTrans>  (report,FD2);

		our+= test_fgemm<FieldI,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans>(report,FI);
		our+= test_fgemm<FieldI,FFLAS::FflasNoTrans, FFLAS::FflasTrans>  (report,FI);
		our+= test_fgemm<FieldI,FFLAS::FflasTrans,   FFLAS::FflasNoTrans>(report,FI);
		our+= test_fgemm<FieldI,FFLAS::FflasTrans,   FFLAS::FflasTrans>  (report,FI);

		our+= test_fgemm<FieldI,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans>(report,FI2);
		our+= test_fgemm<FieldI,FFLAS::FflasNoTrans, FFLAS::FflasTrans>  (report,FI2);
		our+= test_fgemm<FieldI,FFLAS::FflasTrans,   FFLAS::FflasNoTrans>(report,FI2);
		our+= test_fgemm<FieldI,FFLAS::FflasTrans,   FFLAS::FflasTrans>  (report,FI2);

		our+= test_fgemm<FieldF,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans>(report,FF);
		our+= test_fgemm<FieldF,FFLAS::FflasNoTrans, FFLAS::FflasTrans>  (report,FF);
		our+= test_fgemm<FieldF,FFLAS::FflasTrans,   FFLAS::FflasNoTrans>(report,FF);
		our+= test_fgemm<FieldF,FFLAS::FflasTrans,   FFLAS::FflasTrans>  (report,FF);

		our+= test_fgemm<FieldF,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans>(report,FF2);
		our+= test_fgemm<FieldF,FFLAS::FflasNoTrans, FFLAS::FflasTrans>  (report,FF2);
		our+= test_fgemm<FieldF,FFLAS::FflasTrans,   FFLAS::FflasNoTrans>(report,FF2);
		our+= test_fgemm<FieldF,FFLAS::FflasTrans,   FFLAS::FflasTrans>  (report,FF2);

#ifdef __LINBOX_HAVE_INT64
		our+= test_fgemm<FieldU,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans>(report,FU);
		our+= test_fgemm<FieldU,FFLAS::FflasNoTrans, FFLAS::FflasTrans>  (report,FU);
		our+= test_fgemm<FieldU,FFLAS::FflasTrans,   FFLAS::FflasNoTrans>(report,FU);
		our+= test_fgemm<FieldU,FFLAS::FflasTrans,   FFLAS::FflasTrans>  (report,FU);

		our+= test_fgemm<FieldU,FFLAS::FflasNoTrans, FFLAS::FflasNoTrans>(report,FU2);
		our+= test_fgemm<FieldU,FFLAS::FflasNoTrans, FFLAS::FflasTrans>  (report,FU2);
		our+= test_fgemm<FieldU,FFLAS::FflasTrans,   FFLAS::FflasNoTrans>(report,FU2);
		our+= test_fgemm<FieldU,FFLAS::FflasTrans,   FFLAS::FflasTrans>  (report,FU2);
#endif

	}
#ifdef DEBUG
	report << "# \033[1;33m>\033[0m fgemm  passed " << our << "/" << tot << "tests" <<std::endl;
#endif
	if(our != tot) fail = true;

	commentator.stop(MSG_STATUS (!fail), (const char *) 0,"ftrmm et al full tests suite");
#ifdef DEBUG
	if (our != tot)
		report << "# \033[1;31m>\033[0m fgemm failed" << std::endl;
#endif
	return fail ;
	// return false ;

}

#undef _LB_TIME
#undef _LB_MAX_SZ
#undef _LB_ITERS
#undef _LB_DEBUG

