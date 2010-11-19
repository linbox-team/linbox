/* fflas/test-ftrmm.inl
 * Copyright (C) 2010 LinBox
 *
 * Written by Brice Boyer <Brice.Boyer@imag.fr>
 *
 * See COPYING for license information.
 */

//#define DEBUG

#include "linbox/linbox-config.h"
#include "linbox/fflas/fflas.h"
#include "linbox/ffpack/ffpack.h"
#include "linbox/integer.h" // for fflas on integers...

//!@todo !
//#define __LINBOX_HAVE_INT64

#include "linbox/field/modular-balanced-int32.h"
#include "linbox/field/modular-int32.h"
#ifdef __LINBOX_HAVE_INT64
#include "linbox/field/modular-balanced-int64.h"
#endif
#include "linbox/field/modular-balanced-double.h"
#include "linbox/field/modular-double.h"
#include "linbox/field/modular-balanced-float.h"
#include "linbox/field/modular-float.h"
#include "Matio.h"

//#define _LB_TIME
#define _LB_MAX_SZ 50
#define _LB_ITERS 3 
//#define _LB_DEBUG


using namespace LinBox ;

template<class Field, FFLAS::FFLAS_SIDE Side,FFLAS::FFLAS_UPLO UpLo, FFLAS::FFLAS_TRANSPOSE Trans, FFLAS::FFLAS_DIAG Diag >
int test_ftrmm(const Field & F) 
{
	srand(time(NULL));
	size_t M    = random()%_LB_MAX_SZ ;
	size_t N    = random()%_LB_MAX_SZ ; // B is MxN in a ldb*rows table
	size_t ldb  = random()%_LB_MAX_SZ+1 ;
	size_t lda  = random()%_LB_MAX_SZ+1 ; 
	while (ldb<N || ldb<M) ldb = random()%_LB_MAX_SZ;

	size_t Ldim = M; // A is MxM or NxN in a lda x rows table
	if (Side == FFLAS::FflasRight) Ldim = N ; // A = N x lda, else A is M x lda.
	while (lda<Ldim)  lda  = random()%_LB_MAX_SZ ; 

	size_t rows = std::max(N,M)+3; 

#ifdef _LB_DEBUG
	std::cout << "#M x N :"<< M << 'x' << N << std::endl; 
	std::cout << "#lda x ldb :"<< lda << 'x' << ldb << std::endl; 
#endif
	assert(N    <= ldb);
	assert(M    <= ldb);
	assert(Ldim <= lda);
	assert(rows >= Ldim);
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
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 

	/* init A,B,C and more...*/
	for (size_t i = 0 ; i < rows*lda ; ++i) G.random( *(A+i) ) ;
#ifdef _LB_DEBUG
	Element * E = new Element[rows*lda];
	assert(E);
	for (size_t i = 0 ; i < rows*lda ; ++i) *(E+i) = *(A+i);
#endif

	for (size_t i = 0 ; i < rows*ldb  ; ++i) G.random( *(B+i) ) ;

	Element zero ; F.init(zero,0UL);
	for (size_t i = 0 ; i < M*N ; ++i) *(C+i) = zero;
	for (size_t i = 0 ; i < rows*ldb ; ++i) *(D+i) = *(B+i);
	Element alpha ;
	//! @todo F.isInvertible()
	//! @toto InvertibleRandomIter
	Element invalpha ;
	F.init(invalpha,0UL);
	Gn.random(alpha);
	F.inv(invalpha,alpha);
	/* *************** */
	/*  testing ftrmm  */
	/* *************** */

	/*  compute C pedestrially */
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

	/*  checking  */
	int err = 0 ;
	for (size_t i = 0 ; i < M && !err ; ++i)
		for (size_t j = 0 ; j < N && !err ; ++j)
			if (!F.areEqual(*(C+i*N+j),*(B+i*ldb+j))) 
				err = -1  ;
	// checking B has nothing written outside B...
	for (size_t i = M ; i < rows && !err ; ++i)
		for (size_t j = N ; j < ldb  && !err ; ++j)
			if (!F.areEqual(*(D+i*ldb+j),*(B+i*ldb+j))) 
				err = -1  ;
#ifdef _LB_DEBUG
	{
	if (err)
	{
		std::cout<<"#-------------T-------------" <<std::endl;
		std::cout << "T :=" ;
		write_field(F,std::cout,UpLo,Diag,E,Ldim,Ldim,lda,true,true);
		std::cout<<"#-------------M-------------" <<std::endl;
		std::cout << "M :=" ;
		write_field(F,std::cout,D,M,N,ldb,true,true);
		std::cout<<"#------------a--------------" <<std::endl;
		std::cout << "alpha := " << alpha << ':' << std::endl;
		std::cout<<"#------------C--------------" <<std::endl;
		std::cout << "C := " ; write_field(F,std::cout,C,M,N,N,true,true);
		std::cout << "#------------B--------------"                     << std::endl;
		std::cout << "B := " ; write_field(F,std::cout,B,M,N,ldb,true,true);
		std::cout << "N := alpha * " ; 
		if (Side == FFLAS::FflasRight)
			std::cout << "M." ;
		else
			if (Trans == FFLAS::FflasNoTrans)
				std::cout << "T." ;
			else
				std::cout << "LinearAlgebra:-Transpose(T)." ;
		if (Side == FFLAS::FflasLeft)
			std::cout << "M" ;
		else
			if (Trans == FFLAS::FflasNoTrans)
				std::cout << "T" ;
			else
				std::cout << "LinearAlgebra:-Transpose(T)" ;
		std:: cout << "  mod " << F.characteristic() << ':' << std::endl;
		std::cout << "linalg:-iszero(C - N  mod " << F.characteristic() << "),";
		std::cout << "linalg:-iszero(B - N  mod " << F.characteristic() << ");"  <<  std::endl;

	}

	if (err) std::cout << "# \033[1;31m>\033[0m ftrmm fail" ;
	else std::cout << "# \033[1;32m>\033[0m ftrmm success";

	std::cout << " with " << ((Side == FFLAS::FflasLeft)?("left"):("right")) ;
	std::cout << " with " ;
	std::cout << ((UpLo==FFLAS::FflasUpper)?("upper "):("lower ")) ;
	std::cout << ((Diag==FFLAS::FflasUnit)?(""):("non-")) ;
	std::cout << "unit " ;
	std::cout << ((Trans==FFLAS::FflasTrans)?(""):("non-")) ;
	std::cout << "transposed triangular matrices (on field \""  ;
	F.write(std::cout) ;
	std::cout << "\")" << std::endl;
	}
#endif
	/* *************** */
	/*  testing ftrsm  */
	/* *************** */
	int eur = 0 ;

	for (size_t i = 0 ; i < rows*ldb ; ++i) *(B+i) = *(D+i);
	if (Diag == FFLAS::FflasNonUnit)
		for (size_t i = 0 ; i < Ldim ; ++i) Gn.random(*(A+i*(lda+1))) ; // invertible diag !
	FFLAS::ftrmm(F, Side, UpLo, Trans, Diag, M, N, alpha,    A, lda, B, ldb);
	/* revert with ftrsm  */
	FFLAS::ftrsm(F, Side, UpLo, Trans, Diag, M, N, invalpha, A, lda, B, ldb);
	//! @todo check ftrsm fails nicely with non invertible A !

	for (size_t i = 0 ; i < M && !eur ; ++i)
		for (size_t j = 0 ; j < N && !eur ; ++j)
			if (!F.areEqual(*(D+i*ldb+j),*(B+i*ldb+j))) 
				eur = -1  ;

#ifdef _LB_DEBUG
	{
		if (eur) std::cout << "# \033[1;31m>\033[0m ftrsm fail" ;
		else std::cout << "# \033[1;32m>\033[0m ftrsm success";
		std::cout << " with " << ((Side == FFLAS::FflasLeft)?("left"):("right")) ;
		std::cout << " with " ;
		std::cout << ((UpLo==FFLAS::FflasUpper)?("upper "):("lower ")) ;
		std::cout << ((Diag==FFLAS::FflasUnit)?(""):("non-")) ;
		std::cout << "unit " ;
		std::cout << ((Trans==FFLAS::FflasTrans)?(""):("non-")) ;
		std::cout << "transposed triangular matrices (on field \""  ;
		F.write(std::cout) ;
		std::cout << "\")" << std::endl;
		delete[] E ;
	}
#endif

	delete[] A;
	delete[] B;
	delete[] C;
	delete[] D;

	return err+eur ;
}
//!@todo  test \c NULL permutation
template<class Field>
int test_applyP(const Field & F)
{
	srand(time(NULL));
	size_t M    = random()%_LB_MAX_SZ+1 ;
	size_t N    = random()%_LB_MAX_SZ+1 ; 
	size_t lda  = random()%_LB_MAX_SZ+1 ; 
	if (lda<N) std::swap(lda,N);

#ifdef _LB_DEBUG
	std::cout << "#M x N :"<< M << 'x' << N << std::endl; 
	std::cout << "#lda :"<< lda << std::endl; 
#endif
	typedef typename Field::Element Element;

	Element * A = new Element[M*lda];
	assert(A);
	Element * B = new Element[M*lda];
	assert(B);

	typedef typename Field::RandIter RandIter;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 

	/* init A,B,C and more...*/
	for (size_t i = 0 ; i < M*lda ; ++i) G.random( *(A+i) ) ;
	for (size_t i = 0 ; i < M*lda ; ++i)  *(B+i) = *(A+i) ;

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
		std::cout << "# \033[1;31m>\033[0m row applyP failed" << std::endl;

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
		std::cout << "# \033[1;31m>\033[0mcol applyP failed" << std::endl;

	delete[] A;
	delete[] B ;
	return (err+eur);

}


template<class Field, FFLAS::FFLAS_TRANSPOSE At, FFLAS::FFLAS_TRANSPOSE Bt>
int test_fgemm(const Field & F)
{
	srand(time(NULL));
	size_t M    = random()%_LB_MAX_SZ ;
	size_t N    = random()%_LB_MAX_SZ ; 
	size_t K    = random()%_LB_MAX_SZ ; 

	size_t lda  = random()%_LB_MAX_SZ+1 ; 
	size_t ldb  = random()%_LB_MAX_SZ+1 ; 
	size_t ldc  = random()%_LB_MAX_SZ+1 ; 

	// A is M x K
	// B is K x N
	// C is M x N
	
	size_t rowA = M;
	size_t rowB = K;
	size_t rowC = M;

	typedef typename Field::RandIter RandIter;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 
	typedef typename Field::Element Element;

	Element alpha, beta ;
	G.random(alpha);
	G.random(beta);




	Element * A = new Element[rowA*lda];
	assert(A);
	Element * B = new Element[rowB*ldb];
	assert(B);
	Element * C = new Element[rowC*ldc] ; // le résultat, le vrai, le bon.
	assert(C);
	Element * D = new Element[rowC*ldc] ; // backup de C
	assert(D);



	for (size_t i = 0 ; i < rowA*lda ; ++i) G.random( *(A+i) ) ;
	for (size_t i = 0 ; i < rowB*lda ; ++i) G.random( *(B+i) ) ;
	for (size_t i = 0 ; i < rowC*lda ; ++i) G.random( *(D+i) ) ;
	for (size_t i = 0 ; i < rowC*lda ; ++i) *(C+i) = *(D+i)  ;


	if (At == FFLAS::FflasNoTrans)
		if (Bt == FFLAS::FflasNoTrans) // A.B + C
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
		if (Bt == FFLAS::FflasNoTrans) // tA.B + C
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

	FFLAS::fgemm(F,At,Bt,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc);

	FFLAS::fgemm(F,At,Bt,M,N,K,alpha,A,lda,B,ldb,beta,C,ldc,1);
return 0 ;

}

int main() 
{
//        typedef ModularBalanced<float>  FieldF;
	typedef Modular<float>          FieldF;
//        typedef ModularBalanced<double> FieldD;
	typedef Modular<double>         FieldD;
//        typedef ModularBalanced<int32>  FieldI;
	typedef Modular<int32>          FieldI;
	//!@bug : this one completely fails :
	//typedef Modular<Integer>          FieldI;

	srand(time(NULL));
	bool fail = false ;

	// need to be prime for ftrsm
	FieldD FD(13);
	FieldD FD2(65563);
	FieldF FF(13);
	FieldF FF2(1069);
	FieldI FI(13);
	FieldI FI2(106739);
	int tot = 6;
#ifdef __LINBOX_HAVE_INT64
	//        typedef ModularBalanced<int64>  FieldU;
	typedef Modular<int64>          FieldU;
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
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FD);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FD);

		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FD2);

		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FD2);

		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FD2);
		ret+=test_ftrmm<FieldD,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FD2);

		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FF);

		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FF);

		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FF);

		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FF);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FF);

		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FF2);

		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FF2);

		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FF2);

		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FF2);
		ret+=test_ftrmm<FieldF,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FF2);


		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FI);

		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FI);

		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FI);

		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FI);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FI);

		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FI2);

		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FI2);

		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FI2);

		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FI2);
		ret+=test_ftrmm<FieldI,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FI2);

#ifdef __LINBOX_HAVE_INT64
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FU);

		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FU);

		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FU);

		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FU);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FU);

		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FU2);

		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasLeft,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FU2);

		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasUpper,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FU2);

		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasTrans,FFLAS::FflasNonUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasUnit>(FU2);
		ret+=test_ftrmm<FieldU,FFLAS::FflasRight,FFLAS::FflasLower,FFLAS::FflasNoTrans,FFLAS::FflasNonUnit>(FU2);
#endif
	}


#ifdef DEBUG
	std::cout << "# \033[1;33m>\033[0m ftr(s/m)m  passed " << ret << "/" << tot << "tests" <<std::endl;
#endif
	if (ret != tot) fail=true;
	int our = tot = 6*_LB_ITERS*2 ;
#ifdef __LINBOX_HAVE_INT64
	our = tot = tot+2*_LB_ITERS*2 ;
#endif
	//-------//
	/* APPLY */
	//-------//

	for (size_t r = 0 ; r < _LB_ITERS ; ++r)
	{
		our+= test_applyP(FD);
		our+= test_applyP(FD2);
		our+= test_applyP(FI);
		our+= test_applyP(FI2);
		our+= test_applyP(FF);
		our+= test_applyP(FF2);
#ifdef __LINBOX_HAVE_INT64
		our+= test_applyP(FU);
		our+= test_applyP(FU2);
#endif

	}
#ifdef DEBUG
	std::cout << "# \033[1;33m>\033[0m applyP  passed " << our << "/" << tot << "tests" <<std::endl;
#endif
	if(our != tot) fail = true;
	return (fail) ;

}

#undef _LB_TIME
#undef _LB_MAX_SZ
#undef _LB_ITERS
#undef _LB_DEBUG

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
