/* fflas/test-ftrmm.inl
 * Copyright (C) 2010 LinBox
 *
 * Written by Brice Boyer <Brice.Boyer@imag.fr>
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"
#include "linbox/fflas/fflas.h"

//!@todo !
//#define __LINBOX_HAVE_INT64

#include "linbox/field/modular-balanced-int32.h"
#ifdef __LINBOX_HAVE_INT64
#include "linbox/field/modular-balanced-int64.h"
#endif
#include "linbox/field/modular-balanced-double.h"
#include "linbox/field/modular-balanced-float.h"
#include "Matio.h"

//#define _LB_TIME
#define _LB_MAX_SZ 100
#define _LB_ITERS 3 
//#define _LB_DEBUG


using namespace LinBox ;

template<class Field, FFLAS::FFLAS_SIDE Side,FFLAS::FFLAS_UPLO UpLo, FFLAS::FFLAS_TRANSPOSE Trans, FFLAS::FFLAS_DIAG Diag >
int test_ftrmm(const Field & F) 
{
	size_t M    = random()%_LB_MAX_SZ ;
	size_t lda  = random()%_LB_MAX_SZ ; 
	size_t N    = random()%_LB_MAX_SZ ;
	size_t ldb  = random()%_LB_MAX_SZ ;

	if (lda<M) std::swap(lda,M);
	if (ldb<N) std::swap(ldb,N);

//        std::cout << "#M x N :"<< M << 'x' << N << std::endl; 
//        std::cout << "#lda x ldb :"<< lda << 'x' << ldb << std::endl; 

	typedef typename Field::Element Element;

	size_t Ldim = M;
	if (Side == FFLAS::FflasRight){
		Ldim = N ;
		while (lda<N)  lda  = random()%_LB_MAX_SZ ; 
	}

	size_t K = std::max(lda,ldb)+10 ; 
	K *= K ;// size of allocated matrices


	Element * A = new Element[K];
	Element * B = new Element[K];
	Element * C = new Element[M*N] ; // le résultat, le vrai, le bon.

	typedef typename Field::RandIter RandIter;
	RandIter G(F);
	NonzeroRandIter<Field> Gn(F,G); 

	/* init A,B,C and more...*/
	for (size_t i = 0 ; i < K ; ++i)
		G.random(*(A+i) ) ;
	for (size_t i = 0 ; i < K ; ++i)
		G.random(*(B+i) ) ;

	Element zero ; F.init(zero,0UL);
	for (size_t i = 0 ; i < M*N ; ++i)
		*(C+i) = zero;
	Element alpha ;
	Gn.random(alpha);
#ifdef _LB_DEBUG

		std::cout<<"#-------------T-------------" <<std::endl;
		std::cout << "T :=" ;
		write_field(F,std::cout,UpLo,Diag,A,Ldim,Ldim,lda,true,true);
		std::cout<<"#-------------M-------------" <<std::endl;
		std::cout << "M :=" ;
		write_field(F,std::cout,B,M,N,ldb,true,true);
#endif
	/*  compute C */
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
		F.mulin(*(A+i),alpha);

	/*  compute B */
	FFLAS::ftrmm(F,Side,UpLo,Trans,Diag,M,N,alpha,A,lda,B,ldb);

	int ret = 0 ;
	for (size_t i = 0 ; i < M && !ret ; ++i)
		for (size_t j = 0 ; j < N && !ret ; ++j)
			if (!F.areEqual(*(C+i*N+j),*(B+i*ldb+j))) 
				ret = -1  ;
#ifdef _LB_DEBUG
	if (ret)
	{

		std::cout<<"#------------a--------------" <<std::endl;
		std::cout << "alpha := " << alpha << ':' << std::endl;
		std::cout<<"#------------C--------------" <<std::endl;
		std::cout << "C := " ; write_field(F,std::cout,C,M,N,N,true,true);
		std::cout << "#------------B--------------"                     << std::endl;
		std::cout << "B := " ; write_field(F,std::cout,B,M,N,ldb,true,true);
		std::cout << "N := " ; 
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
	if (ret)
		std::cout << "#ftrmm fail" ;
	else
		std::cout << "#ftrmm success";

		std::cout << " with " << ((Side == FFLAS::FflasLeft)?("left"):("right")) ;
		std::cout << " with " ;
		std::cout << ((UpLo==FFLAS::FflasUpper)?("upper "):("lower ")) ;
		std::cout << ((Diag==FFLAS::FflasUnit)?(""):("non-")) ;
		std::cout << "unit " ;
		std::cout << ((Trans==FFLAS::FflasTrans)?(""):("non-")) ;
		std::cout << "transposed triangular matrices."  ;
		std::cout << "on " ;
		F.write(std::cout) ;
		std::cout << std::endl;
#endif



	return ret ;

}




int main() 
{
	typedef ModularBalanced<float>  FieldF;
	typedef ModularBalanced<double> FieldD;
	typedef ModularBalanced<int32>    FieldI;

	FieldD FD(13);
	FieldD FD2(65561);
	FieldF FF(13);
	FieldF FF2(1067);
	FieldI FI(13);
	FieldI FI2(1067);
	int tot = 16*6;
#ifdef __LINBOX_HAVE_INT64
	typedef ModularBalanced<int64>    FieldU;
	FieldU FU(13);
	FieldU FU2(13132131);
	tot += 16*2 ;
#endif
	int ret = tot;
	

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



	std::cout << "trmm  passed " << ret << "/" << tot << std::endl;

}

#undef _LB_TIME
#undef _LB_MAX_SZ
#undef _LB_ITERS
#undef _LB_DEBUG

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
