
/* fflas/fflas_ftrsm.inl
 * Copyright (C) 2005 Clement Pernet
 *
 * Written by Clement Pernet <Clement.Pernet@imag.fr>
 *
 * See COPYING for license information.
 */


//---------------------------------------------------------------------
// ftrsm: TRiangular System solve with matrix
// Computes  B <- alpha.op(A^-1).B,  B <- alpha.B.op(A^-1)
// B is M*N, A is M*M if Side==FflasLeft, N*N if Side==FflasRight
// Warning : unsafe with Trans ==  FflasTrans (debugging in progress)
//---------------------------------------------------------------------
template<class Field>
inline void
FFLAS::ftrsm (const Field& F, const FFLAS_SIDE Side,
	      const FFLAS_UPLO Uplo, 
	      const FFLAS_TRANSPOSE TransA,
	      const FFLAS_DIAG Diag, 
	      const size_t M, const size_t N,
	      const typename Field::Element alpha,
	      typename Field::Element * A, const size_t lda,
	      typename Field::Element * B, const size_t ldb) 
{
	if (!M || !N ) return; 

		
	if ( Side==FflasLeft ){
		if ( Uplo==FflasUpper){
			if (TransA == FflasNoTrans){
				if (Diag == FflasUnit)
					ftrsmLeftUpperNoTransUnit<typename Field::Element> ()(F,M,N,A,lda,B,ldb);
				else
					ftrsmLeftUpperNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
			} else {
				if (Diag == FflasUnit)
					ftrsmLeftUpperTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
				else
					ftrsmLeftUpperTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
			}
		} else {
			if (TransA == FflasNoTrans){
				if (Diag == FflasUnit)
					ftrsmLeftLowerNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
				else
					ftrsmLeftLowerNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
			} else {
				if (Diag == FflasUnit)
					ftrsmLeftLowerTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
				else
					ftrsmLeftLowerTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
			}
		}
	} else {
		if ( Uplo == FflasUpper){
			if (TransA == FflasNoTrans){
				if (Diag == FflasUnit)
					ftrsmRightUpperNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
				else
					ftrsmRightUpperNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
			} else {
				if (Diag == FflasUnit)
					ftrsmRightUpperTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
				else
					ftrsmRightUpperTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
			}
		} else {
			if (TransA == FflasNoTrans){
				if (Diag == FflasUnit)
					ftrsmRightLowerNoTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
				else
					ftrsmRightLowerNoTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
			} else {
				if (Diag == FflasUnit)
					ftrsmRightLowerTransUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
				else
					ftrsmRightLowerTransNonUnit<typename Field::Element>()(F,M,N,A,lda,B,ldb);
			}
		}
	}
	if (!F.isOne(alpha))
		for (size_t i=0; i< M; ++i)
			for (size_t j=0; j<N; ++j)
				F.mulin(*(B+i*ldb+j),alpha);
	
}

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT



#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__GENERIC
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__GENERIC
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT
//==

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT



#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__DOUBLE
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__DOUBLE
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__LEFT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__LEFT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT



#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__UP
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__UP
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT


#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__NOTRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__UNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__NONUNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NONUNIT

#define __FFLAS__FLOAT
#define __FFLAS__RIGHT
#define __FFLAS__LOW
#define __FFLAS__TRANSPOSE
#define __FFLAS__UNIT
#include "fflas_ftrsm_src.inl"
#undef __FFLAS__FLOAT
#undef __FFLAS__RIGHT
#undef __FFLAS__LOW
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__UNIT
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
