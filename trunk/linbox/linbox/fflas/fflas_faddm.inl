/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* fflas/fflas_faddm.inl
 * Copyright (C) 2010 LinBox
 *
 * Written by Brice Boyer <Brice.Boyer@imag.fr>
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_fflas_faddm_H
#define __LINBOX_fflas_faddm_H

#ifdef __LINBOX_HAVE_SSE2
#include <emmintrin.h>
#endif

/** faddm
 * A <- A+op(B)
 * with op(B) = B or B^T
 */

template<class Field>
inline void FFLAS::faddm(const Field & F,
		  const FFLAS_TRANSPOSE transA,
		  const size_t M, const size_t N,
		  const typename Field::Element * A, const size_t lda,
		        typename Field::Element * B, const size_t ldb)
{
	if (!M || !N) return ;
	if (transA ==  FflasNoTrans)
		faddmNoTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb);
	else
		faddmTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb);
	return ;
}

/** faddm
 * C <- op(A)+op(B)
 * with op(B) = B or B^T
 */

template<class Field>
inline void FFLAS::faddm(const Field & F,
		  const FFLAS_TRANSPOSE transA,
		  const FFLAS_TRANSPOSE transB,
		  const size_t M, const size_t N,
		  const typename Field::Element * A, const size_t lda,
		  const typename Field::Element * B, const size_t ldb,
		        typename Field::Element * C, const size_t ldc )
{
	if (!M || !N) return ;
	if (transA ==  FflasNoTrans)
		if (transB ==  FflasNoTrans)
			faddmNoTransNoTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb,C,ldc);
		else
			faddmNoTransTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb,C,ldc);
	else
		if (transB ==  FflasNoTrans)
			faddmTransNoTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb,C,ldc);
		else
			faddmTransTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb,C,ldc);

	return ;
}

/** fsubm
 * A <- A-op(B)
 * with op(B) = B or B^T
 */

template<class Field>
inline void FFLAS::fsubm(const Field & F,
		  const FFLAS_TRANSPOSE transA,
		  const size_t M, const size_t N,
		  const typename Field::Element * A, const size_t lda,
		        typename Field::Element * B, const size_t ldb)
		  
{
	if (!M || !N) return ;
	if (transA ==  FflasNoTrans)
		fsubmNoTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb);
	else
		fsubmTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb);
	return ;
}

/** fsubm
 * C <- op(A)-op(B)
 * with op(B) = B or B^T
 */

template<class Field>
inline void FFLAS::fsubm(const Field & F,
		  const FFLAS_TRANSPOSE transA,
		  const FFLAS_TRANSPOSE transB,
		  const size_t M, const size_t N,
		  const typename Field::Element * A, const size_t lda,
		  const typename Field::Element * B, const size_t ldb,
		        typename Field::Element * C, const size_t ldc )
{
	if (!M || !N) return ;
	if (transA ==  FflasNoTrans)
		if (transB ==  FflasNoTrans)
			fsubmNoTransNoTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb,C,ldc);
		else
			fsubmNoTransTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb,C,ldc);
	else
		if (transB ==  FflasNoTrans)
			fsubmTransNoTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb,C,ldc);
		else
			fsubmTransTrans<typename Field::Element>()(F,M,N,A,lda,B,ldb,C,ldc);

	return ;
}

#undef __FFLAS__FLOAT
#undef __FFLAS__DOUBLE
#undef __FFLAS__GENERIC
#undef __FFLAS__NOTRANSPOSE
#undef __FFLAS__ATRANSPOSE
#undef __FFLAS__ANOTRANSPOSE
#undef __FFLAS__BTRANSPOSE
#undef __FFLAS__BNOTRANSPOSE
#undef __FFLAS__TRANSPOSE
#undef __FFLAS__NOTRANSPOSE

#define __FFLAS__GENERIC
#define __FFLAS__NOTRANSPOSE // no transpose
#include "fflas_faddmin_src.inl"
#undef  __FFLAS__NOTRANSPOSE
#undef  __FFLAS__GENERIC

#define __FFLAS__GENERIC
#define __FFLAS__TRANSPOSE // no transpose
#include "fflas_faddmin_src.inl"
#undef  __FFLAS__TRANSPOSE
#undef  __FFLAS__GENERIC

#define __FFLAS__FLOAT
#define __FFLAS__NOTRANSPOSE // no transpose
#include "fflas_faddmin_src.inl"
#undef  __FFLAS__NOTRANSPOSE
#undef  __FFLAS__FLOAT

#define __FFLAS__FLOAT
#define __FFLAS__TRANSPOSE // no transpose
#include "fflas_faddmin_src.inl"
#undef  __FFLAS__TRANSPOSE
#undef  __FFLAS__FLOAT

#define __FFLAS__DOUBLE
#define __FFLAS__NOTRANSPOSE // no transpose
#include "fflas_faddmin_src.inl"
#undef  __FFLAS__NOTRANSPOSE
#undef  __FFLAS__DOUBLE

#define __FFLAS__DOUBLE
#define __FFLAS__TRANSPOSE // no transpose
#include "fflas_faddmin_src.inl"
#undef  __FFLAS__TRANSPOSE
#undef  __FFLAS__DOUBLE

//
#define __FFLAS__GENERIC
#define __FFLAS__ANOTRANSPOSE // no A transpose
#define __FFLAS__BNOTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BNOTRANSPOSE // no B transpose
#undef  __FFLAS__ANOTRANSPOSE // no A transpose
#undef  __FFLAS__GENERIC

#define __FFLAS__GENERIC
#define __FFLAS__ANOTRANSPOSE // no A transpose
#define __FFLAS__BTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BTRANSPOSE // no B transpose
#undef  __FFLAS__ANOTRANSPOSE // no A transpose
#undef  __FFLAS__GENERIC


#define __FFLAS__GENERIC
#define __FFLAS__ATRANSPOSE // no A transpose
#define __FFLAS__BNOTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BNOTRANSPOSE // no B transpose
#undef  __FFLAS__ATRANSPOSE // no A transpose
#undef  __FFLAS__GENERIC

#define __FFLAS__GENERIC
#define __FFLAS__ATRANSPOSE // no A transpose
#define __FFLAS__BTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BTRANSPOSE // no B transpose
#undef  __FFLAS__ATRANSPOSE // no A transpose
#undef  __FFLAS__GENERIC

//
#define __FFLAS__FLOAT
#define __FFLAS__ANOTRANSPOSE // no A transpose
#define __FFLAS__BNOTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BNOTRANSPOSE // no B transpose
#undef  __FFLAS__ANOTRANSPOSE // no A transpose
#undef  __FFLAS__FLOAT

#define __FFLAS__FLOAT
#define __FFLAS__ANOTRANSPOSE // no A transpose
#define __FFLAS__BTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BTRANSPOSE // no B transpose
#undef  __FFLAS__ANOTRANSPOSE // no A transpose
#undef  __FFLAS__FLOAT


#define __FFLAS__FLOAT
#define __FFLAS__ATRANSPOSE // no A transpose
#define __FFLAS__BNOTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BNOTRANSPOSE // no B transpose
#undef  __FFLAS__ATRANSPOSE // no A transpose
#undef  __FFLAS__FLOAT

#define __FFLAS__FLOAT
#define __FFLAS__ATRANSPOSE // no A transpose
#define __FFLAS__BTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BTRANSPOSE // no B transpose
#undef  __FFLAS__ATRANSPOSE // no A transpose
#undef  __FFLAS__FLOAT

//
#define __FFLAS__DOUBLE
#define __FFLAS__ANOTRANSPOSE // no A transpose
#define __FFLAS__BNOTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BNOTRANSPOSE // no B transpose
#undef  __FFLAS__ANOTRANSPOSE // no A transpose
#undef  __FFLAS__DOUBLE

#define __FFLAS__DOUBLE
#define __FFLAS__ANOTRANSPOSE // no A transpose
#define __FFLAS__BTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BTRANSPOSE // no B transpose
#undef  __FFLAS__ANOTRANSPOSE // no A transpose
#undef  __FFLAS__DOUBLE


#define __FFLAS__DOUBLE
#define __FFLAS__ATRANSPOSE // no A transpose
#define __FFLAS__BNOTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BNOTRANSPOSE // no B transpose
#undef  __FFLAS__ATRANSPOSE // no A transpose
#undef  __FFLAS__DOUBLE

#define __FFLAS__DOUBLE
#define __FFLAS__ATRANSPOSE // no A transpose
#define __FFLAS__BTRANSPOSE // no B transpose
#include "fflas_faddm_src.inl"
#undef  __FFLAS__BTRANSPOSE // no B transpose
#undef  __FFLAS__ATRANSPOSE // no A transpose
#undef  __FFLAS__DOUBLE


#endif // __LINBOX_fflas_faddm_H

