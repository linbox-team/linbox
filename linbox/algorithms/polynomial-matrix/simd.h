/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2013  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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
 */

#ifndef __LINBOX_SIMD__
#define __LINBOX_SIMD__


//#include <immintrin.h>
#include "fflas-ffpack/utils/fflas_intrinsic.h"
#include <iostream>


#ifdef __LINBOX_HAVE_AVX2_INSTRUCTIONS
/* 256 bits CODE HERE */

// define 256 bits simd vector type
typedef __m256i  _vect256_t; 

/* 256 bits CODE HERE */
// C=C+AB using X as temporary
#define VEC256_MADD_32(C,A,B,X)					\
	X= _mm256_mul_epi32(A,B); C= _mm256_add_epi64(C,X);

// C=A*B (2 op 32x32->64)
#define VEC256_MUL_32(C,A,B)			\
	C= _mm256_mul_epu32(A,B); 

// C=A*B (4 op 32x32->32 low product)
#define VEC256_MUL_LO_32(C,A,B)			\
	C=  _mm256_mullo_epi32(A,B);

// C=A+B
#define VEC256_ADD_32(C,A,B)			\
	C= _mm256_add_epi32(A,B);
// C=A+B
#define VEC256_ADD_64(C,A,B)			\
	C= _mm256_add_epi64(A,B);

// A+=B
#define VEC256_ADDIN_64(A,B)			\
	A= _mm256_add_epi64(A,B);

// C=A-B
#define VEC256_SUB_32(C,A,B)			\
	C= _mm256_sub_epi32(A,B);

// C = A mod P using T as temporary  (A must lie in [0 2P[ )
#define VEC256_MOD_P(C,A,P,T)						\
	T=_mm256_cmpgt_epi32(P,A);					\
	C=_mm256_sub_epi32(A,_mm256_andnot_si256(T,P));			
// Rk: __m128i _mm_cmpgt_epi64 (__m128d a, __m128d b) // compare a>b si vrai renvoi 0xFFFFFFFFFFFFFFFF

// C = A >> X
#define VEC256_RSHIFT_64(C,A,X)			\
	C=_mm256_srli_epi64 (A,X);

// C = A << X
#define VEC256_LSHIFT_64(C,A,X)			\
	C=_mm256_slli_epi64 (A,X);

// C = A || X
#define VEC256_OR(C,A,X)				\
	C=_mm256_or_si256(A,X);

// C = unpack_lo32(A,B)
#define VEC256_UNPACK_LO_32(C,A,B)			\
	 C = _mm256_unpacklo_epi32(A,B);

// C = unpack_hi32(A,B)
#define VEC256_UNPACK_HI_32(C,A,B)			\
	 C = _mm256_unpackhi_epi32(A,B);

// C = unpack_lo32(A,B)
#define VEC256_UNPACK_LO_64(C,A,B)			\
	 C = _mm256_unpacklo_epi64(A,B);

// C = unpack_hi32(A,B)
#define VEC256_UNPACK_HI_64(C,A,B)			\
	 C = _mm256_unpackhi_epi64(A,B);

// C = unpack_lo128(A,B)
#define VEC256_UNPACK_LO_128(C,A,B)			\
	C = _mm256_permute2x128_si256(A,B,32);

// C = unpack_hi128(A,B)
#define VEC256_UNPACK_HI_128(C,A,B)			\
	C = _mm256_permute2x128_si256(A,B,49);
 
// C= blend_epi32(A,B,X)			
#define VEC256_BLEND_32(C,A,B,X)		\
	C = _mm256_blend_epi32(A,B,X);

// C =shuffle_32(A,X)				
#define VEC256_SHUFFLE_32(C,A,X)		\
	C = _mm256_shuffle_epi32(A,X);

// C=X[0,1]
#define VEC256_LOAD(C,X)							\
	C=_mm256_load_si256(reinterpret_cast<const __m256i*>(X));
#define VEC256_LOADU(C,X)							\
	C=_mm256_loadu_si256(reinterpret_cast<const __m256i*>(X));
// C[0,1]=X
#define VEC256_STORE(C,X)						\
	_mm256_store_si256(reinterpret_cast<__m256i*>(C),X);
#define VEC256_STOREU(C,X)						\
	_mm256_storeu_si256(reinterpret_cast<__m256i*>(C),X);
#define VEC256_SET_64(C,X)						\
	C=_mm256_set1_epi64x(X);
#define VEC256_SET_32(C,X)						\
	C=_mm256_set1_epi32(X);

#define VEC256_PRINT_32(X)			\
	{uint32_t T[8];	VEC256_STORE(T,X);std::cout<<"[ "<<T[0]<<","<<T[1]<<","<<T[2]<<","<<T[3]<<","<<T[4]<<","<<T[5]<<","<<T[6]<<","<<T[7]<<"]";}

#endif

/* 128 bits CODE HERE */

// define 128 bits simd vector type
typedef __m128i  _vect128_t;


// C=C+AB using X as temporary
#define VEC128_MADD_32(C,A,B,X)				\
	X= _mm_mul_epi32(A,B); C= _mm_add_epi64(C,X);

// C=A*B (2 op 32x32->64)
#define VEC128_MUL_32(C,A,B)			\
	C= _mm_mul_epu32(A,B);

// C=A*B (4 op 32x32->32 low product)
#define VEC128_MUL_LO_32(C,A,B)			\
	C=  _mm_mullo_epi32(A,B);

// C=A+B
#define VEC128_ADD_32(C,A,B)			\
	C= _mm_add_epi32(A,B);

// C=A+B
#define VEC128_ADD_64(C,A,B)			\
	C= _mm_add_epi64(A,B);

// A+=B
#define VEC128_ADDIN_64(A,B)			\
	A= _mm_add_epi64(A,B);

// C=A-B
#define VEC128_SUB_32(C,A,B)			\
	C= _mm_sub_epi32(A,B);

// C = A mod P using T as temporary  (A must lie in [0 2P[ )
#define VEC128_MOD_P(C,A,P,T)						\
	T=_mm_cmplt_epi32(A,P);						\
	C=_mm_sub_epi32(A,_mm_andnot_si128(T,P));			
// Rk: __m128i _mm_cmpgt_epi64 (__m128d a, __m128d b) // compare a>b si vrai renvoi 0xFFFFFFFFFFFFFFFF

// C = A >> X
#define VEC128_RSHIFT_64(C,A,X)			\
	C=_mm_srli_epi64 (A,X);

// C = A << X
#define VEC128_LSHIFT_64(C,A,X)			\
	C=_mm_slli_epi64 (A,X);

// C = A || X
#define VEC128_OR(C,A,X)				\
	C=_mm_or_si128(A,X);

// C = unpack_lo32(A,B)
#define VEC128_UNPACK_LO_32(C,A,B)			\
	 C = _mm_unpacklo_epi32(A,B);

// C = unpack_hi32(A,B)
#define VEC128_UNPACK_HI_32(C,A,B)			\
	 C = _mm_unpackhi_epi32(A,B);

// C = unpack_lo64(A,B)
#define VEC128_UNPACK_LO_64(C,A,B)			\
	 C = _mm_unpacklo_epi64(A,B);

// C = unpack_hi64(A,B)
#define VEC128_UNPACK_HI_64(C,A,B)			\
	 C = _mm_unpackhi_epi64(A,B);

// C =shuffle_32(A,X)				
#define VEC128_SHUFFLE_32(C,A,X)			\
	C = _mm_shuffle_epi32(A,X);

// C=X[0,1]
#define VEC128_LOAD(C,X)							\
	C=_mm_load_si128(reinterpret_cast<const __m128i*>(X));
#define VEC128_LOADU(C,X)							\
	C=_mm_loadu_si128(reinterpret_cast<const __m128i*>(X));
// C[0,1]=X
#define VEC128_STORE(C,X)						\
	_mm_store_si128(reinterpret_cast<__m128i*>(C),X);
#define VEC128_STOREU(C,X)						\
	_mm_storeu_si128(reinterpret_cast<__m128i*>(C),X);
#define VEC128_SET_64(C,X)				\
	C=_mm_set1_epi64x(X);
#define VEC128_SET_32(C,X)				\
	C=_mm_set1_epi32(X);

#define VEC128_PRINT_32(X)				\
	{ uint32_t T[4];VEC128_STORE(T,X); cout<<"[ "<<T[0]<<","<<T[1]<<","<<T[2]<<","<<T[3]<<"]";}

// END OF 128 BITS CODE

 
// C=A*B (4 op 32x32->32 high product) // A and with a mask can be used to remove the last two shift with C
#define VEC128_MUL_HI_32(C,A,B,A1,B1)			\
	VEC128_MUL_32(C,A,B);				\
	VEC128_RSHIFT_64(A1,A,32);			\
	VEC128_RSHIFT_64(B1,B,32);			\
	VEC128_MUL_32(A1,A1,B1);			\
	VEC128_RSHIFT_64(C,C,32);			\
	VEC128_RSHIFT_64(A1,A1,32);			\
	VEC128_LSHIFT_64(A1,A1,32);			\
	VEC128_OR(C,C,A1);			       
#define VEC256_MUL_HI_32(C,A,B,A1,B1)			\
	VEC256_MUL_32(C,A,B);				\
	VEC256_RSHIFT_64(A1,A,32);			\
	VEC256_RSHIFT_64(B1,B,32);			\
	VEC256_MUL_32(A1,A1,B1);			\
	VEC256_RSHIFT_64(C,C,32);			\
	VEC256_RSHIFT_64(A1,A1,32);			\
	VEC256_LSHIFT_64(A1,A1,32);			\
	VEC256_OR(C,C,A1);			
	
// C= A+B mod P using T as temporary
#define VEC128_ADD_MOD(C,A,B,P,T)			\
	VEC128_ADD_32(C,A,B); VEC128_MOD_P(C,C,P,T); 
#define VEC256_ADD_MOD(C,A,B,P,T)			\
	VEC256_ADD_32(C,A,B); VEC256_MOD_P(C,C,P,T); 

// C= A*X mod P using T1, T2 as temporaries
#define VEC128_MUL_MOD(C,A,X,P,Xp,Q,T1,T2)		\
	VEC128_MUL_HI_32(Q,A,Xp,T1,T2);			\
	VEC128_MUL_LO_32(C,A,X);			\
	VEC128_MUL_LO_32(T1,Q,P);			\
	VEC128_SUB_32(C,C,T1);
#define VEC256_MUL_MOD(C,A,X,P,Xp,Q,T1,T2)		\
	VEC256_MUL_HI_32(Q,A,Xp,T1,T2);			\
	VEC256_MUL_LO_32(C,A,X);			\
	VEC256_MUL_LO_32(T1,Q,P);			\
	VEC256_SUB_32(C,C,T1);



#endif // end of file
