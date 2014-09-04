/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/*
 * Copyright (C) 2014  Pascal Giorgi, Romain Lebreton
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
 *            Romain Lebreton <romain.lebreton@lirmm.fr>
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


#ifndef __LINBOX_polynomial_fft_transform_simd_INL
#define __LINBOX_polynomial_fft_transform_simd_INL

#include "linbox/algorithms/polynomial-matrix/simd.h"

namespace LinBox{

	/******************************************************************************************************************
	 ******************************************************************************************************************
	 ***********************************    FFT with SSE CODE   *******************************************************
	 ******************************************************************************************************************
	 ******************************************************************************************************************/
  
	template <class Field>
	inline void FFT_transform<Field>::reduce128_modp(Element* ABCD, const _vect128_t& P) {
		_vect128_t V1,T;
		// V1=[A B C D], V2=[E F G H]
		VEC128_LOAD(V1,ABCD);
		VEC128_MOD_P(V1,V1,P,T);
		VEC128_STORE(ABCD,V1);
	}

	/*-----------------------------------*/
	/*--        Butterflies DIF      ----*/
	/*-----------------------------------*/


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_4x1_SSE(Element* ABCD, Element* EFGH,
								      const Element* alpha,
								      const Element* alphap,
								      const _vect128_t& P, const _vect128_t& P2) 
	{
		_vect128_t V1,V2,V3,V4,W,Wp,T;
		// V1=[A B C D], V2=[E F G H]
		VEC128_LOAD(V1,ABCD);
		VEC128_LOAD(V2,EFGH);
		VEC128_LOAD(W ,alpha);
		VEC128_LOAD(Wp,alphap);	
		// V3 = V1 + V2 mod 
		VEC128_ADD_MOD(V3,V1,V2,P2,T);
		VEC128_STORE(ABCD,V3);	
		// V4 = (V1+(2P-V2))alpha mod 2P
		VEC128_SUB_32(T,V2,P2);
		VEC128_SUB_32(V4,V1,T);
		VEC128_MUL_MOD(T,V4,W,P,Wp,V1,V2,V3);// V3 is the result
		VEC128_STORE(EFGH,T);
	}


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_4x1_SSE_laststep(Element* ABCD, Element* EFGH, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,V5;
		// V1=[A B C D], V2=[E F G H]
		VEC128_LOAD(V1,ABCD);
		VEC128_LOAD(V2,EFGH);
		// V3 = [A C B D], V4 = [E G F H]
		VEC128_SHUFFLE_32(V3,V1,0xD8);
		VEC128_SHUFFLE_32(V4,V2,0xD8);    
		// V1 = [A E C G], V2 = [B F D H]
		VEC128_UNPACK_LO_32(V1,V3,V4);
		VEC128_UNPACK_HI_32(V2,V3,V4);    
		// V3 = V1 + V2 mod 2P
		VEC128_ADD_MOD(V3,V1,V2,P2,V5);	
		// V4 = V1 + (2P - V2) mod 2P
		VEC128_SUB_32(V5,V2,P2);
		VEC128_SUB_32(V2,V1,V5);
		VEC128_MOD_P(V4,V2,P2,V5);
		// V1 = [A C E G], V2 = [B D F H]
		VEC128_SHUFFLE_32(V1,V3,0xD8);
		VEC128_SHUFFLE_32(V2,V4,0xD8);
		// V3 = [A B C D], V4 = [E F G H]
		VEC128_UNPACK_LO_32(V3,V1,V2);
		VEC128_UNPACK_HI_32(V4,V1,V2);    
		// Store
		VEC128_STORE(ABCD,V3);
		VEC128_STORE(EFGH,V4);
	}

	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_4x2_SSE(Element* ABCD, Element* EFGH, Element* IJKL, Element* MNOP,
								      const Element* alpha,  const Element*beta ,  const Element* gamma,
								      const Element* alphap, const Element*betap , const Element* gammap,
								      const _vect128_t& P, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,W,Wp,T1,T2,T3,T4,T5,T6,T7,T8;
	
		// V1=[A B C D], V2=[E F G H], V3=[I J K L], V4=[M N O P]
		VEC128_LOAD(V1,ABCD);
		VEC128_LOAD(V2,IJKL);
		VEC128_LOAD(W ,alpha);
		VEC128_LOAD(Wp,alphap);
		/**************/
		// T1 = V1 + V2 mod 2P
		VEC128_ADD_MOD(T1,V1,V2,P2,T8);
		// T2 = (V1+(2P-V2))alpha mod 2P
		VEC128_SUB_32(T7,V2,P2);
		VEC128_SUB_32(T6,V1,T7);
		VEC128_MUL_MOD(T2,T6,W,P,Wp,T3,T4,T5);
		/**************/
		VEC128_LOAD(V3,EFGH);
		VEC128_LOAD(V4,MNOP);
		VEC128_LOAD(W ,beta);
		VEC128_LOAD(Wp,betap);
		/**************/
		// T3 = V3 + V4 mod 2P
		VEC128_ADD_MOD(T3,V3,V4,P2,T8);
		// T4 = (V3+(2P-V4))beta mod 2P
		VEC128_SUB_32(T7,V4,P2);
		VEC128_SUB_32(T6,V3,T7);
		VEC128_MUL_MOD(T4,T6,W,P,Wp,V1,V2,T8);// T1 is the result
		/**************/
		VEC128_LOAD(W ,gamma);
		VEC128_LOAD(Wp,gammap);
		/**************/
		// V1 = T1 + T3 mod 2P
		VEC128_ADD_MOD(V1,T1,T3,P2,T8);
		// V3 = (T1+(2P-T3))gamma mod 2P
		VEC128_SUB_32(T7,T3,P2);
		VEC128_SUB_32(T6,T1,T7);
		VEC128_MUL_MOD(V3,T6,W,P,Wp,T3,T5,T8);// T1 is the result
		/**************/
		// V2 = T2 + T4 mod 2P
		VEC128_ADD_MOD(V2,T2,T4,P2,T8);
		// V4 = (T2+(2P-T4))gamma mod 2P
		VEC128_SUB_32(T7,T4,P2);
		VEC128_SUB_32(T6,T2,T7);
		VEC128_MUL_MOD(V4,T6,W,P,Wp,T1,T3,T8);// T1 is the result
		/**************/
		VEC128_STORE(ABCD,V1);
		VEC128_STORE(EFGH,V3);
		VEC128_STORE(IJKL,V2);
		VEC128_STORE(MNOP,V4);
	}



	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_4x2_SSE_last2step(Element* ABCD, Element* EFGH,
										const _vect128_t& W,
										const _vect128_t& Wp,
										const _vect128_t& P, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,V5,V6,V7;
		// V1=[A B C D], V2=[E F G H]
		VEC128_LOAD(V1,ABCD);
		VEC128_LOAD(V2,EFGH);
		// V3=[A E B F], V4=[C G D H]
		VEC128_UNPACK_LO_32(V3,V1,V2);
		VEC128_UNPACK_HI_32(V4,V1,V2);
		// V1 = V3 + V4 mod 2P
		// P2 = [2p 2p 2p 2p]
		VEC128_ADD_MOD(V1,V3,V4,P2,V5);
		// V2 = (V3+(2P-V4))alpha mod 2P
		VEC128_SUB_32(V5,V4,P2);
		VEC128_SUB_32(V6,V3,V5);
		VEC128_MOD_P(V2,V6,P2,V2);
		// V4 = [D D H H]
		VEC128_UNPACK_HI_32(V4,V2,V2);		
		// V6 = V4 * Wp mod 2^64
		// Wp = [Wp ? Wp ?]
		VEC128_MUL_32(V7,V4,Wp);
		VEC128_MUL_LO_32(V5,V7,P);
		// At this point V4= [? Q_D*p ? Q_H*p]
		// V5 = [D D H H] * [W W W W] mod 2^32
		VEC128_MUL_LO_32(V6,V4,W);
		VEC128_SUB_32(V4,V6,V5);
		VEC128_SHUFFLE_32(V3,V4,0xDD);
		//At this point, V2 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]
		// At this time I have V1=[A E B F], V2=[C G ? ?], V3=[? ? D H]
		// I need V3 = [A C E G], V4 = [B D F H]
		VEC128_UNPACK_HI_32(V4,V1,V3);
		VEC128_UNPACK_LO_32(V3,V1,V2);		
		// V1 = V3 + V4 mod 2P
		VEC128_ADD_MOD(V1,V3,V4,P2,V5);
		// V2 = V3 + (2P - V4) mod 2P
		VEC128_SUB_32(V5,V4,P2);
		VEC128_SUB_32(V6,V3,V5);
		VEC128_MOD_P(V2,V6,P2,V2);
		// Result in V1 = [A C E G]  and V2 = [B D F H]
		// Transform to V3=[A B C D], V4=[E F G H]
		VEC128_UNPACK_LO_32(V3,V1,V2);
		VEC128_UNPACK_HI_32(V4,V1,V2);
		// Store
		VEC128_STORE(ABCD,V3);
		VEC128_STORE(EFGH,V4);
	}


	/*-----------------------------------*/
	/*--        Butterflies DIT      ----*/
	/*-----------------------------------*/


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p_4x1_SSE(Element* ABCD, Element* EFGH,
								      const Element* alpha,
								      const Element* alphap,
								      const _vect128_t& P, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,W,Wp,T1,T2;
		// V1=[A B C D], V2=[E F G H]
		VEC128_LOAD(V1,ABCD);
		VEC128_LOAD(V2,EFGH);
		VEC128_LOAD(W ,alpha);
		VEC128_LOAD(Wp,alphap);
		// V3 = V1 mod 2P
		VEC128_MOD_P (V3,V1,P2,T1);	
		// V4 = V2 * W mod P
		VEC128_MUL_MOD(V4,V2,W,P,Wp,V1,T1,T2);	
		// V1 = V3 + V4
		VEC128_ADD_32(V1,V3,V4);
		VEC128_STORE(ABCD,V1);	
		// V2 = V3 - (V4 - 2P)
		VEC128_SUB_32(T1,V4,P2);
		VEC128_SUB_32(V2,V3,T1);
		VEC128_STORE(EFGH,V2);
	}

	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p_4x2_SSE_first2step(Element* ABCD, Element* EFGH,
										 const _vect128_t& W,
										 const _vect128_t& Wp,
										 const _vect128_t& P, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,T1,T2,T3,T4;
		// V1=[A B C D], V2=[E F G H]
		VEC128_LOAD(V1,ABCD);
		VEC128_LOAD(V2,EFGH);	
		// T1 = [A C B D], T2 = [E G F H]
		VEC128_SHUFFLE_32(T1,V1,0xD8);
		VEC128_SHUFFLE_32(T2,V2,0xD8);	
		// V1 = [A E C G], V2 = [B F D H]
		VEC128_UNPACK_LO_32(V1,T1,T2);
		VEC128_UNPACK_HI_32(V2,T1,T2);	
		// V3 = V1 + V2
		// Rk: No need for (. mod 2P) since entries are <P
		VEC128_ADD_32(V3,V1,V2);	
		// V4 = V1 + (P - V2)
		// Rk: No need for (. mod 2P) since entries are <P
		VEC128_SUB_32(T1,V2,P);
		VEC128_SUB_32(V4,V1,T1);
		// T1 = [D D H H]
		VEC128_UNPACK_HI_32(T1,V4,V4);	
		// T2 = T1 * Wp mod 2^64
		// Wp = [Wp ? Wp ?]
		VEC128_MUL_32(T2,T1,Wp);
		VEC128_MUL_LO_32(T3,T2,P);
		// At this point T3= [? Q_D*p ? Q_H*p]	
		// T4 = [D D H H] * [W W W W] mod 2^32
		VEC128_MUL_LO_32(T4,T1,W);
		VEC128_SUB_32(T1,T4,T3);
		VEC128_SHUFFLE_32(T2,T1,0XDD);
		//At this point, T2 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]	
		// At this time I have V3=[A E C G], V4=[B F ? ?], T2=[? ? D H]
		// I need V1 = [A B E F], V2 = [C D G H]
		VEC128_UNPACK_LO_32(V1,V3,V4);
		VEC128_UNPACK_HI_32(V2,V3,T2);
		// T1 = V1 + V2
		VEC128_ADD_32(T1,V1,V2);	
		// T2 = V1 - (V2 - 2P)
		VEC128_SUB_32(T3,V2,P2);
		VEC128_SUB_32(T2,V1,T3);	
		// Result in T1 = [A B E F]  and T2 = [C D G H]
		// Transform to V1=[A C B D], V2=[E G F H]
		VEC128_UNPACK_LO_32(V1,T1,T2);
		VEC128_UNPACK_HI_32(V2,T1,T2);	
		// Then T1=[A B C D], T2=[E F G H]
		VEC128_SHUFFLE_32(T1,V1,0xD8);
		VEC128_SHUFFLE_32(T2,V2,0xD8);
		// Store
		VEC128_STORE(ABCD,T1);
		VEC128_STORE(EFGH,T2);	
	}

	/*-----------------------------------*/
	/*--       SSE FFT functions     ----*/
	/*-----------------------------------*/

	template <class Field>
	template <class Polynomial>
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative4x1_SSE (Polynomial &fft) {
		_vect128_t P,P2;
		P  = _mm_set1_epi32(_pl);
		P2 = _mm_set1_epi32(_dpl);
		Element * tab_w = &pow_w [0];
		Element * tab_wp= &pow_wp[0];
		size_t w, f;
		for (w = n >> 1, f = 1; w >= 4; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1){
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < w; j+=4)
					//for (ptrdiff_t j = w-4; j >=0; j-=4)

#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
					Butterfly_DIF_mod2p_4x1_SSE(A0,A4, tab_w+j,tab_wp+j,P,P2);
#undef A0
#undef A4
			//std::cout<<fft<<std::endl;
		}
		// Last two steps
		if (n >= 8) {
			//    if (n >= 1024) {
			_vect128_t W,Wp;
			W = _mm_set1_epi32 ((int)tab_w [1]);
			Wp= _mm_set1_epi32 ((int)tab_wp[1]);

			for (size_t i = 0; i < f; i+=2)
#define A0 &fft[0] +  (i << 2)
#define A4 &fft[0] + ((i << 2)+4)
				Butterfly_DIF_mod2p_4x2_SSE_last2step(A0,A4,W,Wp,P,P2);
			//std::cout<<fft<<std::endl;
#undef A0
#undef A4
		} else {
			for (; w >= 1; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1)
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
		}
	}

	template <class Field>
	template <class Polynomial>
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative4x2_SSE (Polynomial &fft) {
		size_t w, f;
		_vect128_t P,P2;
		P  = _mm_set1_epi32(_pl);
		P2 = _mm_set1_epi32(_dpl);
		Element * tab_w =  &pow_w[0];
		Element * tab_wp= &pow_wp[0];
		for (w = n >> 1, f = 1; w >= 8; tab_w+=w+(w>>1), tab_wp+=w+(w>>1), w >>= 2, f <<= 2)
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < (w >> 1); j+=4) {
#define A0 &fft[0] + (i << 1)    *w+ j          
#define A1 &fft[0] + (i << 1)    *w+(j+(w >> 1))
#define A2 &fft[0] + ((i << 1)+1)*w+ j          
#define A3 &fft[0] + ((i << 1)+1)*w+(j+(w >> 1))

					Butterfly_DIF_mod2p_4x2_SSE(A0, A1, A2, A3,
								    tab_w +j, tab_w +j+(w >> 1), tab_w +j+w,
								    tab_wp+j, tab_wp+j+(w >> 1), tab_wp+j+w,
								    P,P2);
#undef A0
#undef A1
#undef A2
#undef A3
				}

		// Last two steps
		if (n >= 8) {
			if (w == 4) {
				for (size_t i = 0; i < f; i++)
#define A0 &fft[0] +  (i << 1)   *w
#define A4 &fft[0] + ((i << 1)+1)*w
					Butterfly_DIF_mod2p_4x1_SSE(A0,A4, tab_w,tab_wp,P,P2);
#undef A0
#undef A4
				tab_w+=w;
				tab_wp+=w;
				w >>= 1;
				f <<= 1;
			}

			_vect128_t W,Wp;
			W = _mm_set1_epi32 ((int)tab_w [1]);
			Wp= _mm_set1_epi32 ((int)tab_wp[1]);

			for (size_t i = 0; i < f; i+=2)
#define A0 &fft[0] +  (i << 2)
#define A4 &fft[0] + ((i << 2)+4)
				Butterfly_DIF_mod2p_4x2_SSE_last2step(A0,A4,W,Wp,P,P2);
#undef A0
#undef A4
		} else {
			for (; w >= 1; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1)
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
		}
	}

	template <class Field>
	template <class Polynomial>
	void FFT_transform<Field>::FFT_DIT_Harvey_mod4p_iterative4x1_SSE (Polynomial &fft) 
	{
		_vect128_t P,P2;
		VEC128_SET_32(P,_pl);
		VEC128_SET_32(P2,_dpl);
		// Last two steps
		if (n >= 8) {
			_vect128_t W,Wp;
			W = _mm_set1_epi32 ((int)pow_w [n-3]);
			Wp= _mm_set1_epi32 ((int)pow_wp[n-3]);
	    
			for (size_t i = 0; i < n; i+=8)
				Butterfly_DIT_mod4p_4x2_SSE_first2step(&fft[i],&fft[i+4],W,Wp,P,P2);

			Element * tab_w = &pow_w [n-8];
			Element * tab_wp= &pow_wp[n-8];
			for (size_t w = 4, f = n >> 3; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w){
				// w : witdh of butterflies
				// f : # families of butterflies
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j+=4)
#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
						Butterfly_DIT_mod4p_4x1_SSE(A0,A4, tab_w+j,tab_wp+j,P,P2);
		    
#undef A0
#undef A4
		    
			}
		} else {
			Element * tab_w = &pow_w [n-2];
			Element * tab_wp= &pow_wp[n-2];
			for (size_t w = 1, f = n >> 1; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w)
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
		}
	}


	/******************************************************************************************************************
	 ******************************************************************************************************************
	 ***********************************   FFT with AVX 2 CODE  *******************************************************
	 ******************************************************************************************************************
	 ******************************************************************************************************************/


#ifdef __LINBOX_HAVE_AVX2

	template <class Field>
	inline void FFT_transform<Field>::reduce256_modp(Element* ABCD, const _vect256_t& P) {
		_vect256_t V1,T;
		VEC256_LOADU(V1,ABCD);
		VEC256_MOD_P(V1,V1,P,T);
		VEC256_STOREU(ABCD,V1);
	}


	/*---------------------------------------------------*/
	/*--  implementation of DIF with 256-bits AVX    ----*/
	/*---------------------------------------------------*/
 
	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_8x1_AVX(Element* ABCDEFGH, Element* IJKLMNOP,
								      const Element* alpha,
								      const Element* alphap,
								      const _vect256_t& P, const _vect256_t& P2) {
		_vect256_t V1,V2,V3,V4,W,Wp,T;
		// V1=[A B C D E F G H], V2=[I J K L M N O P]
		VEC256_LOADU(V1,ABCDEFGH);
		VEC256_LOADU(V2,IJKLMNOP);
		VEC256_LOADU(W ,alpha);
		VEC256_LOADU(Wp,alphap);
	
		// V3 = V1 + V2 mod 
		VEC256_ADD_MOD(V3,V1,V2,P2,T);
		VEC256_STOREU(ABCDEFGH,V3);
	
		// V4 = (V1+(2P-V2))alpha mod 2P
		VEC256_SUB_32(T,V2,P2);
		VEC256_SUB_32(V4,V1,T);	
		VEC256_MUL_MOD(T,V4,W,P,Wp,V1,V2,V3);// V3 is the result
		VEC256_STOREU(IJKLMNOP,T);
	}


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_8x3_AVX_last3step(Element* ABCDEFGH, Element* IJKLMNOP,
										const _vect256_t& alpha,const _vect256_t& alphap,
										const _vect256_t& beta ,const _vect256_t& betap,
										const _vect256_t& P    ,const _vect256_t& P2) {
		_vect256_t V1,V2,V3,V4,V5,V6,V7,Q;

		// V1=[A B C D E F G H], V2=[I J K L M N O P]
		VEC256_LOADU(V1,ABCDEFGH);
		VEC256_LOADU(V2,IJKLMNOP);

		/* 1st step */
		// V3=[A B C D I J K L] V4=[E F G H M N O P]
		VEC256_UNPACK_LO_128(V3,V1,V2);
		VEC256_UNPACK_HI_128(V4,V1,V2);

		// V1 = V3 + V4 mod 2P
		// P2 = [2p 2p 2p 2p]
		VEC256_ADD_MOD(V1,V3,V4,P2,V5);

		// V2 = (V3+(2P-V4))alpha mod 2P
		VEC256_SUB_32(V5,V4,P2);
		VEC256_SUB_32(V6,V3,V5);
		VEC256_MOD_P(V7,V6,P2,V2);
		VEC256_MUL_MOD(V2,V7,alpha,P,alphap,V3,V4,V5);

		/* 2nd step */

		// V3=[A E B F I M J N] V4=[C G D H K O L P]
		VEC256_UNPACK_LO_32(V3,V1,V2);
		VEC256_UNPACK_HI_32(V4,V1,V2);

		// V1 = V3 + V4 mod 2P
		// P2 = [2p 2p 2p 2p]
		VEC256_ADD_MOD(V1,V3,V4,P2,V5);

		// V2 = (V3+(2P-V4))alpha mod 2P
		// V7 =  (V3+(2P-V4)) mod 2P
		VEC256_SUB_32(V5,V4,P2);
		VEC256_SUB_32(V6,V3,V5);
		VEC256_MOD_P(V7,V6,P2,V2);
	
		// V4 = [D D H H L L P P ]
		VEC256_UNPACK_HI_32(V4,V7,V7);

		// Q = V4 * beta mod 2^64 = [* Qd * Qh * Ql * Qp]
		// with betap= [ betap * betap * betap * betap *]
		VEC256_MUL_32(Q,V4,betap);	
		// V5 = [* Qd.P * Qh.P * Ql.P * Qp.P]
		VEC256_MUL_LO_32(V5,Q,P);	
		// V6 = V4 * beta mod 2^32
		VEC256_MUL_LO_32(V6,V4,beta);
		// V3 = V6 - V5 = [* (D.beta mod p) * (H.beta mod p) * (L.beta mod p) * (P.beta mod p)]
		VEC256_SUB_32(V3,V6,V5);
		// V2=[* * D H * * L P]
		VEC256_SHUFFLE_32(V2,V3,0xDD);
 
		/* 3nd step */
		// At this time I have V1=[A B E F I J M N], V7=[C G * * K O * *], V2=[* * D H * * L P]
		// I need V3 = [A C E G I K M O], V4=[B D F H J L N P]
		VEC256_UNPACK_LO_32(V3,V1,V7);
		VEC256_UNPACK_HI_32(V4,V1,V2);

		// V1 = V3 + V4 mod 2P
		VEC256_ADD_MOD(V1,V3,V4,P2,V5);

		// V2 = V3 + (2P - V4) mod 2P
		VEC256_SUB_32(V5,V4,P2);
		VEC256_SUB_32(V6,V3,V5);
		VEC256_MOD_P(V2,V6,P2,V2);

		// Result in    V1=[A C E G I K M O] V2=[B D F H J L N P]
		// Transform to V3=[A B C D I J K L],V4=[E F G H M N O P]
		VEC256_UNPACK_LO_32(V3,V1,V2);
		VEC256_UNPACK_HI_32(V4,V1,V2);

		// Transform to V1=[A B C D E F G H], V2=[I J K L M N O P]
		VEC256_UNPACK_LO_128(V1,V3,V4);
		VEC256_UNPACK_HI_128(V2,V3,V4);

		// Store
		VEC256_STOREU(ABCDEFGH,V1);
		VEC256_STOREU(IJKLMNOP,V2);
	}



	template <class Field>
	template <class Polynomial>
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative8x1_AVX (Polynomial &fft) {
		_vect256_t P,P2;
		VEC256_SET_32(P,_pl);
		VEC256_SET_32(P2,_dpl);
	
		Element * tab_w = &pow_w [0];
		Element * tab_wp= &pow_wp[0];
		size_t w, f;
		for (w = n >> 1, f = 1; w >= 8; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1){		
			// w : witdh of butterflies
			// f : # families of butterflies
			for (size_t i = 0; i < f; i++)
				for (size_t j = 0; j < w; j+=8)				
#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
					Butterfly_DIF_mod2p_8x1_AVX(A0,A4, tab_w+j,tab_wp+j,P,P2);
		
#undef A0
#undef A4
			//std::cout<<fft<<std::endl;
		}
		// Last three steps
		//if (false) {
		if (n >= 16) {
			_vect256_t alpha,alphap,beta,betap;
			Element tmp[8];		
			tmp[0]=tmp[4]=tab_w[0];
			tmp[1]=tmp[5]=tab_w[1];
			tmp[2]=tmp[6]=tab_w[2];
			tmp[3]=tmp[7]=tab_w[3];
			VEC256_LOADU(alpha,tmp);
			tmp[0]=tmp[4]=tab_wp[0];
			tmp[1]=tmp[5]=tab_wp[1];
			tmp[2]=tmp[6]=tab_wp[2];
			tmp[3]=tmp[7]=tab_wp[3];
			VEC256_LOADU(alphap,tmp);
			VEC256_SET_32(beta,tab_w [5]);
			VEC256_SET_32(betap,tab_wp [5]);
			
			for (size_t i = 0; i < f; i+=2)
#define A0 &fft[0] + (i << 3)
#define A4 &fft[0] + (i << 3)+8
				Butterfly_DIF_mod2p_8x3_AVX_last3step(A0,A4,alpha,alphap,beta,betap,P,P2);
#undef A0
#undef A4
			//std::cout<<fft<<std::endl;
		} else {
			for (; w >= 1; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1)
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						Butterfly_DIF_mod2p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
		}
	}


	/*---------------------------------------------------*/
	/*--  implementation of DIF with 256-bits AVX    ----*/
	/*---------------------------------------------------*/

	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p_8x1_AVX(Element* ABCDEFGH, Element* IJKLMNOP,
								      const Element* alpha,
								      const Element* alphap,
								      const _vect256_t& P, const _vect256_t& P2) {
		_vect256_t V1,V2,V3,V4,W,Wp,T1,T2;
		// V1=[A B C D E F G H], V2=[I J K L M N O P]
		VEC256_LOADU(V1,ABCDEFGH);
		VEC256_LOADU(V2,IJKLMNOP);
		VEC256_LOADU(W ,alpha);
		VEC256_LOADU(Wp,alphap);
	
		// V3 = V1 mod 2P
		VEC256_MOD_P (V3,V1,P2,T1);
	
		// V4 = V2 * W mod P
		VEC256_MUL_MOD(V4,V2,W,P,Wp,V1,T1,T2);
	
		// V1 = V3 + V4
		VEC256_ADD_32(V1,V3,V4);
		VEC256_STOREU(ABCDEFGH,V1);
	
		// V2 = V3 - (V4 - 2P)
		VEC256_SUB_32(T1,V4,P2);
		VEC256_SUB_32(V2,V3,T1);
		VEC256_STOREU(IJKLMNOP,V2);
	}


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p_8x3_AVX_first3step(Element* ABCDEFGH, Element* IJKLMNOP,
										 const _vect256_t& alpha,const _vect256_t& alphap,
										 const _vect256_t& beta ,const _vect256_t& betap,
										 const _vect256_t& P, const _vect256_t& P2) {
		_vect256_t V1,V2,V3,V4,V5,V6,V7,Q;
		// V1=[A B C D E F G H], V2=[I J K L M N O P]
		VEC256_LOADU(V1,ABCDEFGH);
		VEC256_LOADU(V2,IJKLMNOP);


		/*********************************************/
		/* 1st STEP */
		/*********************************************/
		// Transform to V3=[A I C K E M G O], V4=[B J D L F N H P]
		VEC256_UNPACK_LO_32(V6,V1,V2); // V6=[A I B J E M F N]
		VEC256_UNPACK_HI_32(V7,V1,V2); // V7=[C K D L G O H P]
		VEC256_UNPACK_LO_64(V3,V6,V7); // V3=[A I C K E M G O]
		VEC256_UNPACK_HI_64(V4,V6,V7); // V4=[B J D L F N H P]




		// V1 = V3 + V4;       V1 = [A I C K E M G O]
		// Rk: No need for (. mod 2P) since entries are <P
		VEC256_ADD_32(V1,V3,V4);
	
		// V2 = V3 + (P - V4); V2 = [B J D L F N H P]
		// Rk: No need for (. mod 2P) since entries are <P
		VEC256_SUB_32(V6,V4,P);
		VEC256_SUB_32(V2,V3,V6);
	
		/*********************************************/
		/* 2nd STEP */
		/*********************************************/
		// V5 = [D D L L H H P P]
		VEC256_UNPACK_HI_32(V5,V2,V2); 
		// Q = V5 * alpha mod 2^64 = [* Qd * Qh * Ql * Qp]
		// with betap= [ alphap * alphap * alphap * alphap *]
		VEC256_MUL_32(Q,V5,alphap);	
		// V6 = [* Qd.P * Qh.P * Ql.P * Qp.P]
		VEC256_MUL_LO_32(V6,Q,P);	
		// V7 = V5 * alpha mod 2^32
		VEC256_MUL_LO_32(V7,V5,alpha);
		// V3 = V7 - V6 = [* (D.alpha mod p) * (L.alpha mod p) * (H.alpha mod p) * (P.alpha mod p)]
		VEC256_SUB_32(V3,V7,V6);
		// V7=[D L * * H P * *]
		VEC256_SHUFFLE_32(V7,V3,0xFD);
		// V6 = [B J D L F N H P]
		VEC256_UNPACK_LO_64(V6,V2,V7);
		// V3= [A B I J E F M N], V4=[C D K L G H O P]
		VEC256_UNPACK_LO_32(V3,V1,V6);
		VEC256_UNPACK_HI_32(V4,V1,V6);

		// V1 = V3+V4
		VEC256_ADD_32(V1,V3,V4);	
		// V2 = V3 - (V4 - 2P)
		VEC256_SUB_32(V7,V4,P2);
		VEC256_SUB_32(V2,V3,V7);

		/*********************************************/
		/* 3nd STEP */
		/*********************************************/
		// V3= [A B C D I J K L] V4= [E F G H M N O P]
		VEC256_UNPACK_LO_64(V6,V1,V2);
		VEC256_UNPACK_HI_64(V7,V1,V2);
		VEC256_UNPACK_LO_128(V3,V6,V7);
		VEC256_UNPACK_HI_128(V4,V6,V7);

		// V6= V3 mod 2P
		VEC256_MOD_P(V6,V3,P2,V7);

		// V7= V4.beta mod p 
		VEC256_MUL_MOD(V7,V4,beta,P,betap,V1,V2,V5);

		// V1 = V6+V7
		VEC256_ADD_32(V1,V6,V7);

		// V2 = V6 - (V7 - 2P)
		VEC256_SUB_32(V5,V7,P2);
		VEC256_SUB_32(V2,V6,V5);
	
		/*********************************************/
		// V3=[A B C D E F G H] V4=[I J K L M N O P]
		VEC256_UNPACK_LO_128(V3,V1,V2);
		VEC256_UNPACK_HI_128(V4,V1,V2);

		// Store
		VEC256_STOREU(ABCDEFGH,V3);
		VEC256_STOREU(IJKLMNOP,V4);	
	}



	template <class Field>
	template <class Polynomial>
	void FFT_transform<Field>::FFT_DIT_Harvey_mod4p_iterative8x1_AVX (Polynomial &fft) {
		_vect256_t P,P2;
		VEC256_SET_32(P,_pl);
		VEC256_SET_32(P2,_dpl);

		// first three steps	
		if (n >= 16) {
			_vect256_t alpha,alphap,beta,betap;
			VEC256_SET_32(alpha,pow_w[n-3]);
			VEC256_SET_32(alphap,pow_wp[n-3]);
			Element tmp[8];		
			tmp[0]=tmp[4]=pow_w[n-8];
			tmp[1]=tmp[5]=pow_w[n-7];
			tmp[2]=tmp[6]=pow_w[n-6];
			tmp[3]=tmp[7]=pow_w[n-5];
			VEC256_LOADU(beta,tmp);
			tmp[0]=tmp[4]=pow_wp[n-8];
			tmp[1]=tmp[5]=pow_wp[n-7];
			tmp[2]=tmp[6]=pow_wp[n-6];
			tmp[3]=tmp[7]=pow_wp[n-5];
			VEC256_LOADU(betap,tmp);
			for (size_t i = 0; i < n; i+=16)
				Butterfly_DIT_mod4p_8x3_AVX_first3step(&fft[i],&fft[i+8],alpha,alphap,beta,betap,P,P2);
			Element * tab_w = &pow_w [n-16];
			Element * tab_wp= &pow_wp[n-16];
			for (size_t w = 8, f = n >> 4; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w){
				// w : witdh of butterflies
				// f : # families of butterflies
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j+=8)
#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
						Butterfly_DIT_mod4p_8x1_AVX(A0,A4, tab_w+j,tab_wp+j,P,P2);
#undef A0
#undef A4
			}
		} else {
			Element * tab_w = &pow_w [n-2];
			Element * tab_wp= &pow_wp[n-2];
			for (size_t w = 1, f = n >> 1; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w)
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j++)
						Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
		}
	}



} // enf of namespace LinBox

#endif // end of AVX2 section

#endif //end of file
