/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#include "fflas-ffpack/fflas/fflas_simd.h"

//#include "linbox/algorithms/polynomial-matrix/simd.h"

#ifndef additional_modular_simd_functions
#define additional_modular_simd_functions

#define Simd_vect typename Simd::vect_t

template <class Simd>
inline Simd_vect reduce (const Simd_vect a, const Simd_vect p) {
	Simd_vect t = Simd::greater(p,a);
	return Simd::sub(a, Simd::vandnot(p,t));
}

template <class Simd>
inline Simd_vect add_mod (const Simd_vect a, const Simd_vect b, const Simd_vect p) {
	Simd_vect c = Simd::add(a,b);
	return reduce<Simd>(c, p);
}

template <class Simd>
inline Simd_vect mul_mod (const Simd_vect a, const Simd_vect b, const Simd_vect p, const Simd_vect bp) {
	Simd_vect q = Simd::mulhi(a,bp);
	Simd_vect c = Simd::mullo(a,b);
	Simd_vect t = Simd::mullo(q,p);
	return Simd::sub(c,t);
}
#undef Simd_vect
#endif

namespace LinBox {

	
	/******************************************************************************************************************
	 ******************************************************************************************************************
	 ***********************************    FFT with SSE CODE   *******************************************************
	 ******************************************************************************************************************
	 ******************************************************************************************************************/

	template <class Field>
	inline void FFT_transform<Field>::reduce128_modp(uint32_t* ABCD, const _vect128_t& P) {
		_vect128_t V1;
		// V1=[A B C D], V2=[E F G H]
		V1 = Simd128<uint32_t>::load(ABCD);
		V1 = reduce<Simd128<uint32_t> >(V1, P);
		Simd128<uint32_t>::store(ABCD,V1);
	}

	/*-----------------------------------*/
	/*--        Butterflies DIF      ----*/
	/*-----------------------------------*/


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_4x1_SSE(uint32_t* ABCD, uint32_t* EFGH,
																  const uint32_t* alpha,
																  const uint32_t* alphap,
																  const _vect128_t& P, const _vect128_t& P2)
	{
		_vect128_t V1,V2,V3,V4,W,Wp,T;
		// V1=[A B C D], V2=[E F G H]
		V1 = Simd128<uint32_t>::load(ABCD);
		V2 = Simd128<uint32_t>::load(EFGH);
		W  = Simd128<uint32_t>::load(alpha);
		Wp = Simd128<uint32_t>::load(alphap);
		// V3 = V1 + V2 mod 2P
		V3 = add_mod<Simd128<uint32_t> >(V1,V2,P2);
		Simd128<uint32_t>::store(ABCD,V3);
		// V4 = (V1+(2P-V2))alpha mod 2P
		T = Simd128<uint32_t>::sub(V2,P2);
		V4 = Simd128<uint32_t>::sub(V1,T);
		T = mul_mod<Simd128<uint32_t> >(V4,W,P,Wp);// T is the result
		Simd128<uint32_t>::store(EFGH,T);
	}

/*
	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_4x1_SSE_laststep(uint32_t* ABCD, uint32_t* EFGH, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,V5;
		// V1=[A B C D], V2=[E F G H]
		V1 = Simd128<uint32_t>::load(ABCD);
		V2 = Simd128<uint32_t>::load(EFGH);
		// V3 = [A C B D], V4 = [E G F H]
		V3 = Simd128<uint32_t>::shuffle<0xD8>(V1);
		V4 = Simd128<uint32_t>::shuffle<0xD8>(V2);
		// V1 = [A E C G], V2 = [B F D H]
		V1 = Simd128<uint32_t>::unpacklo(V3,V4);
		V2 = Simd128<uint32_t>::unpackhi(V3,V4);
		// V3 = V1 + V2 mod 2P
		V3 = add_mod<Simd128<uint32_t> >(V1,V2,P2);
		// V4 = V1 + (2P - V2) mod 2P
		V5 = Simd128<uint32_t>::sub(V2,P2);
		V2 = Simd128<uint32_t>::sub(V1,V5);
		V4 = reduce<Simd128<uint32_t> >(V2, P2);
		// V1 = [A C E G], V2 = [B D F H]
		V1 = Simd128<uint32_t>::shuffle<0xD8>(V3);
		V2 = Simd128<uint32_t>::shuffle<0xD8>(V4);
		// V3 = [A B C D], V4 = [E F G H]
		V3 = Simd128<uint32_t>::unpacklo(V1,V2);
		V4 = Simd128<uint32_t>::unpackhi(V1,V2);
		// Store
		Simd128<uint32_t>::store(ABCD,V3);
		Simd128<uint32_t>::store(EFGH,V4);
	}
*/

	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_4x2_SSE(uint32_t* ABCD, uint32_t* EFGH, uint32_t* IJKL, uint32_t* MNOP,
																  const uint32_t* alpha,  const uint32_t*beta ,  const uint32_t* gamma,
																  const uint32_t* alphap, const uint32_t*betap , const uint32_t* gammap,
																  const _vect128_t& P, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,W,Wp,T1,T2,T3,T4,T6,T7;

		// V1=[A B C D], V2=[E F G H], V3=[I J K L], V4=[M N O P]
		V1 = Simd128<uint32_t>::load(ABCD);
		V2 = Simd128<uint32_t>::load(IJKL);
		W  = Simd128<uint32_t>::load(alpha);
		Wp = Simd128<uint32_t>::load(alphap);
		/**************/
		// T1 = V1 + V2 mod 2P
		T1 = add_mod<Simd128<uint32_t> >(V1,V2,P2);
		// T2 = (V1+(2P-V2))alpha mod 2P
		T7 = Simd128<uint32_t>::sub(V2,P2);
		T6 = Simd128<uint32_t>::sub(V1,T7);
		T2 = mul_mod<Simd128<uint32_t> >(T6,W,P,Wp);
		/**************/
		V3 = Simd128<uint32_t>::load(EFGH);
		V4 = Simd128<uint32_t>::load(MNOP);
		W  = Simd128<uint32_t>::load(beta);
		Wp = Simd128<uint32_t>::load(betap);
		/**************/
		// T3 = V3 + V4 mod 2P
		T3 = add_mod<Simd128<uint32_t> >(V3,V4,P2);
		// T4 = (V3+(2P-V4))beta mod 2P
		T7 = Simd128<uint32_t>::sub(V4,P2);
		T6 = Simd128<uint32_t>::sub(V3,T7);
		T4 = mul_mod<Simd128<uint32_t> >(T6,W,P,Wp);// T1 is the result
		/**************/
		W  = Simd128<uint32_t>::load(gamma);
		Wp = Simd128<uint32_t>::load(gammap);
		/**************/
		// V1 = T1 + T3 mod 2P
		V1 = add_mod<Simd128<uint32_t> >(T1,T3,P2);
		// V3 = (T1+(2P-T3))gamma mod 2P
		T7 = Simd128<uint32_t>::sub(T3,P2);
		T6 = Simd128<uint32_t>::sub(T1,T7);
		V3 = mul_mod<Simd128<uint32_t> >(T6,W,P,Wp);// T1 is the result
		/**************/
		// V2 = T2 + T4 mod 2P
		V2 = add_mod<Simd128<uint32_t> >(T2,T4,P2);
		// V4 = (T2+(2P-T4))gamma mod 2P
		T7 = Simd128<uint32_t>::sub(T4,P2);
		T6 = Simd128<uint32_t>::sub(T2,T7);
		V4 = mul_mod<Simd128<uint32_t> >(T6,W,P,Wp);// T1 is the result
		/**************/
		Simd128<uint32_t>::store(ABCD,V1);
		Simd128<uint32_t>::store(EFGH,V3);
		Simd128<uint32_t>::store(IJKL,V2);
		Simd128<uint32_t>::store(MNOP,V4);
	}



	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_4x2_SSE_last2step(uint32_t* ABCD, uint32_t* EFGH,
																			const _vect128_t& W,
																			const _vect128_t& Wp,
																			const _vect128_t& P, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,V5,V6,V7;
		// V1=[A B C D], V2=[E F G H]
		V1 = Simd128<uint32_t>::load(ABCD);
		V2 = Simd128<uint32_t>::load(EFGH);
		// V3=[A E B F], V4=[C G D H]
		V3 = Simd128<uint32_t>::unpacklo(V1,V2);
		V4 = Simd128<uint32_t>::unpackhi(V1,V2);
		// V1 = V3 + V4 mod 2P
		// P2 = [2p 2p 2p 2p]
		V1 = add_mod<Simd128<uint32_t> >(V3,V4,P2);
		// V2 = (V3+(2P-V4))alpha mod 2P
		V5 = Simd128<uint32_t>::sub(V4,P2);
		V6 = Simd128<uint32_t>::sub(V3,V5);
		V2 = reduce<Simd128<uint32_t> >(V6, P2);
		// V4 = [D D H H]
		V4 = Simd128<uint32_t>::unpackhi(V2,V2);
		// V6 = V4 * Wp mod 2^64
		// Wp = [Wp ? Wp ?]
		V7 = Simd128<uint64_t>::mulx(V4,Wp);
		V5 = Simd128<uint32_t>::mullo(V7,P);
		// At this point V4= [? Q_D*p ? Q_H*p]
		// V5 = [D D H H] * [W W W W] mod 2^32
		V6 = Simd128<uint32_t>::mullo(V4,W);
		V4 = Simd128<uint32_t>::sub(V6,V5);
		V3 = Simd128<uint32_t>::shuffle<0xDD>(V4);
		//At this point, V2 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]
		// At this time I have V1=[A E B F], V2=[C G ? ?], V3=[? ? D H]
		// I need V3 = [A C E G], V4 = [B D F H]
		V4 = Simd128<uint32_t>::unpackhi(V1,V3);
		V3 = Simd128<uint32_t>::unpacklo(V1,V2);
		// V1 = V3 + V4 mod 2P
		V1 = add_mod<Simd128<uint32_t> >(V3,V4,P2);
		// V2 = V3 + (2P - V4) mod 2P
		V5 = Simd128<uint32_t>::sub(V4,P2);
		V6 = Simd128<uint32_t>::sub(V3,V5);
		V2 = reduce<Simd128<uint32_t> >(V6, P2);
		// Result in V1 = [A C E G]  and V2 = [B D F H]
		// Transform to V3=[A B C D], V4=[E F G H]
		V3 = Simd128<uint32_t>::unpacklo(V1,V2);
		V4 = Simd128<uint32_t>::unpackhi(V1,V2);
		// Store
		Simd128<uint32_t>::store(ABCD,V3);
		Simd128<uint32_t>::store(EFGH,V4);
	}


	/*-----------------------------------*/
	/*--        Butterflies DIT      ----*/
	/*-----------------------------------*/


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p_4x1_SSE(uint32_t* ABCD, uint32_t* EFGH,
																  const uint32_t* alpha,
																  const uint32_t* alphap,
																  const _vect128_t& P, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,W,Wp,T1;
		// V1=[A B C D], V2=[E F G H]
		V1 = Simd128<uint32_t>::load(ABCD);
		V2 = Simd128<uint32_t>::load(EFGH);
		W  = Simd128<uint32_t>::load(alpha);
		Wp = Simd128<uint32_t>::load(alphap);
		// V3 = V1 mod 2P
		V3 = reduce<Simd128<uint32_t> >(V1, P2);
		// V4 = V2 * W mod P
		V4 = mul_mod<Simd128<uint32_t> >(V2,W,P,Wp);
		// V1 = V3 + V4
		V1 = Simd128<uint32_t>::add(V3,V4);
		Simd128<uint32_t>::store(ABCD,V1);
		// V2 = V3 - (V4 - 2P)
		T1 = Simd128<uint32_t>::sub(V4,P2);
		V2 = Simd128<uint32_t>::sub(V3,T1);
		Simd128<uint32_t>::store(EFGH,V2);
	}

	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p_4x2_SSE_first2step(uint32_t* ABCD, uint32_t* EFGH,
																			 const _vect128_t& W,
																			 const _vect128_t& Wp,
																			 const _vect128_t& P, const _vect128_t& P2) {
		_vect128_t V1,V2,V3,V4,T1,T2,T3,T4;
		// V1=[A B C D], V2=[E F G H]
		V1 = Simd128<uint32_t>::load(ABCD);
		V2 = Simd128<uint32_t>::load(EFGH);
		// T1 = [A C B D], T2 = [E G F H]
		T1 = Simd128<uint32_t>::shuffle<0xD8>(V1);
		T2 = Simd128<uint32_t>::shuffle<0xD8>(V2);
		// V1 = [A E C G], V2 = [B F D H]
		V1 = Simd128<uint32_t>::unpacklo(T1,T2);
		V2 = Simd128<uint32_t>::unpackhi(T1,T2);
		// V3 = V1 + V2
		// Rk: No need for (. mod 2P) since entries are <P
		V3 = Simd128<uint32_t>::add(V1,V2);
		// V4 = V1 + (P - V2)
		// Rk: No need for (. mod 2P) since entries are <P
		T1 = Simd128<uint32_t>::sub(V2,P);
		V4 = Simd128<uint32_t>::sub(V1,T1);
		// T1 = [D D H H]
		T1 = Simd128<uint32_t>::unpackhi(V4,V4);
		// T2 = T1 * Wp mod 2^64
		// Wp = [Wp ? Wp ?]
		T2 = Simd128<uint64_t>::mulx(T1,Wp);
		T3 = Simd128<uint32_t>::mullo(T2,P);
		// At this point T3= [? Q_D*p ? Q_H*p]
		// T4 = [D D H H] * [W W W W] mod 2^32
		T4 = Simd128<uint32_t>::mullo(T1,W);
		T1 = Simd128<uint32_t>::sub(T4,T3);
		T2 = Simd128<uint32_t>::shuffle<0xDD>(T1);
		//At this point, T2 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]
		// At this time I have V3=[A E C G], V4=[B F ? ?], T2=[? ? D H]
		// I need V1 = [A B E F], V2 = [C D G H]
		V1 = Simd128<uint32_t>::unpacklo(V3,V4);
		V2 = Simd128<uint32_t>::unpackhi(V3,T2);
		// T1 = V1 + V2
		T1 = Simd128<uint32_t>::add(V1,V2);
		// T2 = V1 - (V2 - 2P)
		T3 = Simd128<uint32_t>::sub(V2,P2);
		T2 = Simd128<uint32_t>::sub(V1,T3);
		// Result in T1 = [A B E F]  and T2 = [C D G H]
		// Transform to V1=[A C B D], V2=[E G F H]
		V1 = Simd128<uint32_t>::unpacklo(T1,T2);
		V2 = Simd128<uint32_t>::unpackhi(T1,T2);
		// Then T1=[A B C D], T2=[E F G H]
		T1 = Simd128<uint32_t>::shuffle<0xD8>(V1);
		T2 = Simd128<uint32_t>::shuffle<0xD8>(V2);
		// Store
		Simd128<uint32_t>::store(ABCD,T1);
		Simd128<uint32_t>::store(EFGH,T2);
	}

	/*-----------------------------------*/
	/*--       SSE FFT functions     ----*/
	/*-----------------------------------*/

	template <class Field>
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative4x1_SSE (uint32_t *fft) {
		_vect128_t P,P2;
		P  = Simd128<uint32_t>::set1(_pl);
		P2 = Simd128<uint32_t>::set1(_dpl);
		uint32_t * tab_w = &pow_w [0];
		uint32_t * tab_wp= &pow_wp[0];
		size_t w, f;
		for (w = n >> 1, f = 1; w >= 4; tab_w+=w, tab_wp+=w, w >>= 1, f <<= 1){
				// w : witdh of butterflies
				// f : # families of butterflies
				for (size_t i = 0; i < f; i++)
					for (size_t j = 0; j < w; j+=4)

#define A0 &fft[0] +  (i << 1)   *w+ j
#define A4 &fft[0] + ((i << 1)+1)*w+ j
						Butterfly_DIF_mod2p_4x1_SSE(A0,A4, tab_w+j,tab_wp+j,P,P2);
#undef A0
#undef A4
				//std::cout<<fft<<std::endl;
			}
		// Last two steps
		if (n >= 8) {
				_vect128_t W,Wp;
				W = Simd128<uint32_t>::set1 ((int)tab_w [1]);
				Wp= Simd128<uint32_t>::set1 ((int)tab_wp[1]);

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
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative4x2_SSE (uint32_t *fft) {
		size_t w, f;
		_vect128_t P,P2;
		P  = Simd128<uint32_t>::set1(_pl);
		P2 = Simd128<uint32_t>::set1(_dpl);
		uint32_t * tab_w =  &pow_w[0];
		uint32_t * tab_wp= &pow_wp[0];
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
				W = Simd128<uint32_t>::set1 ((int)tab_w [1]);
				Wp= Simd128<uint32_t>::set1 ((int)tab_wp[1]);

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
	void FFT_transform<Field>::FFT_DIT_Harvey_mod4p_iterative4x1_SSE (uint32_t *fft)
	{
		_vect128_t P,P2;
		P = Simd128<uint32_t>::set1(_pl);
		P2 = Simd128<uint32_t>::set1(_dpl);
		// Last two steps
		if (n >= 8) {
				_vect128_t W,Wp;
				W = Simd128<uint32_t>::set1 ((int)pow_w [n-3]);
				Wp= Simd128<uint32_t>::set1 ((int)pow_wp[n-3]);

				for (size_t i = 0; i < n; i+=8)
					Butterfly_DIT_mod4p_4x2_SSE_first2step(&fft[i],&fft[i+4],W,Wp,P,P2);

				uint32_t * tab_w = &pow_w [n-8];
				uint32_t * tab_wp= &pow_wp[n-8];
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
				uint32_t * tab_w = &pow_w [n-2];
				uint32_t * tab_wp= &pow_wp[n-2];
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


#ifdef __LINBOX_HAVE_AVX2_INSTRUCTIONS

	template <class Field>
	inline void FFT_transform<Field>::reduce256_modp(uint32_t* ABCD, const _vect256_t& P) {
		_vect256_t V1;
		V1 = Simd256<uint32_t>::loadu(ABCD);
		V1 = reduce<Simd256<uint32_t> >(V1, P);
		Simd256<uint32_t>::storeu(ABCD,V1);
	}


	/*---------------------------------------------------*/
	/*--  implementation of DIF with 256-bits AVX    ----*/
	/*---------------------------------------------------*/

	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_8x1_AVX(uint32_t* ABCDEFGH, uint32_t* IJKLMNOP,
																  const uint32_t* alpha,
																  const uint32_t* alphap,
																  const _vect256_t& P, const _vect256_t& P2) {
		_vect256_t V1,V2,V3,V4,W,Wp,T;
		// V1=[A B C D E F G H], V2=[I J K L M N O P]
		V1 = Simd256<uint32_t>::loadu(ABCDEFGH);
		V2 = Simd256<uint32_t>::loadu(IJKLMNOP);
		W  = Simd256<uint32_t>::loadu(alpha);
		Wp = Simd256<uint32_t>::loadu(alphap);

		// V3 = V1 + V2 mod

		V3 = add_mod<Simd256<uint32_t> >(V1,V2,P2);

		Simd256<uint32_t>::storeu(ABCDEFGH,V3);

		// V4 = (V1+(2P-V2))alpha mod 2P
		T = Simd256<uint32_t>::sub(V2,P2);
		V4 = Simd256<uint32_t>::sub(V1,T);
		T = mul_mod<Simd256<uint32_t> >(V4,W,P,Wp);// T is the result
		Simd256<uint32_t>::storeu(IJKLMNOP,T);
	}


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIF_mod2p_8x3_AVX_last3step(uint32_t* ABCDEFGH, uint32_t* IJKLMNOP,
																			const _vect256_t& alpha,const _vect256_t& alphap,
																			const _vect256_t& beta ,const _vect256_t& betap,
																			const _vect256_t& P    ,const _vect256_t& P2) {
		_vect256_t V1,V2,V3,V4,V5,V6,V7,Q;

		// V1=[A B C D E F G H], V2=[I J K L M N O P]
		V1 = Simd256<uint32_t>::loadu(ABCDEFGH);
		V2 = Simd256<uint32_t>::loadu(IJKLMNOP);

		/* 1st step */
		// V3=[A B C D I J K L] V4=[E F G H M N O P]
		V3 = Simd256<uint64_t>::unpacklo128(V1,V2);
		V4 = Simd256<uint64_t>::unpackhi128(V1,V2);

		// V1 = V3 + V4 mod 2P
		// P2 = [2p 2p 2p 2p]
		V1 = add_mod<Simd256<uint32_t> >(V3,V4,P2);

		// V2 = (V3+(2P-V4))alpha mod 2P
		V5 = Simd256<uint32_t>::sub(V4,P2);
		V6 = Simd256<uint32_t>::sub(V3,V5);
		V7 = reduce<Simd256<uint32_t> >(V6, P2);
		V2 = mul_mod<Simd256<uint32_t> >(V7,alpha,P,alphap);

		/* 2nd step */

		// V3=[A E B F I M J N] V4=[C G D H K O L P]
		V3 = Simd256<uint32_t>::unpacklo_twice(V1,V2);
		V4 = Simd256<uint32_t>::unpackhi_twice(V1,V2);

		// V1 = V3 + V4 mod 2P
		// P2 = [2p 2p 2p 2p]
		V1 = add_mod<Simd256<uint32_t> >(V3,V4,P2);

		// V2 = (V3+(2P-V4))alpha mod 2P
		// V7 =  (V3+(2P-V4)) mod 2P
		V5 = Simd256<uint32_t>::sub(V4,P2);
		V6 = Simd256<uint32_t>::sub(V3,V5);
		V7 = reduce<Simd256<uint32_t> >(V6, P2);

		// V4 = [D D H H L L P P ]
		V4 = Simd256<uint32_t>::unpackhi_twice(V7,V7);

		// Q = V4 * beta mod 2^64 = [* Qd * Qh * Ql * Qp]
		// with betap= [ betap * betap * betap * betap *]
		Q = Simd256<uint64_t>::mulx(V4,betap);
		// V5 = [* Qd.P * Qh.P * Ql.P * Qp.P]
		V5 = Simd256<uint32_t>::mullo(Q,P);
		// V6 = V4 * beta mod 2^32
		V6 = Simd256<uint32_t>::mullo(V4,beta);
		// V3 = V6 - V5 = [* (D.beta mod p) * (H.beta mod p) * (L.beta mod p) * (P.beta mod p)]
		V3 = Simd256<uint32_t>::sub(V6,V5);
		// V2=[* * D H * * L P]
		V2 = Simd256<uint32_t>::shuffle_twice<0xDD>(V3);

		/* 3nd step */
		// At this time I have V1=[A B E F I J M N], V7=[C G * * K O * *], V2=[* * D H * * L P]
		// I need V3 = [A C E G I K M O], V4=[B D F H J L N P]
		V3 = Simd256<uint32_t>::unpacklo_twice(V1,V7);
		V4 = Simd256<uint32_t>::unpackhi_twice(V1,V2);

		// V1 = V3 + V4 mod 2P
		V1 = add_mod<Simd256<uint32_t> >(V3,V4,P2);

		// V2 = V3 + (2P - V4) mod 2P
		V5 = Simd256<uint32_t>::sub(V4,P2);
		V6 = Simd256<uint32_t>::sub(V3,V5);
		V2 = reduce<Simd256<uint32_t> >(V6, P2);

		// Result in    V1=[A C E G I K M O] V2=[B D F H J L N P]
		// Transform to V3=[A B C D I J K L],V4=[E F G H M N O P]
		V3 = Simd256<uint32_t>::unpacklo_twice(V1,V2);
		V4 = Simd256<uint32_t>::unpackhi_twice(V1,V2);

		// Transform to V1=[A B C D E F G H], V2=[I J K L M N O P]
		V1 = Simd256<uint64_t>::unpacklo128(V3,V4);
		V2 = Simd256<uint64_t>::unpackhi128(V3,V4);

		// Store
		Simd256<uint32_t>::storeu(ABCDEFGH,V1);
		Simd256<uint32_t>::storeu(IJKLMNOP,V2);
	}



	template <class Field>
	void FFT_transform<Field>::FFT_DIF_Harvey_mod2p_iterative8x1_AVX (uint32_t *fft) {
		_vect256_t P,P2;
		P = Simd256<uint32_t>::set1(_pl);
		P2 = Simd256<uint32_t>::set1(_dpl);

		uint32_t * tab_w = &pow_w [0];
		uint32_t * tab_wp= &pow_wp[0];
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
		if (n >= 16) {
				_vect256_t alpha,alphap,beta,betap;
				uint32_t tmp[8];
				tmp[0]=tmp[4]=tab_w[0];
				tmp[1]=tmp[5]=tab_w[1];
				tmp[2]=tmp[6]=tab_w[2];
				tmp[3]=tmp[7]=tab_w[3];
				alpha = Simd256<uint32_t>::loadu(tmp);
				tmp[0]=tmp[4]=tab_wp[0];
				tmp[1]=tmp[5]=tab_wp[1];
				tmp[2]=tmp[6]=tab_wp[2];
				tmp[3]=tmp[7]=tab_wp[3];
				alphap = Simd256<uint32_t>::loadu(tmp);
				beta = Simd256<uint32_t>::set1(tab_w [5]);
				betap = Simd256<uint32_t>::set1(tab_wp [5]);

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
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p_8x1_AVX(uint32_t* ABCDEFGH, uint32_t* IJKLMNOP,
																  const uint32_t* alpha,
																  const uint32_t* alphap,
																  const _vect256_t& P, const _vect256_t& P2) {
		_vect256_t V1,V2,V3,V4,W,Wp,T1;
		// V1=[A B C D E F G H], V2=[I J K L M N O P]
		V1 = Simd256<uint32_t>::loadu(ABCDEFGH);
		V2 = Simd256<uint32_t>::loadu(IJKLMNOP);
		W  = Simd256<uint32_t>::loadu(alpha);
		Wp = Simd256<uint32_t>::loadu(alphap);

		// V3 = V1 mod 2P
		V3 = reduce<Simd256<uint32_t> >(V1, P2);

		// V4 = V2 * W mod P
		V4 = mul_mod<Simd256<uint32_t> >(V2,W,P,Wp);

		// V1 = V3 + V4
		V1 = Simd256<uint32_t>::add(V3,V4);
		Simd256<uint32_t>::storeu(ABCDEFGH,V1);

		// V2 = V3 - (V4 - 2P)
		T1 = Simd256<uint32_t>::sub(V4,P2);
		V2 = Simd256<uint32_t>::sub(V3,T1);
		Simd256<uint32_t>::storeu(IJKLMNOP,V2);
	}


	template <class Field>
	inline void FFT_transform<Field>::Butterfly_DIT_mod4p_8x3_AVX_first3step(uint32_t* ABCDEFGH, uint32_t* IJKLMNOP,
																			 const _vect256_t& alpha,const _vect256_t& alphap,
																			 const _vect256_t& beta ,const _vect256_t& betap,
																			 const _vect256_t& P, const _vect256_t& P2) {
		_vect256_t V1,V2,V3,V4,V5,V6,V7,Q;
		// V1=[A B C D E F G H], V2=[I J K L M N O P]
		V1 = Simd256<uint32_t>::loadu(ABCDEFGH);
		V2 = Simd256<uint32_t>::loadu(IJKLMNOP);


		/*********************************************/
		/* 1st STEP */
		/*********************************************/
		// Transform to V3=[A I C K E M G O], V4=[B J D L F N H P]
		V6 = Simd256<uint32_t>::unpacklo_twice(V1,V2); // V6=[A I B J E M F N]
		V7 = Simd256<uint32_t>::unpackhi_twice(V1,V2); // V7=[C K D L G O H P]
		V3 = Simd256<uint64_t>::unpacklo_twice(V6,V7); // V3=[A I C K E M G O]
		V4 = Simd256<uint64_t>::unpackhi_twice(V6,V7); // V4=[B J D L F N H P]




		// V1 = V3 + V4;       V1 = [A I C K E M G O]
		// Rk: No need for (. mod 2P) since entries are <P
		V1 = Simd256<uint32_t>::add(V3,V4);

		// V2 = V3 + (P - V4); V2 = [B J D L F N H P]
		// Rk: No need for (. mod 2P) since entries are <P
		V6 = Simd256<uint32_t>::sub(V4,P);
		V2 = Simd256<uint32_t>::sub(V3,V6);

		/*********************************************/
		/* 2nd STEP */
		/*********************************************/
		// V5 = [D D L L H H P P]
		V5 = Simd256<uint32_t>::unpackhi_twice(V2,V2);
		// Q = V5 * alpha mod 2^64 = [* Qd * Qh * Ql * Qp]
		// with betap= [ alphap * alphap * alphap * alphap *]
		Q = Simd256<uint64_t>::mulx(V5,alphap);
		// V6 = [* Qd.P * Qh.P * Ql.P * Qp.P]
		V6 = Simd256<uint32_t>::mullo(Q,P);
		// V7 = V5 * alpha mod 2^32
		V7 = Simd256<uint32_t>::mullo(V5,alpha);
		// V3 = V7 - V6 = [* (D.alpha mod p) * (L.alpha mod p) * (H.alpha mod p) * (P.alpha mod p)]
		V3 = Simd256<uint32_t>::sub(V7,V6);
		// V7=[D L * * H P * *]
		V7 = Simd256<uint32_t>::shuffle_twice<0xFD>(V3);
		// V6 = [B J D L F N H P]
		V6 = Simd256<uint64_t>::unpacklo_twice(V2,V7);
		// V3= [A B I J E F M N], V4=[C D K L G H O P]
		V3 = Simd256<uint32_t>::unpacklo_twice(V1,V6);
		V4 = Simd256<uint32_t>::unpackhi_twice(V1,V6);

		// V1 = V3+V4
		V1 = Simd256<uint32_t>::add(V3,V4);
		// V2 = V3 - (V4 - 2P)
		V7 = Simd256<uint32_t>::sub(V4,P2);
		V2 = Simd256<uint32_t>::sub(V3,V7);

		/*********************************************/
		/* 3nd STEP */
		/*********************************************/
		// V3= [A B C D I J K L] V4= [E F G H M N O P]
		V6 = Simd256<uint64_t>::unpacklo_twice(V1,V2);
		V7 = Simd256<uint64_t>::unpackhi_twice(V1,V2);
		V3 = Simd256<uint64_t>::unpacklo128(V6,V7);
		V4 = Simd256<uint64_t>::unpackhi128(V6,V7);

		// V6= V3 mod 2P
		V6 = reduce<Simd256<uint32_t> >(V3, P2);

		// V7= V4.beta mod p
		V7 = mul_mod<Simd256<uint32_t> >(V4,beta,P,betap);

		// V1 = V6+V7
		V1 = Simd256<uint32_t>::add(V6,V7);

		// V2 = V6 - (V7 - 2P)
		V5 = Simd256<uint32_t>::sub(V7,P2);
		V2 = Simd256<uint32_t>::sub(V6,V5);

		/*********************************************/
		// V3=[A B C D E F G H] V4=[I J K L M N O P]
		V3 = Simd256<uint64_t>::unpacklo128(V1,V2);
		V4 = Simd256<uint64_t>::unpackhi128(V1,V2);

		// Store
		Simd256<uint32_t>::storeu(ABCDEFGH,V3);
		Simd256<uint32_t>::storeu(IJKLMNOP,V4);
	}



	template <class Field>
	void FFT_transform<Field>::FFT_DIT_Harvey_mod4p_iterative8x1_AVX (uint32_t *fft) {
		_vect256_t P,P2;
		P = Simd256<uint32_t>::set1(_pl);
		P2 = Simd256<uint32_t>::set1(_dpl);

		// first three steps
		if (n >= 16) {
				_vect256_t alpha,alphap,beta,betap;
				alpha = Simd256<uint32_t>::set1(pow_w[n-3]);
				alphap = Simd256<uint32_t>::set1(pow_wp[n-3]);
				uint32_t tmp[8];
				tmp[0]=tmp[4]=pow_w[n-8];
				tmp[1]=tmp[5]=pow_w[n-7];
				tmp[2]=tmp[6]=pow_w[n-6];
				tmp[3]=tmp[7]=pow_w[n-5];
				beta = Simd256<uint32_t>::loadu(tmp);
				tmp[0]=tmp[4]=pow_wp[n-8];
				tmp[1]=tmp[5]=pow_wp[n-7];
				tmp[2]=tmp[6]=pow_wp[n-6];
				tmp[3]=tmp[7]=pow_wp[n-5];
				betap = Simd256<uint32_t>::loadu(tmp);
				for (uint64_t i = 0; i < n; i+=16)
					Butterfly_DIT_mod4p_8x3_AVX_first3step(&fft[i],&fft[i+8],alpha,alphap,beta,betap,P,P2);
				uint32_t * tab_w = &pow_w [n-16];
				uint32_t * tab_wp= &pow_wp[n-16];
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
				uint32_t * tab_w = &pow_w [n-2];
				uint32_t * tab_wp= &pow_wp[n-2];
				for (size_t w = 1, f = n >> 1; f >= 1; w <<= 1, f >>= 1, tab_w-=w, tab_wp-=w)
					for (size_t i = 0; i < f; i++)
						for (size_t j = 0; j < w; j++)
							Butterfly_DIT_mod4p(fft[(i << 1)*w+j], fft[((i << 1)+1)*w+j], tab_w[j], tab_wp[j]);
			}
	}




#endif // end of AVX2 section

} // enf of namespace LinBox

#endif //end of file
