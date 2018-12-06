/*
 * Copyright (C) 2016 Romain Lebreton, Pascal Giorgi
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


#ifndef __LINBOX_polynomial_fft_butterflies_H
#define __LINBOX_polynomial_fft_butterflies_H

#include <iostream>
#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include "fflas-ffpack/fflas/fflas_simd.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-init.h"
#include "linbox/algorithms/polynomial-matrix/simd-additional-functions.h"

#define IS_INTEGRAL \
    typename std::enable_if<std::is_integral<typename Field::Element>::value>::type
#define IS_FLOATING \
    typename std::enable_if<std::is_floating_point<typename Field::Element>::value>::type

namespace LinBox {

	template<typename Field, typename simd = Simd<typename Field::Element>, uint8_t byn = simd::vect_size, typename Enable=void>
	class FFT_butterflies : public FFT_init<Field> {
	public:
		FFT_butterflies(const FFT_init<Field>& f_i) : FFT_init<Field>(f_i) {
			std::cerr<<"Not implemented !\n";
		}
	}; // FFT_butterflies

	template<typename Field>
	class FFT_butterflies<Field, NoSimd<typename Field::Element>, 1, IS_INTEGRAL> : public FFT_init<Field> {
	public:

		using Element = typename Field::Element;

		FFT_butterflies(const FFT_init<Field>& f_i) : FFT_init<Field>(f_i) {}

		inline void Butterfly_DIT_mod4p (Element& A, Element& B, const Element& alpha, const Element& alphap) {
			using Compute_t = typename Field::Compute_t;
			// Harvey's algorithm
			// 0 <= A,B < 4*p, p < 2^32 / 4
			// alphap = Floor(alpha * 2^ 32 / p])

			// TODO : replace by substract if greater
			if (A >= this->_dpl) A -= this->_dpl;

			// TODO : replace by mul_mod_shoup
			Element tmp = ((Element) alphap * (Compute_t)B) >> (8*sizeof(Element));
			tmp = alpha * B - tmp * this->_pl;

			// TODO : replace by add_r and sub_r
			B = A + (this->_dpl - tmp);
			//        B &= 0XFFFFFFFF;
			A += tmp;
		}

		inline void Butterfly_DIF_mod2p (Element& A, Element& B, const Element& alpha, const Element& alphap) {
			//std::cout<<A<<" $$ "<<B<<"("<<alpha<<","<<alphap<<" ) -> ";
			using Compute_t = typename Field::Compute_t;
			// Harvey's algorithm
			// 0 <= A,B < 2*p, p < 2^32 / 4
			// alphap = Floor(alpha * 2^ 32 / p])

			Element tmp = A;

			A += B;

			if (A >= this->_dpl) A -= this->_dpl;

			B = tmp + (this->_dpl - B);

			tmp = ((Element) alphap * (Compute_t) B) >> (8*sizeof(Element));
			B = alpha * B - tmp * this->_pl;
			//B &= 0xFFFFFFFF;
			//std::cout<<A<<" $$ "<<B<<"\n ";
		}

	}; // FFT_butterflies<Field, 1>

	// ATTENTION à tous les uint64_t, SimdComp restants !!!!

	template<typename Field, typename simd>
	class FFT_butterflies<Field, simd, 4, IS_INTEGRAL> : public FFT_init<Field> {
	public:

		using Element = typename Field::Element;
		using vect_t = typename simd::vect_t;
		using SimdComp = typename SimdCompute_t<simd,Field>::Compute_t;

		FFT_butterflies(const FFT_init<Field>& f_i) : FFT_init<Field>(f_i) {
			linbox_check(simd::vect_size == 4);
		}

		// TODO include P, P2 in precomp
		// TODO : Same functions Butterfly_DIT_mod4p Butterfly_DIF_mod2p in FFT_butterflies<Field, 8>
		inline void Butterfly_DIT_mod4p (Element* ABCD, Element* EFGH,
										 const Element* alpha, const Element* alphap,
										 const vect_t& P, const vect_t& P2) {
			vect_t V1,V2,V3,V4,W,Wp,T1;
			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = MemoryOp<Element,simd>::load(ABCD);
			V2 = MemoryOp<Element,simd>::load(EFGH);
			W  = MemoryOp<Element,simd>::load(alpha);
			Wp = MemoryOp<Element,simd>::load(alphap);

			// V3 = V1 mod 2P
			V3 = reduce<simd>(V1, P2);

			// V4 = V2 * W mod P
			V4 = mul_mod<simd>(V2,W,P,Wp);

			// V1 = V3 + V4
			V1 = simd::add(V3,V4);
			MemoryOp<Element,simd>::store(ABCD,V1);

			// V2 = V3 - (V4 - 2P)
			T1 = simd::sub(V4,P2);
			V2 = simd::sub(V3,T1);
			MemoryOp<Element,simd>::store(EFGH,V2);
		}

		inline void Butterfly_DIT_mod4p_firststeps (Element* ABCD, Element* EFGH,
													const vect_t& W,
													const vect_t& Wp,
													const vect_t& P, const vect_t& P2) {
			// First 2 steps
			// First step
			vect_t V1,V2,V3,V4,T1,T2,T3; // T4
			// T1=[A B C D], T2=[E F G H]
			T1 = MemoryOp<Element,simd>::load(ABCD);
			T2 = MemoryOp<Element,simd>::load(EFGH);

			// V1=[AECG], V2=[BFDH]
			MemoryOp<Element,simd>::unpacklohi_twice4(V1,V2,T1,T2);

			// V3 = V1 + V2
			// Rk: No need for (. mod 2P) since entries are <P
			V3 = simd::add(V1,V2);
			// V4 = V1 + (P - V2)
			// Rk: No need for (. mod 2P) since entries are <P
			T1 = simd::sub(V2,P);
			V4 = simd::sub(V1,T1);

			//TODO: Useless?
			//MemoryOp<Element,simd>::unpacklohi4(V1,V2,V3,V4);

			// Second step
			// T1 = [D D H H]
			T1 = MemoryOp<Element,simd>::unpackhi4(V4,V4);

			T2 = mul_mod_half<simd, SimdComp>(T1, W, P, Wp);

			T2 = simd::template shuffle<0xDD>(T2);
			//T2 = simd::template shuffle<0xDD>(T2);

			//At this point, T2 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]

			// At this time I have V3=[A E C G], V4=[B F ? ?], T2=[? ? D H]
			// I need V1 = [A B E F], V2 = [C D G H]
			// This is not refactored in MemoryOp::... because of different arguments (V3,V4) and (V3,T2)
			V1 = MemoryOp<Element,simd>::unpacklo4(V3,V4);
			V2 = MemoryOp<Element,simd>::unpackhi4(V3,T2);

			// T1 = V1 + V2
			T1 = simd::add(V1,V2);
			// T2 = V1 - (V2 - 2P)
			T3 = simd::sub(V2,P2);
			T2 = simd::sub(V1,T3);

			MemoryOp<Element,simd>::unpacklohi2(V1,V2,T1,T2);

			// Store
			MemoryOp<Element,simd>::store(ABCD,V1);
			MemoryOp<Element,simd>::store(EFGH,V2);
		}

		inline void Butterfly_DIF_mod2p (Element* ABCD, Element* EFGH,
										 const Element* alpha, const Element* alphap,
										 const vect_t& P, const vect_t& P2) {
			vect_t V1,V2,V3,V4,W,Wp,T;
			// V1=[A B C D], V2=[E F G H]
			V1 = MemoryOp<Element,simd>::load(ABCD);
			V2 = MemoryOp<Element,simd>::load(EFGH);
			W  = MemoryOp<Element,simd>::load(alpha);
			Wp = MemoryOp<Element,simd>::load(alphap);
			// V3 = V1 + V2 mod
			V3 = add_mod<simd >(V1,V2,P2);
			MemoryOp<Element,simd>::store(ABCD,V3);
			// V4 = (V1+(2P-V2))alpha mod 2P
			T = simd::sub(V2,P2);
			V4 = simd::sub(V1,T);

			T = mul_mod<simd >(V4,W,P,Wp);// T is the result
			MemoryOp<Element,simd>::store(EFGH,T);
		}

		inline void Butterfly_DIF_mod2p_laststeps(Element* ABCD, Element* EFGH,
												  const vect_t& W,
												  const vect_t& Wp,
												  const vect_t& P, const vect_t& P2) {
			vect_t V1,V2,V3,V4,V5,V6; // V7
			// V1=[A B C D], V2=[E F G H]
			V1 = MemoryOp<Element,simd>::load(ABCD);
			V2 = MemoryOp<Element,simd>::load(EFGH);

			/* 1st step */
			// V3=[A E B F], V4=[C G D H]
			MemoryOp<Element,simd>::unpacklohi4(V3,V4,V1,V2);

			// V1 = V3 + V4 mod 2P
			// P2 = [2p 2p 2p 2p]
			V1 = add_mod<simd >(V3,V4,P2);
			// V2 = (V3+(2P-V4))alpha mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V2 = reduce<simd >(V6, P2);
			// V4 = [D D H H]
			V4 = MemoryOp<Element,simd>::unpackhi4(V2,V2);

			// V3 = [* D * H]
			V3 = mul_mod_half<simd, SimdComp>(V4, W, P, Wp);

			//At this point, V3 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]
			V3 = simd::template shuffle<0xDD>(V3);
			//V3 = simd::template shuffle<0xDD>(V3); // 0xDD = [3 1 3 1]_base4

			// At this time I have V1=[A E B F], V2=[C G ? ?], V3=[? ? D H]
			// I need V3 = [A C E G], V4 = [B D F H]
			// This is not refactored in MemoryOp::... because of different arguments (V1,V3) and (V1,V2)
			V4 = MemoryOp<Element,simd>::unpackhi4(V1,V3);
			V3 = MemoryOp<Element,simd>::unpacklo4(V1,V2);

			/* 2nd step */
			// V1 = V3 + V4 mod 2P
			V1 = add_mod<simd >(V3,V4,P2);
			// V2 = V3 + (2P - V4) mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V2 = reduce<simd >(V6, P2);
			// Result in V1 = [A C E G]  and V2 = [B D F H]
			// Transform to V3=[A B C D], V4=[E F G H]
			MemoryOp<Element,simd>::unpacklohi4(V3,V4,V1,V2);
			// Store
			MemoryOp<Element,simd>::store(ABCD,V3);
			MemoryOp<Element,simd>::store(EFGH,V4);
		}

	}; // FFT_butterflies<Field, 4>


	template<typename Field, typename simd>
	class FFT_butterflies<Field, simd, 8, IS_INTEGRAL> : public FFT_init<Field> {
	public:

		using Element = typename Field::Element;
		using vect_t = typename simd::vect_t;
		using SimdComp = typename SimdCompute_t<simd,Field>::Compute_t;

		FFT_butterflies(const FFT_init<Field>& f_i) : FFT_init<Field>(f_i) {
			linbox_check(simd::vect_size == 8);
		}

		// TODO include P, P2 in precomp
		inline void Butterfly_DIT_mod4p (Element* ABCDEFGH, Element* IJKLMNOP,
										 const Element* alpha, const Element* alphap,
										 const vect_t& P, const vect_t& P2) {
			vect_t V1,V2,V3,V4,W,Wp,T1;
			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = MemoryOp<Element,simd>::load(ABCDEFGH);
			V2 = MemoryOp<Element,simd>::load(IJKLMNOP);
			W  = MemoryOp<Element,simd>::load(alpha);
			Wp = MemoryOp<Element,simd>::load(alphap);

			// V3 = V1 mod 2P
			V3 = reduce<simd>(V1, P2);

			// V4 = V2 * W mod P
			V4 = mul_mod<simd>(V2,W,P,Wp);

			// V1 = V3 + V4
			V1 = simd::add(V3,V4);
			MemoryOp<Element,simd>::store(ABCDEFGH,V1);

			// V2 = V3 - (V4 - 2P)
			T1 = simd::sub(V4,P2);
			V2 = simd::sub(V3,T1);
			MemoryOp<Element,simd>::store(IJKLMNOP,V2);
		}

		inline void Butterfly_DIT_mod4p_firststeps (Element* ABCDEFGH, Element* IJKLMNOP,
													const vect_t& alpha,const vect_t& alphap,
													const vect_t& beta ,const vect_t& betap,
													const vect_t& P    ,const vect_t& P2) {
			// First 3 steps
			vect_t V1,V2,V3,V4,V5,V6,V7,Q;
			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = MemoryOp<Element,simd>::load(ABCDEFGH);
			V2 = MemoryOp<Element,simd>::load(IJKLMNOP);

			/*********************************************/
			/* 1st STEP */
			/*********************************************/
			// Transform to V3=[A I C K E M G O], V4=[B J D L F N H P]
			MemoryOp<Element,simd>::unpacklohi_twice8(V6,V7,V1,V2);
			MemoryOp<Element,simd>::unpacklohi_twice4(V3,V4,V6,V7);

			// V1 = V3 + V4;       V1 = [A I C K E M G O]
			// Rk: No need for (. mod 2P) since entries are <P
			V1 = simd::add(V3,V4);

			// V2 = V3 + (P - V4); V2 = [B J D L F N H P]
			// Rk: No need for (. mod 2P) since entries are <P
			V6 = simd::sub(V4,P);
			V2 = simd::sub(V3,V6);

			/*********************************************/
			/* 2nd STEP */
			/*********************************************/
			// V5 = [D D L L H H P P]
			V5 = MemoryOp<Element,simd>::unpackhi_twice8(V2,V2);

			// V3 = [* D * L * H * P]
			V3 = mul_mod_half<simd,SimdComp>(V5,alpha,P,alphap);

			// V7 = [D L D L H P H P]
			V7 = MemoryOp<Element,simd>::shuffletwice8_DD(V3);
			//V7 = simd::template shuffle_twice<0xDD>(V3); // 0xDD = 221 = [3 1 3 1]_base4

			// V3= [A B I J E F M N], V4=[C D K L G H O P]
			V3 = MemoryOp<Element,simd>::unpacklo_twice8(V1,V2);
			V4 = MemoryOp<Element,simd>::unpackhi_twice8(V1,V7);

			// V1 = V3+V4
			V1 = simd::add(V3,V4);
			// V2 = V3 - (V4 - 2P)
			V7 = simd::sub(V4,P2);
			V2 = simd::sub(V3,V7);

			/*********************************************/
			/* 3nd STEP */
			/*********************************************/
			// V3= [A B C D I J K L] V4= [E F G H M N O P]
			MemoryOp<Element,simd>::unpacklohi_twice4(V6,V7,V1,V2);
			MemoryOp<Element,simd>::unpacklohi2(V3,V4,V6,V7);

			// V6= V3 mod 2P
			V6 = reduce<simd >(V3, P2);

			// V7= V4.beta mod p
			V7 = mul_mod<simd >(V4,beta,P,betap);

			// V1 = V6+V7
			V1 = simd::add(V6,V7);

			// V2 = V6 - (V7 - 2P)
			V5 = simd::sub(V7,P2);
			V2 = simd::sub(V6,V5);

			/*********************************************/
			// V3=[A B C D E F G H] V4=[I J K L M N O P]
			MemoryOp<Element,simd>::unpacklohi2(V3,V4,V1,V2);

			// Store
			MemoryOp<Element,simd>::store(ABCDEFGH,V3);
			MemoryOp<Element,simd>::store(IJKLMNOP,V4);
		}

		inline void Butterfly_DIF_mod2p (Element* ABCDEFGH, Element* IJKLMNOP,
										 const Element* alpha, const Element* alphap,
										 const vect_t& P, const vect_t& P2) {
			vect_t V1,V2,V3,V4,W,Wp,T;
			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = MemoryOp<Element,simd>::load(ABCDEFGH);
			V2 = MemoryOp<Element,simd>::load(IJKLMNOP);
			W  = MemoryOp<Element,simd>::load(alpha);
			Wp = MemoryOp<Element,simd>::load(alphap);

			// V3 = V1 + V2 mod

			V3 = add_mod<simd >(V1,V2,P2);

			MemoryOp<Element,simd>::store(ABCDEFGH,V3);

			// V4 = (V1+(2P-V2))alpha mod 2P
			T = simd::sub(V2,P2);
			V4 = simd::sub(V1,T);
			T = mul_mod<simd >(V4,W,P,Wp);// T is the result
			MemoryOp<Element,simd>::store(IJKLMNOP,T);
		}

		inline void Butterfly_DIF_mod2p_laststeps(Element* ABCDEFGH, Element* IJKLMNOP,
												  const vect_t& alpha,const vect_t& alphap,
												  const vect_t& beta ,const vect_t& betap,
												  const vect_t& P, const vect_t& P2) {
			// Last 3 steps
			vect_t V1,V2,V3,V4,V5,V6,V7,Q;

			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = MemoryOp<Element,simd>::load(ABCDEFGH);
			V2 = MemoryOp<Element,simd>::load(IJKLMNOP);

			/* 1st step */
			// V3=[A B C D I J K L] V4=[E F G H M N O P]
			MemoryOp<Element,simd>::unpacklohi2(V3,V4,V1,V2);

			// V1 = V3 + V4 mod 2P
			// P2 = [2p 2p 2p 2p]
			V1 = add_mod<simd >(V3,V4,P2);

			// V2 = (V3+(2P-V4))alpha mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V7 = reduce<simd >(V6, P2);
			V2 = mul_mod<simd >(V7,alpha,P,alphap);

			/* 2nd step */

			// V3=[A E B F I M J N] V4=[C G D H K O L P]
			MemoryOp<Element,simd>::unpacklohi_twice8(V3,V4,V1,V2);

			// V1 = V3 + V4 mod 2P
			// P2 = [2p 2p 2p 2p]
			V1 = add_mod<simd >(V3,V4,P2);

			// V2 = (V3+(2P-V4))alpha mod 2P
			// V7 =  (V3+(2P-V4)) mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V7 = reduce<simd >(V6, P2);

			// V4 = [D D H H L L P P ]
			V4 = MemoryOp<Element,simd>::unpackhi_twice8(V7,V7);

			// V3 = [ * D * H * L * P]
			V3 = mul_mod_half<simd,SimdComp>(V4,beta,P,betap);

			// V2=[D H D H L P L P] but only [* * D H * * L P] matters
			V2 = MemoryOp<Element,simd>::shuffletwice8_DD(V3);
			//V2 = simd::template shuffle_twice<0xDD>(V3); // 0xDD = 221 = [3 1 3 1]_base4

			/* 3rd step */
			// At this time I have V1=[A B E F I J M N], V7=[C G * * K O * *], V2=[* * D H * * L P]
			// I need V3 = [A C E G I K M O], V4=[B D F H J L N P]
			V3 = MemoryOp<Element,simd>::unpacklo_twice8(V1,V7);
			V4 = MemoryOp<Element,simd>::unpackhi_twice8(V1,V2);

			// V1 = V3 + V4 mod 2P
			V1 = add_mod<simd >(V3,V4,P2);

			// V2 = V3 + (2P - V4) mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V2 = reduce<simd >(V6, P2);

			// Result in    V1=[A C E G I K M O] V2=[B D F H J L N P]
			// Transform to V3=[A B C D I J K L],V4=[E F G H M N O P]
			MemoryOp<Element,simd>::unpacklohi_twice8(V3,V4,V1,V2);

			// Transform to V1=[A B C D E F G H], V2=[I J K L M N O P]
			MemoryOp<Element,simd>::unpacklohi2(V1,V2,V3,V4);

			// Store
			MemoryOp<Element,simd>::store(ABCDEFGH,V1);
			MemoryOp<Element,simd>::store(IJKLMNOP,V2);


		}


	}; // FFT_butterflies<Field, 8>

    template<typename Field>
    class FFT_butterflies<Field, NoSimd<typename Field::Element>, 1, IS_FLOATING> : public FFT_init<Field> {
    public:

        using Element = typename Field::Element;

        FFT_butterflies(const FFT_init<Field>& f_i) : FFT_init<Field>(f_i) {}

        inline void Butterfly_DIF (Element& A, Element& B, const Element& alpha)
        {
            Element tmp;
            this->fld->assign(tmp, A);
            this->fld->addin(A, B);
            this->fld->sub(B, tmp, B);
            this->fld->mulin(B, alpha);
        }

        inline void Butterfly_DIT (Element& A, Element& B, const Element& alpha)
        {
            Element tmp;
            this->fld->mul(tmp, alpha, B);
            this->fld->sub(B, A, tmp);
            this->fld->addin(A, tmp);
        }

    }; // FFT_butterflies<Field, 1, IS_FLOATING>

	// ATTENTION à tous les uint64_t, SimdComp restants !!!!

	template<typename Field, typename simd>
	class FFT_butterflies<Field, simd, 4, IS_FLOATING> : public FFT_init<Field> {
	public:

		using Element = typename Field::Element;
		using vect_t = typename simd::vect_t;
		using SimdComp = typename SimdCompute_t<simd,Field>::Compute_t;

		FFT_butterflies(const FFT_init<Field>& f_i) : FFT_init<Field>(f_i) {
			linbox_check(simd::vect_size == 4);
		}

		// TODO include P, P2 in precomp
		// TODO : Same functions Butterfly_DIT_mod4p Butterfly_DIF_mod2p in FFT_butterflies<Field, 8>
		inline void Butterfly_DIT (Element* ABCD, Element* EFGH, const Element* alpha, const vect_t& P, const vect_t& U) {
			vect_t V1,V2,W,T1;

			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = MemoryOp<Element,simd>::load(ABCD);
			V2 = MemoryOp<Element,simd>::load(EFGH);
			W  = MemoryOp<Element,simd>::load(alpha);

			//Element tmp;
            T1 = mul_mod<simd>(W, V2, P, U);
            V2 = sub_mod<simd>(V1, T1, P);
            V1 = add_mod<simd>(V1, T1, P);

            // this->fld->mul(tmp, alpha, B);
            // this->fld->sub(B, A, tmp);
            // this->fld->addin(A, tmp);

			MemoryOp<Element,simd>::store(ABCD,V1);
			MemoryOp<Element,simd>::store(EFGH,V2);
		}

		inline void Butterfly_DIT_firststeps (Element* ABCD, Element* EFGH, const vect_t& W, const vect_t& P, const vect_t& U) {
			// First 2 steps
			// First step
			vect_t V1,V2,V3,V4,T1,T2;
			// T1=[A B C D], T2=[E F G H]
			T1 = MemoryOp<Element,simd>::load(ABCD);
			T2 = MemoryOp<Element,simd>::load(EFGH);

			// V1=[AECG], V2=[BFDH]
			MemoryOp<Element,simd>::unpacklohi_twice4(V1,V2,T1,T2);

			V3 = add_mod<simd>(V1,V2,P);
			V4 = sub_mod<simd>(V1,V2,P);

			//TODO: Useless?
			// V1 = [ABEF], V2= [CDGH]
			//MemoryOp<Element,simd>::unpacklohi4(V1,V2,V3,V4);

			// Second step
			//TODO : Instead do [CDGH] * [1 a 1 a] ?
			// T1 = [D D H H]
			T1 = MemoryOp<Element,simd>::unpackhi4(V4,V4);

			T2 = mul_mod<simd>(T1, W, P, U);

			T2 = simd::template shuffle<0xDD>(T2);

			//At this point, T2 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]

			// At this time I have V3=[A E C G], V4=[B F ? ?], T2=[? ? D H]
			// V1 = [A B E F], V2 = [C D G H]
			V1 = MemoryOp<Element,simd>::unpacklo4(V3,V4);
			V2 = MemoryOp<Element,simd>::unpackhi4(V3,T2);

			T1 = add_mod<simd>(V1,V2,P);
			T2 = sub_mod<simd>(V1,V2,P);

			MemoryOp<Element,simd>::unpacklohi2(V1,V2,T1,T2);

			// Store
			MemoryOp<Element,simd>::store(ABCD,V1);
			MemoryOp<Element,simd>::store(EFGH,V2);
		}

		inline void Butterfly_DIF (Element* ABCD, Element* EFGH, const Element* alpha, const vect_t& P, const vect_t& U) {

			vect_t V1,V2,W,T;
			// V1=[A B C D], V2=[E F G H]
			V1 = MemoryOp<Element,simd>::load(ABCD);
			V2 = MemoryOp<Element,simd>::load(EFGH);
			W  = MemoryOp<Element,simd>::load(alpha);

			T = V1;
			V1 = add_mod<simd>(V1, V2, P);
			V2 = sub_mod<simd>(T, V2, P);
			V2 = mul_mod<simd>(V2, W, P, U);

            // Element tmp;
            // this->fld->assign(tmp, A);
            // this->fld->addin(A, B);
            // this->fld->sub(B, tmp, B);
            // this->fld->mulin(B, alpha);

			MemoryOp<Element,simd>::store(ABCD,V1);
			MemoryOp<Element,simd>::store(EFGH,V2);

		}

		inline void Butterfly_DIF_laststeps(Element* ABCD, Element* EFGH, const vect_t& W, const vect_t& P, const vect_t& U) {
			vect_t V1,V2,V3,V4;
			// V1=[A B C D], V2=[E F G H]
			V1 = MemoryOp<Element,simd>::load(ABCD);
			V2 = MemoryOp<Element,simd>::load(EFGH);

			/* 1st step */
			// V3=[A E B F], V4=[C G D H]
			MemoryOp<Element,simd>::unpacklohi4(V3,V4,V1,V2);

			V1 = add_mod<simd>(V3,V4,P);
			V2 = sub_mod<simd>(V3,V4,P);

			// V4 = [D D H H]
			V4 = MemoryOp<Element,simd>::unpackhi4(V2,V2);

			V3 = mul_mod<simd>(V4, W, P, U);

			//At this point, V3 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]
			V3 = simd::template shuffle<0xDD>(V3);
			//V3 = simd::template shuffle<0xDD>(V3); // 0xDD = [3 1 3 1]_base4

			// At this time I have V1=[A E B F], V2=[C G ? ?], V3=[? ? D H]
			// I need V3 = [A C E G], V4 = [B D F H]
			V4 = MemoryOp<Element,simd>::unpackhi4(V1,V3);
			V3 = MemoryOp<Element,simd>::unpacklo4(V1,V2);

			/* 2nd step */
			V1 = add_mod<simd>(V3,V4,P);
			V2 = sub_mod<simd>(V3,V4,P);

			// Result in V1 = [A C E G]  and V2 = [B D F H]
			// Transform to V3=[A B C D], V4=[E F G H]
			MemoryOp<Element,simd>::unpacklohi4(V3,V4,V1,V2);
			// Store
			MemoryOp<Element,simd>::store(ABCD,V3);
			MemoryOp<Element,simd>::store(EFGH,V4);
		}

	}; // FFT_butterflies<Field, 4>

}

#undef IS_INTEGRAL
#undef IS_FLOATING

#endif // __LINBOX_polynomial_fft_butterflies_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
