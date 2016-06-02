/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#ifndef additional_modular_simd_functions
#define additional_modular_simd_functions

#define Simd_vect typename Simd::vect_t

template <class Simd>
inline Simd_vect reduce (const Simd_vect& a, const Simd_vect& p) {
	Simd_vect t = Simd::greater(p,a);
	return Simd::sub(a, Simd::vandnot(p,t));
}

template <class Simd>
inline Simd_vect add_mod (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p) {
	Simd_vect c = Simd::add(a,b);
	return reduce<Simd>(c, p);
}

template <class Simd>
inline Simd_vect mul_mod (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p, const Simd_vect& bp) {
	Simd_vect q = Simd::mulhi(a,bp);
	Simd_vect c = Simd::mullo(a,b);
	Simd_vect t = Simd::mullo(q,p);
	return Simd::sub(c,t);
}

/*
 * a = [a0, a0, a2, a2, ...]
 * b = [?, b0, ?, b2, ...] with bp its shoup mul_mod precomputation [b0p ? b2p ?, ... ]
 * Return [?, (a0*b0) mod p, ?, (a2*b2) mod p, ... ]
 */
template <class Simd, class SimdCompute_t>
inline Simd_vect mul_mod_half (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p, const Simd_vect& bp) {
#if 1
	return mul_mod<Simd>(a, b , p, bp);
#else
	// T2 = a * bp mod 2^64 (for Modular<Element = uint32, Compute_t = uint64>)
	// bp = [b0p ? b2p ?, ... ] is enough
	Simd_vect T2 = SimdCompute_t::mulx(a,bp);
	Simd_vect T3 = Simd::mullo(T2,p);
	// At this point T3= [? quo(D)*p ? quo(H)*p] mod 2^32
	// T4 = [D D H H] * [?, b0, ?, b2] mod 2^32
	T2 = Simd::mullo(a,b);
	return Simd::sub(T2,T3);
#endif
}


/* First swizzle operation
 * [ABCD], [EFGH] => [AECG], [BFDH]
*/
template<class T, class Simd = Simd<T>>
struct swizzle_1 {
	static inline void swizzle (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
};

template<class T>
struct swizzle_1<T, Simd128<T>> {
	using simd = Simd128<T>;
	using simd_vect = typename simd::vect_t;

	static inline void swizzle (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
		simd_vect a1 = simd::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
		simd_vect b1 = simd::template shuffle<0xD8>(b);
		s1 = simd::unpacklo(a1,b1);
		s2 = simd::unpackhi(a1,b1);
	}
};

template<class T>
struct swizzle_1<T, Simd256<T>> {
	using simd = Simd256<T>;
	using simd_vect = typename simd::vect_t;

	static inline void swizzle (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
		s1 = simd::unpacklo_twice(a,b);
		s2 = simd::unpackhi_twice(a,b);
	}
};

/* Second swizzle operation
 * [ABCD], [EFGH] => [AEBF], [CGDH]
*/
template<class T, class Simd = Simd<T>>
struct swizzle_2 {
	// The implementation is already generic (unpacklo and _hi are generic) but not optimized
	static inline void swizzle (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
};

template<class T>
struct swizzle_2<T, Simd128<T>> {
	using simd = Simd128<T>;
	using simd_vect = typename simd::vect_t;

	static inline void swizzle (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
		s1 = simd::unpacklo(a,b);
		s2 = simd::unpackhi(a,b);
	}
};

template<class T>
struct swizzle_2<T, Simd256<T>> {
	using simd = Simd256<T>;
	using simd_vect = typename simd::vect_t;

	static inline void swizzle (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
		// Old simd::unpacklohi(s1,s2,a,b);
		simd_vect a1 = simd::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4 so a -> [a0,a2,a1,a3]
		simd_vect b1 = simd::template shuffle<0xD8>(b); // 0xD8 = 3120 base_4
		s1 = simd::unpacklo_twice(a1, b1);
		s2 = simd::unpackhi_twice(a1, b1);
	}
};

/* Third swizzle operation
 * [A B E F], [C D G H] => [ABCD], [EFGH]
*/
template<class T, class Simd = Simd<T>>
struct swizzle_3 {
	static inline void swizzle (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
};

template<class T>
struct swizzle_3<T, Simd128<T>> {
	using simd = Simd128<T>;
	using simd_vect = typename simd::vect_t;

	static inline void swizzle (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
		using Simd128_64 = Simd128<uint64_t>;
		s1 = Simd128_64::unpacklo(a,b);
		s2 = Simd128_64::unpackhi(a,b);
	}
};

template<class T>
struct swizzle_3<T, Simd256<T>> {
	using simd = Simd256<T>;
	using simd_vect = typename simd::vect_t;

	static inline void swizzle (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
		s1 = simd::unpacklo128(a,b);
		s2 = simd::unpackhi128(a,b);
	}
};

#undef Simd_vect
#endif

template <typename simd, typename Field>
struct SimdCompute_t {};

template <typename Field>
struct SimdCompute_t<Simd128<typename Field::Element>, Field> {
	using Compute_t = Simd128<typename Field::Compute_t>;
};

template <typename Field>
struct SimdCompute_t<Simd256<typename Field::Element>, Field> {
	using Compute_t = Simd256<typename Field::Compute_t>;
};

namespace LinBox {

	// TODO : template by the number of steps

	template<typename Field, typename simd = Simd<typename Field::Element>, uint8_t byn = simd::vect_size>
	class FFT_butterflies : public FFT_init<Field> {
	public:
		FFT_butterflies(const FFT_init<Field>& f_i) : FFT_init<Field>(f_i) {
			std::cerr<<"Not implemented !\n";
		}
	}; // FFT_butterflies

	template<typename Field>
	class FFT_butterflies<Field, NoSimd<typename Field::Element>, 1> : public FFT_init<Field> {
	public:

		using Element = typename Field::Element;

		FFT_butterflies(const FFT_init<Field>& f_i) : FFT_init<Field>(f_i) {}

		inline void Butterfly_DIT_mod4p (Element& A, Element& B, const Element& alpha, const Element& alphap) {
			// Harvey's algorithm
			// 0 <= A,B < 4*p, p < 2^32 / 4
			// alphap = Floor(alpha * 2^ 32 / p])

			// TODO : replace by substract if greater
			if (A >= this->_dpl) A -= this->_dpl;

			// TODO : replace by mul_mod_shoup
			Element tmp = ((uint32_t) alphap * (uint64_t)B) >> 32;
			tmp = (uint64_t)alpha * B - tmp * this->_pl;

			// TODO : replace by add_r and sub_r
			B = A + (this->_dpl - tmp);
			//        B &= 0XFFFFFFFF;
			A += tmp;
		}

		inline void Butterfly_DIF_mod2p (Element& A, Element& B, const Element& alpha, const Element& alphap) {
			//std::cout<<A<<" $$ "<<B<<"("<<alpha<<","<<alphap<<" ) -> ";
			// Harvey's algorithm
			// 0 <= A,B < 2*p, p < 2^32 / 4
			// alphap = Floor(alpha * 2^ 32 / p])

			Element tmp = A;

			A += B;

			if (A >= this->_dpl) A -= this->_dpl;

			B = tmp + (this->_dpl - B);

			tmp = ((uint32_t) alphap * (uint64_t)B) >> 32;
			B = (uint64_t)alpha * B - tmp * this->_pl;
			//B &= 0xFFFFFFFF;
			//std::cout<<A<<" $$ "<<B<<"\n ";
		}

	}; // FFT_butterflies<Field, 1>

	// ATTENTION à tous les uint64_t, SimdComp restants !!!!

	template<typename Field, typename simd>
	class FFT_butterflies<Field, simd, 4> : public FFT_init<Field> {
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
			V1 = simd::load(ABCD);
			V2 = simd::load(EFGH);
			W  = simd::load(alpha);
			Wp = simd::load(alphap);

			// V3 = V1 mod 2P
			V3 = reduce<simd>(V1, P2);

			// V4 = V2 * W mod P
			V4 = mul_mod<simd>(V2,W,P,Wp);

			// V1 = V3 + V4
			V1 = simd::add(V3,V4);
			simd::storeu(ABCD,V1);

			// V2 = V3 - (V4 - 2P)
			T1 = simd::sub(V4,P2);
			V2 = simd::sub(V3,T1);
			simd::storeu(EFGH,V2);
		}

		inline void Butterfly_DIT_mod4p_firststeps (Element* ABCD, Element* EFGH,
													const vect_t& W,
													const vect_t& Wp,
													const vect_t& P, const vect_t& P2) {
			// First 2 steps
			// First step
			vect_t V1,V2,V3,V4,T1,T2,T3,T4;
			// V1=[A B C D], V2=[E F G H]
			V1 = simd::load(ABCD);
			V2 = simd::load(EFGH);

			// V1=[AECG], V2=[BFDH]
			swizzle_1<Element, simd>::swizzle(V1,V2,V1,V2);

			// V3 = V1 + V2
			// Rk: No need for (. mod 2P) since entries are <P
			V3 = simd::add(V1,V2);
			// V4 = V1 + (P - V2)
			// Rk: No need for (. mod 2P) since entries are <P
			T1 = simd::sub(V2,P);
			V4 = simd::sub(V1,T1);


			// Second step
			// T1 = [D D H H]
			T1 = simd::unpackhi(V4,V4);

			T2 = mul_mod_half<simd, SimdComp>(T1, W, P, Wp);

			T2 = simd::template shuffle<0xDD>(T2);
			//At this point, T2 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]

			// At this time I have V3=[A E C G], V4=[B F ? ?], T2=[? ? D H]
			// I need V1 = [A B E F], V2 = [C D G H]
			V1 = simd::unpacklo(V3,V4);
			V2 = simd::unpackhi(V3,T2);
			// T1 = V1 + V2
			T1 = simd::add(V1,V2);
			// T2 = V1 - (V2 - 2P)
			T3 = simd::sub(V2,P2);
			T2 = simd::sub(V1,T3);

			swizzle_3<Element, simd>::swizzle(V1,V2,T1,T2);

			// Store
			simd::store(ABCD,V1);
			simd::store(EFGH,V2);
		}

		inline void Butterfly_DIF_mod2p (Element* ABCD, Element* EFGH,
										 const Element* alpha, const Element* alphap,
										 const vect_t& P, const vect_t& P2) {
			vect_t V1,V2,V3,V4,W,Wp,T;
			// V1=[A B C D], V2=[E F G H]
			V1 = simd::load(ABCD);
			V2 = simd::load(EFGH);
			W  = simd::load(alpha);
			Wp = simd::load(alphap);
			// V3 = V1 + V2 mod
			V3 = add_mod<simd >(V1,V2,P2);
			simd::storeu(ABCD,V3);
			// V4 = (V1+(2P-V2))alpha mod 2P
			T = simd::sub(V2,P2);
			V4 = simd::sub(V1,T);
			T = mul_mod<simd >(V4,W,P,Wp);// T is the result
			simd::storeu(EFGH,T);
		}

		inline void Butterfly_DIF_mod2p_laststeps(Element* ABCD, Element* EFGH,
												  const vect_t& W,
												  const vect_t& Wp,
												  const vect_t& P, const vect_t& P2) {
			vect_t V1,V2,V3,V4,V5,V6,V7;
			// V1=[A B C D], V2=[E F G H]
			V1 = simd::load(ABCD);
			V2 = simd::load(EFGH);

			swizzle_2<Element,simd>::swizzle(V3,V4,V1,V2);
//			// V3=[A E B F], V4=[C G D H]
//			V3 = simd::unpacklo(V1,V2);
//			V4 = simd::unpackhi(V1,V2);

			// V1 = V3 + V4 mod 2P
			// P2 = [2p 2p 2p 2p]
			V1 = add_mod<simd >(V3,V4,P2);
			// V2 = (V3+(2P-V4))alpha mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V2 = reduce<simd >(V6, P2);
			// V4 = [D D H H]
			V4 = simd::unpackhi(V2,V2);

			V3 = mul_mod_half<simd, SimdComp>(V4, W, P, Wp);
			V3 = simd::template shuffle<0xDD>(V3);
			//At this point, V2 = [D*Wmodp H*Wmodp D*Wmodp H*Wmodp]

			// At this time I have V1=[A E B F], V2=[C G ? ?], V3=[? ? D H]
			// I need V3 = [A C E G], V4 = [B D F H]
			V4 = simd::unpackhi(V1,V3);
			V3 = simd::unpacklo(V1,V2);

			// V1 = V3 + V4 mod 2P
			V1 = add_mod<simd >(V3,V4,P2);
			// V2 = V3 + (2P - V4) mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V2 = reduce<simd >(V6, P2);
			// Result in V1 = [A C E G]  and V2 = [B D F H]
			// Transform to V3=[A B C D], V4=[E F G H]
			V3 = simd::unpacklo(V1,V2);
			V4 = simd::unpackhi(V1,V2);
			// Store
			simd::store(ABCD,V3);
			simd::store(EFGH,V4);
		}

	}; // FFT_butterflies<Field, 4>


	template<typename Field, typename simd>
	class FFT_butterflies<Field, simd, 8> : public FFT_init<Field> {
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
			V1 = simd::load(ABCDEFGH);
			V2 = simd::load(IJKLMNOP);
			W  = simd::load(alpha);
			Wp = simd::load(alphap);

			// V3 = V1 mod 2P
			V3 = reduce<simd>(V1, P2);

			// V4 = V2 * W mod P
			V4 = mul_mod<simd>(V2,W,P,Wp);

			// V1 = V3 + V4
			V1 = simd::add(V3,V4);
			simd::storeu(ABCDEFGH,V1);

			// V2 = V3 - (V4 - 2P)
			T1 = simd::sub(V4,P2);
			V2 = simd::sub(V3,T1);
			simd::storeu(IJKLMNOP,V2);
		}

		inline void Butterfly_DIT_mod4p_firststeps (Element* ABCDEFGH, Element* IJKLMNOP,
													const vect_t& alpha,const vect_t& alphap,
													const vect_t& beta ,const vect_t& betap,
													const vect_t& P    ,const vect_t& P2) {
			// First 3 steps
			vect_t V1,V2,V3,V4,V5,V6,V7,Q;
			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = simd::load(ABCDEFGH);
			V2 = simd::load(IJKLMNOP);

			/*********************************************/
			/* 1st STEP */
			/*********************************************/
			// Transform to V3=[A I C K E M G O], V4=[B J D L F N H P]
			V6 = simd::unpacklo_twice(V1,V2); // V6=[A I B J E M F N]
			V7 = simd::unpackhi_twice(V1,V2); // V7=[C K D L G O H P]
			V3 = SimdComp::unpacklo_twice(V6,V7); // V3=[A I C K E M G O]
			V4 = SimdComp::unpackhi_twice(V6,V7); // V4=[B J D L F N H P]

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
			V5 = simd::unpackhi_twice(V2,V2);

			V3 = mul_mod_half<simd,SimdComp>(V5,alpha,P,alphap);

			// V7=[D L * * H P * *]
			V7 = simd::template shuffle_twice<0xDD>(V3); // 0xDD = 221

			// V6 = [B J D L F N H P]
			V6 = SimdComp::unpacklo_twice(V2,V7);
			// V3= [A B I J E F M N], V4=[C D K L G H O P]
			V3 = simd::unpacklo_twice(V1,V6);
			V4 = simd::unpackhi_twice(V1,V6);

			// V1 = V3+V4
			V1 = simd::add(V3,V4);
			// V2 = V3 - (V4 - 2P)
			V7 = simd::sub(V4,P2);
			V2 = simd::sub(V3,V7);

			/*********************************************/
			/* 3nd STEP */
			/*********************************************/
			// V3= [A B C D I J K L] V4= [E F G H M N O P]
			V6 = SimdComp::unpacklo_twice(V1,V2);
			V7 = SimdComp::unpackhi_twice(V1,V2);
			V3 = SimdComp::unpacklo128(V6,V7);
			V4 = SimdComp::unpackhi128(V6,V7);

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
			V3 = SimdComp::unpacklo128(V1,V2);
			V4 = SimdComp::unpackhi128(V1,V2);

			// Store
			simd::storeu(ABCDEFGH,V3);
			simd::storeu(IJKLMNOP,V4);
		}

		inline void Butterfly_DIF_mod2p (Element* ABCDEFGH, Element* IJKLMNOP,
										 const Element* alpha, const Element* alphap,
										 const vect_t& P, const vect_t& P2) {
			vect_t V1,V2,V3,V4,W,Wp,T;
			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = simd::load(ABCDEFGH);
			V2 = simd::load(IJKLMNOP);
			W  = simd::load(alpha);
			Wp = simd::load(alphap);

			// V3 = V1 + V2 mod

			V3 = add_mod<simd >(V1,V2,P2);

			simd::storeu(ABCDEFGH,V3);

			// V4 = (V1+(2P-V2))alpha mod 2P
			T = simd::sub(V2,P2);
			V4 = simd::sub(V1,T);
			T = mul_mod<simd >(V4,W,P,Wp);// T is the result
			simd::storeu(IJKLMNOP,T);
		}

		inline void Butterfly_DIF_mod2p_laststeps(Element* ABCDEFGH, Element* IJKLMNOP,
												  const vect_t& alpha,const vect_t& alphap,
												  const vect_t& beta ,const vect_t& betap,
												  const vect_t& P, const vect_t& P2) {
			// Last 3 steps
			vect_t V1,V2,V3,V4,V5,V6,V7,Q;

			// V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = simd::load(ABCDEFGH);
			V2 = simd::load(IJKLMNOP);

			/* 1st step */
			// V3=[A B C D I J K L] V4=[E F G H M N O P]
			V3 = SimdComp::unpacklo128(V1,V2);
			V4 = SimdComp::unpackhi128(V1,V2);

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
			V3 = simd::unpacklo_twice(V1,V2);
			V4 = simd::unpackhi_twice(V1,V2);

			// V1 = V3 + V4 mod 2P
			// P2 = [2p 2p 2p 2p]
			V1 = add_mod<simd >(V3,V4,P2);

			// V2 = (V3+(2P-V4))alpha mod 2P
			// V7 =  (V3+(2P-V4)) mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V7 = reduce<simd >(V6, P2);

			// V4 = [D D H H L L P P ]
			V4 = simd::unpackhi_twice(V7,V7);

			V3 = mul_mod_half<simd,SimdComp>(V4,beta,P,betap);

			// V2=[* * D H * * L P]
			V2 = simd::template shuffle_twice<0xDD>(V3); // 0xDD = 221

			/* 3rd step */
			// At this time I have V1=[A B E F I J M N], V7=[C G * * K O * *], V2=[* * D H * * L P]
			// I need V3 = [A C E G I K M O], V4=[B D F H J L N P]
			V3 = simd::unpacklo_twice(V1,V7);
			V4 = simd::unpackhi_twice(V1,V2);

			// V1 = V3 + V4 mod 2P
			V1 = add_mod<simd >(V3,V4,P2);

			// V2 = V3 + (2P - V4) mod 2P
			V5 = simd::sub(V4,P2);
			V6 = simd::sub(V3,V5);
			V2 = reduce<simd >(V6, P2);

			// Result in    V1=[A C E G I K M O] V2=[B D F H J L N P]
			// Transform to V3=[A B C D I J K L],V4=[E F G H M N O P]
			V3 = simd::unpacklo_twice(V1,V2);
			V4 = simd::unpackhi_twice(V1,V2);

			// Transform to V1=[A B C D E F G H], V2=[I J K L M N O P]
			V1 = SimdComp::unpacklo128(V3,V4);
			V2 = SimdComp::unpackhi128(V3,V4);

			// Store
			simd::storeu(ABCDEFGH,V1);
			simd::storeu(IJKLMNOP,V2);


		}


	}; // FFT_butterflies<Field, 8>

}

#endif // __LINBOX_polynomial_fft_butterflies_H
