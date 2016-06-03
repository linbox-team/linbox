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

namespace LinBox {

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

	/* Memory operations
	*/
	template<class T, class Simd = Simd<T>>
	struct MemoryOp {

		// Call load /store  (16 bits alignement)        if Simd128
		static inline Simd_vect load (const T* const p);

		// Call loadu/storeu (no alignement requirement) if Simd256
		static inline void store(T *p, Simd_vect v);

		static inline Simd_vect unpacklo2 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo4 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo8 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo16 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpackhi2 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi4 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi8 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi16 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpacklo_twice2 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo_twice4 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo_twice8 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo_twice16 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpackhi_twice2 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi_twice4 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi_twice8 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi_twice16 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpacklohi_twice2 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi_twice4 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi_twice8 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi_twice16 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpacklohi2 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi4 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi8 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi16 (Simd_vect& s1, Simd_vect& s2, const Simd_vect& a, const Simd_vect& b);

	}; // MemoryOp

	template<class T>
	struct MemoryOp<T, Simd128<T>> {
		using simd = Simd128<T>;
		using simd_vect = typename simd::vect_t;

		/**************/
		/* load/store */
		/**************/
		static inline simd_vect load (const T* const p) {return simd::load(p);}
		static inline void store(T *p, simd_vect v) {return simd::store(p, v);}

		/********************/
		/* unpacklo         */
		/********************/
		static inline simd_vect unpacklo2 (const simd_vect& a, const simd_vect& b) {return Simd128<uint64_t>::unpacklo(a,b); }
		static inline simd_vect unpacklo4 (const simd_vect& a, const simd_vect& b) {return Simd128<uint32_t>::unpacklo(a,b); }
		static inline simd_vect unpacklo8 (const simd_vect& a, const simd_vect& b) {return Simd128<uint16_t>::unpacklo(a,b); }

		/********************/
		/* unpackhi         */
		/********************/
		static inline simd_vect unpackhi2 (const simd_vect& a, const simd_vect& b) {return Simd128<uint64_t>::unpackhi(a,b); }
		static inline simd_vect unpackhi4 (const simd_vect& a, const simd_vect& b) {return Simd128<uint32_t>::unpackhi(a,b); }
		static inline simd_vect unpackhi8 (const simd_vect& a, const simd_vect& b) {return Simd128<uint16_t>::unpackhi(a,b); }

		/**************/
		/* unpacklohi */
		/**************/
		static inline void unpacklohi2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_64 = Simd128<uint64_t>;
			s1 = simd128_64::unpacklo(a, b);
			s2 = simd128_64::unpackhi(a, b);
		}

		static inline void unpacklohi4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_32 = Simd128<uint32_t>;
			s1 = simd128_32::unpacklo(a, b);
			s2 = simd128_32::unpackhi(a, b);
		}

		static inline void unpacklohi8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_16 = Simd128<uint16_t>;
			s1 = simd128_16::unpacklo(a, b);
			s2 = simd128_16::unpackhi(a, b);
		}

		/********************/
		/* unpacklo_twice   */
		/********************/
		static inline simd_vect unpacklo_twice2 (const simd_vect& a, const simd_vect& b) { return unpacklo2(a,b); }

		static inline simd_vect unpacklo_twice4 (const simd_vect& a, const simd_vect& b) {
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			return simd128_32::unpacklo(a1,b1);
		}

		static inline simd_vect unpacklo_twice8 (const simd_vect& a, const simd_vect& b) {
			using simd128_16 = Simd128<uint16_t>;
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			return simd128_16::unpacklo(a1,b1);
		}

		/********************/
		/* unpackhi_twice   */
		/********************/
		static inline simd_vect unpackhi_twice2 (const simd_vect& a, const simd_vect& b) { return unpackhi2(a,b); }

		static inline simd_vect unpackhi_twice4 (const simd_vect& a, const simd_vect& b) {
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			return simd128_32::unpackhi(a1,b1);
		}

		static inline simd_vect unpackhi_twice8 (const simd_vect& a, const simd_vect& b) {
			using simd128_16 = Simd128<uint16_t>;
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			return simd128_16::unpackhi(a1,b1);
		}

		/********************/
		/* unpacklohi_twice */
		/********************/
		static inline void unpacklohi_twice2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			unpacklohi2(s1, s2, a, b);
		}

		static inline void unpacklohi_twice4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			s1 = simd128_32::unpacklo(a1,b1);
			s2 = simd128_32::unpackhi(a1,b1);
		}

		static inline void unpacklohi_twice8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_16 = Simd128<uint16_t>;
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			s1 = simd128_16::unpacklo(a1,b1);
			s2 = simd128_16::unpackhi(a1,b1);
		}
	}; // MemoryOp<T, Simd128<T>>


	template<class T>
	struct MemoryOp<T, Simd256<T>> {
		using simd = Simd256<T>;
		using simd_vect = typename simd::vect_t;

		/**************/
		/* load/store */
		/**************/
		static inline simd_vect load (const T* const p) {return simd::loadu(p);}
		static inline void store(T *p, simd_vect v) {return simd::storeu(p, v);}

		/********************/
		/* unpacklo         */
		/********************/
		static inline simd_vect unpacklo2 (const simd_vect& a, const simd_vect& b) {return simd::unpacklo128(a, b); }

		static inline simd_vect unpacklo4 (const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_64::unpacklo_twice(a1,b1);
		}

		static inline simd_vect unpacklo8 (const simd_vect& a, const simd_vect& b) {
			using simd256_32 = Simd256<uint32_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_32::unpacklo_twice(a1, b1);
		}

		static inline simd_vect unpacklo16 (const simd_vect& a, const simd_vect& b) {
			using simd256_16 = Simd256<uint16_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_16::unpacklo_twice(a1, b1);
		}

		/********************/
		/* unpackhi         */
		/********************/
		static inline simd_vect unpackhi2 (const simd_vect& a, const simd_vect& b) {return simd::unpackhi128(a, b); }

		static inline simd_vect unpackhi4 (const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_64::unpackhi_twice(a1,b1);
		}

		static inline simd_vect unpackhi8 (const simd_vect& a, const simd_vect& b) {
			using simd256_32 = Simd256<uint32_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_32::unpackhi_twice(a1, b1);
		}

		static inline simd_vect unpackhi16 (const simd_vect& a, const simd_vect& b) {
			using simd256_16 = Simd256<uint16_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_16::unpackhi_twice(a1, b1);
		}

		/**************/
		/* unpacklohi */
		/**************/
		static inline void unpacklohi2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd::unpacklo128(a, b);
			s2 = simd::unpackhi128(a, b);
		}

		static inline void unpacklohi4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			s1 = simd256_64::unpacklo_twice(a1, b1);
			s2 = simd256_64::unpackhi_twice(a1, b1);
		}

		static inline void unpacklohi8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_32 = Simd256<uint32_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			s1 = simd256_32::unpacklo_twice(a1, b1);
			s2 = simd256_32::unpackhi_twice(a1, b1);
		}

		static inline void unpacklohi16 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_16 = Simd256<uint16_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			s1 = simd256_16::unpacklo_twice(a1, b1);
			s2 = simd256_16::unpackhi_twice(a1, b1);
		}

		/********************/
		/* unpacklo_twice   */
		/********************/
		static inline simd_vect unpacklo_twice2 (const simd_vect& a, const simd_vect& b) { return unpacklo2(a,b); }

		static inline simd_vect unpacklo_twice4 (const simd_vect& a, const simd_vect& b) { return Simd256<uint64_t>::unpacklo_twice(a, b); }

		static inline simd_vect unpacklo_twice8 (const simd_vect& a, const simd_vect& b) { return Simd256<uint32_t>::unpacklo_twice(a, b); }

		static inline simd_vect unpacklo_twice16 (const simd_vect& a, const simd_vect& b) { return Simd256<uint16_t>::unpacklo_twice(a, b); }

		/********************/
		/* unpackhi_twice   */
		/********************/
		static inline simd_vect unpackhi_twice2 (const simd_vect& a, const simd_vect& b) { return unpackhi2(a,b); }

		static inline simd_vect unpackhi_twice4 (const simd_vect& a, const simd_vect& b) { return Simd256<uint64_t>::unpackhi_twice(a, b); }

		static inline simd_vect unpackhi_twice8 (const simd_vect& a, const simd_vect& b) { return Simd256<uint32_t>::unpackhi_twice(a, b); }

		static inline simd_vect unpackhi_twice16 (const simd_vect& a, const simd_vect& b) { return Simd256<uint16_t>::unpackhi_twice(a, b); }

		/********************/
		/* unpacklohi_twice */
		/********************/
		static inline void unpacklohi_twice2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			unpacklohi2(s1, s2, a, b);
		}

		static inline void unpacklohi_twice4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<uint64_t>;
			s1 = simd256_64::unpacklo_twice(a, b);
			s2 = simd256_64::unpackhi_twice(a, b);
		}

		static inline void unpacklohi_twice8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_32 = Simd256<uint32_t>;
			s1 = simd256_32::unpacklo_twice(a, b);
			s2 = simd256_32::unpackhi_twice(a, b);
		}

		static inline void unpacklohi_twice16 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_16 = Simd256<uint16_t>;
			s1 = simd256_16::unpacklo_twice(a, b);
			s2 = simd256_16::unpackhi_twice(a, b);
		}

	};// MemoryOp<T, Simd256<T>>

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
			vect_t V1,V2,V3,V4,T1,T2,T3,T4;
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

			MemoryOp<Element,simd>::unpacklohi4(V1,V2,V3,V4);

			// Second step
			// T1 = [D D H H]
			T1 = MemoryOp<Element,simd>::unpackhi4(V4,V4);

			T2 = mul_mod_half<simd, SimdComp>(T1, W, P, Wp);

			T2 = simd::template shuffle<0xDD>(T2);
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
			vect_t V1,V2,V3,V4,V5,V6,V7;
			// V1=[A B C D], V2=[E F G H]
			V1 = MemoryOp<Element,simd>::load(ABCD);
			V2 = MemoryOp<Element,simd>::load(EFGH);

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
			V3 = simd::template shuffle<0xDD>(V3); // 0xDD = [3 1 3 1]_base4

			// At this time I have V1=[A E B F], V2=[C G ? ?], V3=[? ? D H]
			// I need V3 = [A C E G], V4 = [B D F H]
			// This is not refactored in MemoryOp::... because of different arguments (V1,V3) and (V1,V2)
			V4 = MemoryOp<Element,simd>::unpackhi4(V1,V3);
			V3 = MemoryOp<Element,simd>::unpacklo4(V1,V2);

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
			V7 = simd::template shuffle_twice<0xDD>(V3); // 0xDD = 221 = [3 1 3 1]_base4

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

			// V2=[* * D H * * L P]
			V2 = simd::template shuffle_twice<0xDD>(V3); // 0xDD = 221 = [3 1 3 1]_base4

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
			MemoryOp<Element,simd>::unpacklohi_twice8(V3,V4,V1,V2);

			// Transform to V1=[A B C D E F G H], V2=[I J K L M N O P]
			MemoryOp<Element,simd>::unpacklohi2(V1,V2,V3,V4);

			// Store
			MemoryOp<Element,simd>::store(ABCDEFGH,V1);
			MemoryOp<Element,simd>::store(IJKLMNOP,V2);


		}


	}; // FFT_butterflies<Field, 8>

}

#endif // __LINBOX_polynomial_fft_butterflies_H
