/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/*
 * Copyright (C) 2016 Romain Lebreton
 *
 * Written by Romain Lebreton <romain.lebreton@lirmm.fr>
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

#ifndef __LINBOX_simd_additional_functions_H
#define __LINBOX_simd_additional_functions_H

#include <iostream>
#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include "fflas-ffpack/fflas/fflas_simd.h"

namespace LinBox {


#define Simd_vect typename Simd::vect_t

	/*
	 * Generic arithmetic operation
	 */
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
#if 0
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


	/*
	 * Generic memory operations
	*/
	template<class T, class Simd = Simd<T>>
	struct MemoryOp {

		// Call load /store  (16 bits alignement)        if Simd128
		static inline Simd_vect load (const T* const p);

		// Call loadu/storeu (no alignement requirement) if Simd256
		static inline void store(T *p, Simd_vect v);

		static inline Simd_vect shuffletwice8_DD (Simd_vect& s1);

		static inline Simd_vect unpacklo2 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo4 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo8 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo16 (const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpackhi2 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi4 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi8 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi16 (const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpacklo_twice2 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo_twice4 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo_twice8 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklo_twice16 (const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpackhi_twice2 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi_twice4 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi_twice8 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpackhi_twice16 (const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpacklohi_twice2 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi_twice4 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi_twice8 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi_twice16 (const Simd_vect& a, const Simd_vect& b);

		static inline Simd_vect unpacklohi2 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi4 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi8 (const Simd_vect& a, const Simd_vect& b);
		static inline Simd_vect unpacklohi16 (const Simd_vect& a, const Simd_vect& b);

	}; // MemoryOp

#undef Simd_vect

#if defined(__FFLASFFPACK_USE_SIMD)
	template<class T>
	struct MemoryOp<T, Simd128<T>> {
		using simd = Simd128<T>;
		using simd_vect = typename simd::vect_t;

		/**************/
		/* load/store */
		/**************/
		static inline simd_vect load (const T* const p) {return simd::load(p);}
		static inline void store(T *p, simd_vect v) {return simd::store(p, v);}

		/*********************/
		/* Specific shuffles */
		/*********************/
		static inline simd_vect shuffletwice8_DD (simd_vect& s1) {
			using simd128_16 = Simd128<uint16_t>;
			using simd128_64 = Simd128<uint64_t>;
			simd_vect s2 = simd128_64::sll(s1,16);
			return simd128_16::template blend<0x44>(s1,s2); // 0x44 = [0 1 0 0 0 1 0 0]_base2
		}

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
#endif

#if defined(__FFLASFFPACK_USE_AVX2)
	template<class T>
	struct MemoryOp<T, Simd256<T>> {
		using simd = Simd256<T>;
		using simd_vect = typename simd::vect_t;

		/**************/
		/* load/store */
		/**************/
		static inline simd_vect load (const T* const p) {return simd::loadu(p);}
		static inline void store(T *p, simd_vect v) {return simd::storeu(p, v);}

		/*********************/
		/* Specific shuffles */
		/*********************/
		static inline simd_vect shuffletwice8_DD (simd_vect& s1) {
			using simd256_32 = Simd256<uint32_t>;
			return simd256_32::template shuffle_twice<0xDD>(s1);
		}

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
#endif

	template <typename simd, typename Field>
	struct SimdCompute_t {};

#if defined(__FFLASFFPACK_USE_SIMD)
	template <typename Field>
	struct SimdCompute_t<Simd128<typename Field::Element>, Field> {
		using Compute_t = Simd128<typename Field::Compute_t>;
	};
#endif

#if defined(__FFLASFFPACK_USE_AVX)
	template <typename Field>
	struct SimdCompute_t<Simd256<typename Field::Element>, Field> {
		using Compute_t = Simd256<typename Field::Compute_t>;
	};
#endif

}

#endif // __LINBOX_simd_additional_functions_H
