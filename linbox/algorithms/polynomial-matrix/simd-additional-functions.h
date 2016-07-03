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

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define INLINE __attribute__((always_inline)) inline
#else
#define INLINE inline
#endif

#if defined(__GNUC__) || defined(__clang__) || defined(__INTEL_COMPILER)
#define CONST __attribute__((const))
#else
#define CONST
#endif


namespace LinBox {


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


#define Simd_vect typename Simd::vect_t

	/*
	 * Generic memory operations
	*/
	template<class T, class Simd = Simd<T>, bool is_integral = std::is_integral<T>::value>
	struct MemoryOp {

		// Call load /store  (16 bits alignement)        if Simd128
		static INLINE Simd_vect load (const T* const p);

		// Call loadu/storeu (no alignement requirement) if Simd256
		static INLINE void store(T *p, Simd_vect v);

		static INLINE Simd_vect shuffletwice8_DD (Simd_vect& s1);

		static INLINE Simd_vect unpacklo2 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklo4 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklo8 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklo16 (const Simd_vect& a, const Simd_vect& b);

		static INLINE Simd_vect unpackhi2 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpackhi4 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpackhi8 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpackhi16 (const Simd_vect& a, const Simd_vect& b);

		static INLINE Simd_vect unpacklo_twice2 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklo_twice4 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklo_twice8 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklo_twice16 (const Simd_vect& a, const Simd_vect& b);

		static INLINE Simd_vect unpackhi_twice2 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpackhi_twice4 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpackhi_twice8 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpackhi_twice16 (const Simd_vect& a, const Simd_vect& b);

		static INLINE Simd_vect unpacklohi_twice2 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklohi_twice4 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklohi_twice8 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklohi_twice16 (const Simd_vect& a, const Simd_vect& b);

		static INLINE Simd_vect unpacklohi2 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklohi4 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklohi8 (const Simd_vect& a, const Simd_vect& b);
		static INLINE Simd_vect unpacklohi16 (const Simd_vect& a, const Simd_vect& b);

	}; // MemoryOp

#undef Simd_vect

#if defined(__FFLASFFPACK_USE_SIMD)
	template<class T>
	struct MemoryOp<T, Simd128<T>, true> {
		using simd = Simd128<T>;
		using simd_vect = typename simd::vect_t;
		using simd_2 = Simd128<uint64_t>;
		using simd_4 = Simd128<uint32_t>;
		using simd_8 = Simd128<uint16_t>;

		/**************/
		/* load/store */
		/**************/
		static INLINE simd_vect load (const T* const p) {
#ifndef NDEBUG
			assert(simd::valid(p));
#endif
			return simd::load(p);
		}
		static INLINE void store(T *p, simd_vect v) {
#ifndef NDEBUG
			assert(simd::valid(p));
#endif
			return simd::store(p, v);
		}

		/*********************/
		/* Specific shuffles */
		/*********************/
		static INLINE simd_vect shuffletwice8_DD (simd_vect& s1) {
			simd_vect s2 = simd_2::sll(s1,16);
			return simd_8::template blend<0x44>(s1,s2); // 0x44 = [0 1 0 0 0 1 0 0]_base2
		}

		/********************/
		/* unpacklo         */
		/********************/
		static INLINE simd_vect unpacklo2 (const simd_vect& a, const simd_vect& b) {return simd_2::unpacklo(a,b); }
		static INLINE simd_vect unpacklo4 (const simd_vect& a, const simd_vect& b) {return simd_4::unpacklo(a,b); }
		static INLINE simd_vect unpacklo8 (const simd_vect& a, const simd_vect& b) {return simd_8::unpacklo(a,b); }

		/********************/
		/* unpackhi         */
		/********************/
		static INLINE simd_vect unpackhi2 (const simd_vect& a, const simd_vect& b) {return simd_2::unpackhi(a,b); }
		static INLINE simd_vect unpackhi4 (const simd_vect& a, const simd_vect& b) {return simd_4::unpackhi(a,b); }
		static INLINE simd_vect unpackhi8 (const simd_vect& a, const simd_vect& b) {return simd_8::unpackhi(a,b); }

		/**************/
		/* unpacklohi */
		/**************/
		static INLINE void unpacklohi2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd_2::unpacklo(a, b);
			s2 = simd_2::unpackhi(a, b);
		}

		static INLINE void unpacklohi4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd_4::unpacklo(a, b);
			s2 = simd_4::unpackhi(a, b);
		}

		static INLINE void unpacklohi8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd_8::unpacklo(a, b);
			s2 = simd_8::unpackhi(a, b);
		}

		/********************/
		/* unpacklo_twice   */
		/********************/
		static INLINE simd_vect unpacklo_twice2 (const simd_vect& a, const simd_vect& b) { return unpacklo2(a,b); }

		static INLINE simd_vect unpacklo_twice4 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_4::unpacklo(a1,b1);
		}

		static INLINE simd_vect unpacklo_twice8 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_8::unpacklo(a1,b1);
		}

		/********************/
		/* unpackhi_twice   */
		/********************/
		static INLINE simd_vect unpackhi_twice2 (const simd_vect& a, const simd_vect& b) { return unpackhi2(a,b); }

		static INLINE simd_vect unpackhi_twice4 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_4::unpackhi(a1,b1);
		}

		static INLINE simd_vect unpackhi_twice8 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_8::unpackhi(a1,b1);
		}

		/********************/
		/* unpacklohi_twice */
		/********************/
		static INLINE void unpacklohi_twice2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			unpacklohi2(s1, s2, a, b);
		}

		static INLINE void unpacklohi_twice4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			s1 = simd_4::unpacklo(a1,b1);
			s2 = simd_4::unpackhi(a1,b1);
		}

		static INLINE void unpacklohi_twice8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			s1 = simd_8::unpacklo(a1,b1);
			s2 = simd_8::unpackhi(a1,b1);
		}
	}; // MemoryOp<T, Simd128<T>>

	template<class T>
	struct MemoryOp<T, Simd128<T>, false> {
		using simd = Simd128<T>;
		using simd_vect = typename simd::vect_t;
		using simd_2 = Simd128<double>;
		using simd_4 = Simd128<float>;

		/**************/
		/* load/store */
		/**************/
		static INLINE simd_vect load (const T* const p) {
#ifndef NDEBUG
			assert(simd::valid(p));
#endif
			return simd::load(p);
		}
		static INLINE void store(T *p, simd_vect v) {
#ifndef NDEBUG
			assert(simd::valid(p));
#endif
			return simd::store(p, v);
		}

		/********************/
		/* unpacklo         */
		/********************/
		static INLINE __m128d unpacklo2 (const __m128d& a, const __m128d& b) {return simd_2::unpacklo(a,b); }
		static INLINE __m128 unpacklo2 (const __m128& a, const __m128& b) {
			return _mm_castpd_ps(unpacklo2(_mm_castps_pd(a),_mm_castps_pd(b)));
		}
		static INLINE __m128 unpacklo4 (const __m128& a, const __m128& b) {return simd_4::unpacklo(a,b); }

		/********************/
		/* unpackhi         */
		/********************/
		static INLINE __m128d unpackhi2 (const __m128d& a, const __m128d& b) {return simd_2::unpackhi(a,b); }
		static INLINE __m128 unpackhi2 (const __m128& a, const __m128& b) {
			return _mm_castpd_ps(unpackhi2(_mm_castps_pd(a),_mm_castps_pd(b)));
		}
		static INLINE __m128 unpackhi4 (const __m128& a, const __m128& b) {return simd_4::unpackhi(a,b); }

		/**************/
		/* unpacklohi */
		/**************/
		static INLINE void unpacklohi2 (__m128d& s1, __m128d& s2, const __m128d& a, const __m128d& b) {
			s1 = simd_2::unpacklo(a, b);
			s2 = simd_2::unpackhi(a, b);
		}
		static INLINE void unpacklohi2 (__m128& s1, __m128& s2, const __m128& a, const __m128& b) {
			s1 = _mm_castpd_ps(simd_2::unpacklo(_mm_castps_pd(a), _mm_castps_pd(b)));
			s2 = _mm_castpd_ps(simd_2::unpackhi(_mm_castps_pd(a), _mm_castps_pd(b)));
		}
		static INLINE void unpacklohi4 (__m128& s1, __m128& s2, const __m128& a, const __m128& b) {
			s1 = simd_4::unpacklo(a, b);
			s2 = simd_4::unpackhi(a, b);
		}

		/********************/
		/* unpacklo_twice   */
		/********************/
		static INLINE simd_vect unpacklo_twice2 (const simd_vect& a, const simd_vect& b) { return unpacklo2(a,b); }
		static INLINE __m128 unpacklo_twice4 (const __m128& a, const __m128& b) {
			__m128 a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			__m128 b1 = simd_4::template shuffle<0xD8>(b);
			return simd_4::unpacklo(a1,b1);
		}

		/********************/
		/* unpackhi_twice   */
		/********************/
		static INLINE simd_vect unpackhi_twice2 (const simd_vect& a, const simd_vect& b) { return unpackhi2(a,b); }

		static INLINE __m128 unpackhi_twice4 (const __m128& a, const __m128& b) {
			__m128 a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			__m128 b1 = simd_4::template shuffle<0xD8>(b);
			return simd_4::unpackhi(a1,b1);
		}

		/********************/
		/* unpacklohi_twice */
		/********************/
		static INLINE void unpacklohi_twice2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			unpacklohi2(s1, s2, a, b);
		}

		static INLINE void unpacklohi_twice4 (__m128& s1, __m128& s2, const __m128& a, const __m128& b) {
			__m128 a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			__m128 b1 = simd_4::template shuffle<0xD8>(b);
			s1 = simd_4::unpacklo(a1,b1);
			s2 = simd_4::unpackhi(a1,b1);
		}

	}; // MemoryOp<T, Simd128<T>>

#endif

#if defined(__FFLASFFPACK_USE_AVX2)
	template<class T>
	struct MemoryOp<T, Simd256<T>, true> {
		using simd = Simd256<T>;
		using simd_vect = typename simd::vect_t;
		using simd_4 = Simd256<uint64_t>;
		using simd_8 = Simd256<uint32_t>;
		using simd_16 = Simd256<uint16_t>;

		/**************/
		/* load/store */
		/**************/
		static INLINE simd_vect load (const T* const p) {
#ifndef NDEBUG
			assert(simd::valid(p));
#endif
			return simd::load(p);
			// {return simd::loadu(p);}
		}
		static INLINE void store(T *p, simd_vect v) {
#ifndef NDEBUG
			assert(simd::valid(p));
#endif
			return simd::store(p, v);
			// {return simd::storeu(p, v);}
		}

		/*********************/
		/* Specific shuffles */
		/*********************/
		static INLINE simd_vect shuffletwice8_DD (simd_vect& s1) {
			return simd_8::template shuffle_twice<0xDD>(s1);
		}

		/********************/
		/* unpacklo         */
		/********************/
		static INLINE simd_vect unpacklo2 (const simd_vect& a, const simd_vect& b) {return simd::unpacklo128(a, b); }

		static INLINE simd_vect unpacklo4 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_4::unpacklo_twice(a1,b1);
		}

		static INLINE simd_vect unpacklo8 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_8::unpacklo_twice(a1, b1);
		}

		static INLINE simd_vect unpacklo16 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_16::unpacklo_twice(a1, b1);
		}

		/********************/
		/* unpackhi         */
		/********************/
		static INLINE simd_vect unpackhi2 (const simd_vect& a, const simd_vect& b) {return simd::unpackhi128(a, b); }

		static INLINE simd_vect unpackhi4 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_4::unpackhi_twice(a1,b1);
		}

		static INLINE simd_vect unpackhi8 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_8::unpackhi_twice(a1, b1);
		}

		static INLINE simd_vect unpackhi16 (const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			return simd_16::unpackhi_twice(a1, b1);
		}

		/**************/
		/* unpacklohi */
		/**************/
		static INLINE void unpacklohi2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd::unpacklo128(a, b);
			s2 = simd::unpackhi128(a, b);
		}

		static INLINE void unpacklohi4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			s1 = simd_4::unpacklo_twice(a1, b1);
			s2 = simd_4::unpackhi_twice(a1, b1);
		}

		static INLINE void unpacklohi8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			s1 = simd_8::unpacklo_twice(a1, b1);
			s2 = simd_8::unpackhi_twice(a1, b1);
		}

		static INLINE void unpacklohi16 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			simd_vect a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd_4::template shuffle<0xD8>(b);
			s1 = simd_16::unpacklo_twice(a1, b1);
			s2 = simd_16::unpackhi_twice(a1, b1);
		}

		/********************/
		/* unpacklo_twice   */
		/********************/
		static INLINE simd_vect unpacklo_twice2 (const simd_vect& a, const simd_vect& b) { return unpacklo2(a,b); }

		static INLINE simd_vect unpacklo_twice4 (const simd_vect& a, const simd_vect& b) { return simd_4::unpacklo_twice(a, b); }

		static INLINE simd_vect unpacklo_twice8 (const simd_vect& a, const simd_vect& b) { return simd_8::unpacklo_twice(a, b); }

		static INLINE simd_vect unpacklo_twice16 (const simd_vect& a, const simd_vect& b) { return simd_16::unpacklo_twice(a, b); }

		/********************/
		/* unpackhi_twice   */
		/********************/
		static INLINE simd_vect unpackhi_twice2 (const simd_vect& a, const simd_vect& b) { return unpackhi2(a,b); }

		static INLINE simd_vect unpackhi_twice4 (const simd_vect& a, const simd_vect& b) { return simd_4::unpackhi_twice(a, b); }

		static INLINE simd_vect unpackhi_twice8 (const simd_vect& a, const simd_vect& b) { return simd_8::unpackhi_twice(a, b); }

		static INLINE simd_vect unpackhi_twice16 (const simd_vect& a, const simd_vect& b) { return simd_16::unpackhi_twice(a, b); }

		/********************/
		/* unpacklohi_twice */
		/********************/
		static INLINE void unpacklohi_twice2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			unpacklohi2(s1, s2, a, b);
		}

		static INLINE void unpacklohi_twice4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd_4::unpacklo_twice(a, b);
			s2 = simd_4::unpackhi_twice(a, b);
		}

		static INLINE void unpacklohi_twice8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd_8::unpacklo_twice(a, b);
			s2 = simd_8::unpackhi_twice(a, b);
		}

		static INLINE void unpacklohi_twice16 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd_16::unpacklo_twice(a, b);
			s2 = simd_16::unpackhi_twice(a, b);
		}

	};// MemoryOp<T, Simd256<T>>

	template<class T>
	struct MemoryOp<T, Simd256<T>, false> {
		using simd = Simd256<T>;
		using simd_vect = typename simd::vect_t;
		using simd_4 = Simd256<double>;
		using simd_8 = Simd256<float>;

		/**************/
		/* load/store */
		/**************/
		static INLINE simd_vect load (const T* const p) {
#ifndef NDEBUG
			assert(simd::valid(p));
#endif
			return simd::load(p);
			// {return simd::loadu(p);}
		}
		static INLINE void store(T *p, simd_vect v) {
#ifndef NDEBUG
			assert(simd::valid(p));
#endif
			return simd::store(p, v);
			// {return simd::storeu(p, v);}
		}

		/*********************/
		/* Specific shuffles */
		/*********************/
		static INLINE __m256 shuffletwice8_DD (__m256& s1) {
			return simd_8::template shuffle_twice<0xDD>(s1);
		}

		/********************/
		/* unpacklo         */
		/********************/
		static INLINE simd_vect unpacklo2 (const simd_vect& a, const simd_vect& b) {return simd::unpacklo128(a, b); }

		static INLINE __m256d unpacklo4 (const __m256d& a, const __m256d& b) {
			__m256d a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			__m256d b1 = simd_4::template shuffle<0xD8>(b);
			return simd_4::unpacklo_twice(a1,b1);
		}

		static INLINE __m256 unpacklo4 (const __m256& a, const __m256& b) {
			__m256d a1 = simd_4::template shuffle<0xD8>(_mm256_castps_pd(a)); // 0xD8 = 3120 base_4
			__m256d b1 = simd_4::template shuffle<0xD8>(_mm256_castps_pd(b));
			return _mm256_castpd_ps(simd_4::unpacklo_twice(a1,b1));
		}

		static INLINE __m256 unpacklo8 (const __m256& a, const __m256& b) {
			__m256 a1 = _mm256_castpd_ps(simd_4::template shuffle<0xD8>(_mm256_castps_pd(a))); // 0xD8 = 3120 base_4
			__m256 b1 = _mm256_castpd_ps(simd_4::template shuffle<0xD8>(_mm256_castps_pd(b)));
			return simd_8::unpacklo_twice(a1, b1);
		}

		/********************/
		/* unpackhi         */
		/********************/
		static INLINE simd_vect unpackhi2 (const simd_vect& a, const simd_vect& b) {return simd::unpackhi128(a, b); }

		static INLINE __m256d unpackhi4 (const __m256d& a, const __m256d& b) {
			__m256d a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			__m256d b1 = simd_4::template shuffle<0xD8>(b);
			return simd_4::unpackhi_twice(a1,b1);
		}

		static INLINE __m256 unpackhi4 (const __m256& a, const __m256& b) {
			__m256d a1 = simd_4::template shuffle<0xD8>(_mm256_castps_pd(a)); // 0xD8 = 3120 base_4
			__m256d b1 = simd_4::template shuffle<0xD8>(_mm256_castps_pd(b));
			return _mm256_castpd_ps(simd_4::unpackhi_twice(a1,b1));
		}

		static INLINE __m256 unpackhi8 (const __m256& a, const __m256& b) {
			__m256 a1 = _mm256_castpd_ps(simd_4::template shuffle<0xD8>(_mm256_castps_pd(a))); // 0xD8 = 3120 base_4
			__m256 b1 = _mm256_castpd_ps(simd_4::template shuffle<0xD8>(_mm256_castps_pd(b)));
			return simd_8::unpackhi_twice(a1, b1);
		}

		/**************/
		/* unpacklohi */
		/**************/
		static INLINE void unpacklohi2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd::unpacklo128(a, b);
			s2 = simd::unpackhi128(a, b);
		}

		static INLINE void unpacklohi4 (__m256d& s1, __m256d& s2, const __m256d& a, const __m256d& b) {
			__m256d a1 = simd_4::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			__m256d b1 = simd_4::template shuffle<0xD8>(b);
			s1 = simd_4::unpacklo_twice(a1, b1);
			s2 = simd_4::unpackhi_twice(a1, b1);
		}

		static INLINE void unpacklohi4 (__m256& s1, __m256& s2, const __m256& a, const __m256& b) {
			__m256d a1 = simd_4::template shuffle<0xD8>(_mm256_castps_pd(a)); // 0xD8 = 3120 base_4
			__m256d b1 = simd_4::template shuffle<0xD8>(_mm256_castps_pd(b));
			s1 = _mm256_castpd_ps(simd_4::unpacklo_twice(a1, b1));
			s2 = _mm256_castpd_ps(simd_4::unpackhi_twice(a1, b1));
		}

		static INLINE void unpacklohi8 (__m256& s1, __m256& s2, const __m256& a, const __m256& b) {
			__m256 a1 = _mm256_castpd_ps(simd_4::template shuffle<0xD8>(_mm256_castps_pd(a))); // 0xD8 = 3120 base_4
			__m256 b1 = _mm256_castpd_ps(simd_4::template shuffle<0xD8>(_mm256_castps_pd(b)));
			s1 = simd_8::unpacklo_twice(a1, b1);
			s2 = simd_8::unpackhi_twice(a1, b1);
		}

		/********************/
		/* unpacklo_twice   */
		/********************/
		static INLINE simd_vect unpacklo_twice2 (const simd_vect& a, const simd_vect& b) { return unpacklo2(a,b); }

		static INLINE __m256d unpacklo_twice4 (const __m256d& a, const __m256d& b) {
			return simd_4::unpacklo_twice(a, b);
		}

		static INLINE __m256 unpacklo_twice4 (const __m256& a, const __m256& b) {
			return _mm256_castpd_ps(simd_4::unpacklo_twice(_mm256_castps_pd(a), _mm256_castps_pd(b)));
		}

		static INLINE __m256 unpacklo_twice8 (const __m256& a, const __m256& b) { return simd_8::unpacklo_twice(a, b); }

		/********************/
		/* unpackhi_twice   */
		/********************/
		static INLINE simd_vect unpackhi_twice2 (const simd_vect& a, const simd_vect& b) { return unpackhi2(a,b); }

		static INLINE __m256d unpackhi_twice4 (const __m256d& a, const __m256d& b) {
			return simd_4::unpackhi_twice(a, b);
		}

		static INLINE __m256 unpackhi_twice4 (const __m256& a, const __m256& b) {
			return _mm256_castpd_ps(simd_4::unpackhi_twice(_mm256_castps_pd(a), _mm256_castps_pd(b)));
		}

		static INLINE __m256 unpackhi_twice8 (const __m256& a, const __m256& b) { return simd_8::unpackhi_twice(a, b); }

		/********************/
		/* unpacklohi_twice */
		/********************/
		static INLINE void unpacklohi_twice2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			unpacklohi2(s1, s2, a, b);
		}

		static INLINE void unpacklohi_twice4 (__m256d& s1, __m256d& s2, const __m256d& a, const __m256d& b) {
			s1 = simd_4::unpacklo_twice(a, b);
			s2 = simd_4::unpackhi_twice(a, b);
		}

		static INLINE void unpacklohi_twice4 (__m256& s1, __m256& s2, const __m256& a, const __m256& b) {
			s1 = _mm256_castpd_ps(simd_4::unpacklo_twice(_mm256_castps_pd(a), _mm256_castps_pd(b)));
			s2 = _mm256_castpd_ps(simd_4::unpackhi_twice(_mm256_castps_pd(a), _mm256_castps_pd(b)));
		}

		static INLINE void unpacklohi_twice8 (__m256& s1, __m256& s2, const __m256& a, const __m256& b) {
			s1 = simd_8::unpacklo_twice(a, b);
			s2 = simd_8::unpackhi_twice(a, b);
		}

	};// MemoryOp<T, Simd256<T>>

#endif

#define Simd_vect typename Simd::vect_t

	/*
	 * Generic arithmetic operation
	 */

#define IS_INTEGRAL(Simd) std::enable_if<std::is_integral<typename Simd::scalar_t>::value>::type* = nullptr
#define IS_FLOATINGPOINT(Simd) std::enable_if<std::is_floating_point<typename Simd::scalar_t>::value>::type* = nullptr

	template <class Simd, typename IS_INTEGRAL(Simd)>
	INLINE Simd_vect reduce (const Simd_vect& a, const Simd_vect& p) {
		Simd_vect t = Simd::greater(p,a);
		return Simd::sub(a, Simd::vandnot(p,t));
	}

	template <class Simd, typename IS_FLOATINGPOINT(Simd)>
	INLINE Simd_vect reduce (const Simd_vect& a, const Simd_vect& p) {
		Simd_vect amp = Simd::sub(a,p);
		return Simd::blendv(amp, a, amp);
	}

	template <class Element, class Simd>
	INLINE void reduce (Element* a, const Simd_vect& p) {
		Simd_vect V1;
		V1 = MemoryOp<Element, Simd>::load(a);
		V1 = reduce<Simd>(V1, p);
		MemoryOp<Element, Simd>::store(a,V1);
	}

	template <class Simd>
	INLINE Simd_vect add_mod (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p) {
		Simd_vect c = Simd::add(a,b);
		return reduce<Simd>(c, p);
	}

	template <class Simd, typename IS_INTEGRAL(Simd)>
	INLINE Simd_vect mul_mod (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p,
							  const Simd_vect& bp, const Simd_vect& u = Simd::set1((typename Simd::scalar_t) 1)) {
		Simd_vect q = Simd::mulhi(a,bp);
		Simd_vect c = Simd::mullo(a,b);
		Simd_vect t = Simd::mullo(q,p);
		return Simd::sub(c,t);
	}

	template <class Simd, typename IS_FLOATINGPOINT(Simd)>
	INLINE Simd_vect mul_mod (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p,
							  const Simd_vect& bp, const Simd_vect& u = Simd::set1((typename Simd::scalar_t) 1)) {
		Simd_vect h = Simd::mul(a, b);
		Simd_vect l = Simd::fmsub(h, a, b);
		Simd_vect q = Simd::floor(Simd::mul(h, u));
		Simd_vect d = Simd::fnmadd(h, q, p);
		Simd_vect g = Simd::add(d, l);
		// if(g > p) g -= p;
		Simd_vect t = Simd::sub(g,p);
		g = Simd::blendv(t, g, t);
		// if(g < 0) g += p;
		t = Simd::add(g,p);
		return Simd::blendv(g, t, g);
	}

#define IS_INTEGRAL_AND_COMPUTET_INT128(Field)                              \
	std::enable_if<(std::is_integral<typename Field::Element>::value) &     \
	(std::is_same<typename Field::Compute_t,uint128_t>::value) >::type* = nullptr

#define IS_INTEGRAL_AND_COMPUTET_NOT_INT128(Field)                                  \
	std::enable_if<(std::is_integral<typename Field::Element>::value) &             \
	!(std::is_same<typename Field::Compute_t,uint128_t>::value) >::type* = nullptr

	/*
	 * a = [a0, a0, a2, a2, ...]
	* b = [?, b0, ?, b2, ...] with bp its shoup mul_mod precomputation [b0p ? b2p ?, ... ]
	* Return [?, (a0*b0) mod p, ?, (a2*b2) mod p, ... ]
	*/
	// Special template if COmpute_t == uint128_t since Simd128<uint128> and Simd256<uint128> do not exist
	template <typename Field, class Simd, typename IS_INTEGRAL_AND_COMPUTET_INT128(Field)>
	INLINE Simd_vect mul_mod_half (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p,
								   const Simd_vect& bp, const Simd_vect& u = Simd::set1((typename Simd::scalar_t) 1)) {
		return mul_mod<Simd>(a, b , p, bp);
	}

	template <typename Field, class Simd, typename IS_INTEGRAL_AND_COMPUTET_NOT_INT128(Field)>
	INLINE Simd_vect mul_mod_half (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p,
								   const Simd_vect& bp, const Simd_vect& u = Simd::set1((typename Simd::scalar_t) 1)) {
		using SimdComp = typename SimdCompute_t<Simd,Field>::Compute_t;
		// T2 = a * bp mod 2^64 (for Modular<Element = uint32, Compute_t = uint64>)
		// bp = [b0p ? b2p ?, ... ] is enough
		Simd_vect T2 = SimdComp::mulx(a,bp);
		Simd_vect T3 = Simd::mullo(T2,p);
		// At this point T3= [? quo(D)*p ? quo(H)*p] mod 2^32
		// T4 = [D D H H] * [?, b0, ?, b2] mod 2^32
		T2 = Simd::mullo(a,b);
		return Simd::sub(T2,T3);
	}

	template <typename Field, class Simd, typename IS_FLOATINGPOINT(Simd)>
	INLINE Simd_vect mul_mod_half (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p,
								   const Simd_vect& bp, const Simd_vect& u = Simd::set1((typename Simd::scalar_t) 1)) {
		return mul_mod<Simd>(a, b , p, bp, u);
	}

#undef IS_NOT_COMPUTET_INT128
#undef IS_COMPUTET_INT128
#undef IS_FLOATINGPOINT
#undef IS_INTEGRAL
#undef Simd_vect

}

#endif // __LINBOX_simd_additional_functions_H
