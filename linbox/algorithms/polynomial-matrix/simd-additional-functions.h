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

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
	template <typename Field>
	struct SimdCompute_t<Simd128<typename Field::Element>, Field> {
		using Compute_t = Simd128<typename Field::Compute_t>;
	};
#endif

#if defined(__FFLASFFPACK_HAVE_AVX_INSTRUCTIONS)
	template <typename Field>
	struct SimdCompute_t<Simd256<typename Field::Element>, Field> {
		using Compute_t = Simd256<typename Field::Compute_t>;
	};
#endif


#define Simd_vect typename Simd::vect_t

	/*
	 * Generic memory operations
	*/
	template<class T, class Simd = Simd<T>>
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

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
	template<class T>
	struct MemoryOp<T, Simd128<T>> {
		using simd = Simd128<T>;
		using simd_vect = typename simd::vect_t;

		/**************/
		/* load/store */
		/**************/
		static INLINE simd_vect load (const T* const p) {return simd::load(p);}
		static INLINE void store(T *p, simd_vect v) {return simd::store(p, v);}

		/*********************/
		/* Specific shuffles */
		/*********************/
		static INLINE simd_vect shuffletwice8_DD (simd_vect& s1) {
			using simd128_16 = Simd128<uint16_t>;
			using simd128_64 = Simd128<uint64_t>;
			//			std::cout << "Test shuffletwice8_DD :\n"; FFLAS::print<simd128_16>(std::cout,s1);
			simd_vect s2 = simd128_64::sll(s1,16);
			//			std::cout << "\n"; FFLAS::print<simd128_16>(std::cout,s2);
			//			std::cout << "\n"; FFLAS::print<simd128_16>(std::cout,simd128_16::template blend<0x44>(s1,s2)); std::cout << "\n\n";
			return simd128_16::template blend<0x44>(s1,s2); // 0x44 = [0 1 0 0 0 1 0 0]_base2
		}

		/********************/
		/* unpacklo         */
		/********************/
		static INLINE simd_vect unpacklo2 (const simd_vect& a, const simd_vect& b) {return Simd128<uint64_t>::unpacklo(a,b); }
		static INLINE simd_vect unpacklo4 (const simd_vect& a, const simd_vect& b) {return Simd128<uint32_t>::unpacklo(a,b); }
		static INLINE simd_vect unpacklo8 (const simd_vect& a, const simd_vect& b) {return Simd128<uint16_t>::unpacklo(a,b); }

		/********************/
		/* unpackhi         */
		/********************/
		static INLINE simd_vect unpackhi2 (const simd_vect& a, const simd_vect& b) {return Simd128<uint64_t>::unpackhi(a,b); }
		static INLINE simd_vect unpackhi4 (const simd_vect& a, const simd_vect& b) {return Simd128<uint32_t>::unpackhi(a,b); }
		static INLINE simd_vect unpackhi8 (const simd_vect& a, const simd_vect& b) {return Simd128<uint16_t>::unpackhi(a,b); }

		/**************/
		/* unpacklohi */
		/**************/
		static INLINE void unpacklohi2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_64 = Simd128<uint64_t>;
			s1 = simd128_64::unpacklo(a, b);
			s2 = simd128_64::unpackhi(a, b);
		}

		static INLINE void unpacklohi4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_32 = Simd128<uint32_t>;
			s1 = simd128_32::unpacklo(a, b);
			s2 = simd128_32::unpackhi(a, b);
		}

		static INLINE void unpacklohi8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_16 = Simd128<uint16_t>;
			s1 = simd128_16::unpacklo(a, b);
			s2 = simd128_16::unpackhi(a, b);
		}

		/********************/
		/* unpacklo_twice   */
		/********************/
		static INLINE simd_vect unpacklo_twice2 (const simd_vect& a, const simd_vect& b) { return unpacklo2(a,b); }

		static INLINE simd_vect unpacklo_twice4 (const simd_vect& a, const simd_vect& b) {
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			return simd128_32::unpacklo(a1,b1);
		}

		static INLINE simd_vect unpacklo_twice8 (const simd_vect& a, const simd_vect& b) {
			using simd128_16 = Simd128<uint16_t>;
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			return simd128_16::unpacklo(a1,b1);
		}

		/********************/
		/* unpackhi_twice   */
		/********************/
		static INLINE simd_vect unpackhi_twice2 (const simd_vect& a, const simd_vect& b) { return unpackhi2(a,b); }

		static INLINE simd_vect unpackhi_twice4 (const simd_vect& a, const simd_vect& b) {
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			return simd128_32::unpackhi(a1,b1);
		}

		static INLINE simd_vect unpackhi_twice8 (const simd_vect& a, const simd_vect& b) {
			using simd128_16 = Simd128<uint16_t>;
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			return simd128_16::unpackhi(a1,b1);
		}

		/********************/
		/* unpacklohi_twice */
		/********************/
		static INLINE void unpacklohi_twice2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			unpacklohi2(s1, s2, a, b);
		}

		static INLINE void unpacklohi_twice4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			s1 = simd128_32::unpacklo(a1,b1);
			s2 = simd128_32::unpackhi(a1,b1);
		}

		static INLINE void unpacklohi_twice8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd128_16 = Simd128<uint16_t>;
			using simd128_32 = Simd128<uint32_t>;
			simd_vect a1 = simd128_32::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd128_32::template shuffle<0xD8>(b);
			s1 = simd128_16::unpacklo(a1,b1);
			s2 = simd128_16::unpackhi(a1,b1);
		}
	}; // MemoryOp<T, Simd128<T>>
#endif

#if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
	template<class T>
	struct MemoryOp<T, Simd256<T>> {
		using simd = Simd256<T>;
		using simd_vect = typename simd::vect_t;

		/**************/
		/* load/store */
		/**************/
		static INLINE simd_vect load (const T* const p) {return simd::loadu(p);}
		static INLINE void store(T *p, simd_vect v) {return simd::storeu(p, v);}

		/*********************/
		/* Specific shuffles */
		/*********************/
		static INLINE simd_vect shuffletwice8_DD (simd_vect& s1) {
			using simd256_32 = Simd256<uint32_t>;
			return simd256_32::template shuffle_twice<0xDD>(s1);
		}

		/********************/
		/* unpacklo         */
		/********************/
		static INLINE simd_vect unpacklo2 (const simd_vect& a, const simd_vect& b) {return simd::unpacklo128(a, b); }

		template<typename V = simd_vect, typename std::enable_if<std::is_same<V, __m256d>::value>::type* = nullptr>
		static INLINE simd_vect unpacklo4 (const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<double>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_64::unpacklo_twice(a1,b1);
		}

		template<typename V = simd_vect, typename std::enable_if<std::is_same<V, __m256i>::value>::type* = nullptr>
		static INLINE simd_vect unpacklo4 (const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_64::unpacklo_twice(a1,b1);
		}

		static INLINE simd_vect unpacklo8 (const simd_vect& a, const simd_vect& b) {
			using simd256_32 = Simd256<uint32_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_32::unpacklo_twice(a1, b1);
		}

		static INLINE simd_vect unpacklo16 (const simd_vect& a, const simd_vect& b) {
			using simd256_16 = Simd256<uint16_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_16::unpacklo_twice(a1, b1);
		}

		/********************/
		/* unpackhi         */
		/********************/
		static INLINE simd_vect unpackhi2 (const simd_vect& a, const simd_vect& b) {return simd::unpackhi128(a, b); }

		template<typename V = simd_vect, typename std::enable_if<std::is_same<V, __m256d>::value>::type* = nullptr>
		static INLINE simd_vect unpackhi4 (const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<double>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_64::unpackhi_twice(a1,b1);
		}

		template<typename V = simd_vect, typename std::enable_if<std::is_same<V, __m256i>::value>::type* = nullptr>
		static INLINE simd_vect unpackhi4 (const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_64::unpackhi_twice(a1,b1);
		}

		static INLINE simd_vect unpackhi8 (const simd_vect& a, const simd_vect& b) {
			using simd256_32 = Simd256<uint32_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_32::unpackhi_twice(a1, b1);
		}

		static INLINE simd_vect unpackhi16 (const simd_vect& a, const simd_vect& b) {
			using simd256_16 = Simd256<uint16_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			return simd256_16::unpackhi_twice(a1, b1);
		}

		/**************/
		/* unpacklohi */
		/**************/
		static INLINE void unpacklohi2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			s1 = simd::unpacklo128(a, b);
			s2 = simd::unpackhi128(a, b);
		}

		template<typename V = simd_vect, typename std::enable_if<std::is_same<V, __m256i>::value>::type* = nullptr>
		static INLINE void unpacklohi4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			s1 = simd256_64::unpacklo_twice(a1, b1);
			s2 = simd256_64::unpackhi_twice(a1, b1);
		}

		template<typename V = simd_vect, typename std::enable_if<std::is_same<V, __m256d>::value>::type* = nullptr>
		static INLINE void unpacklohi4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<double>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			s1 = simd256_64::unpacklo_twice(a1, b1);
			s2 = simd256_64::unpackhi_twice(a1, b1);
		}

		static INLINE void unpacklohi8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_32 = Simd256<uint32_t>;
			using simd256_64 = Simd256<uint64_t>;
			simd_vect a1 = simd256_64::template shuffle<0xD8>(a); // 0xD8 = 3120 base_4
			simd_vect b1 = simd256_64::template shuffle<0xD8>(b);
			s1 = simd256_32::unpacklo_twice(a1, b1);
			s2 = simd256_32::unpackhi_twice(a1, b1);
		}

		static INLINE void unpacklohi16 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
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
		static INLINE simd_vect unpacklo_twice2 (const simd_vect& a, const simd_vect& b) { return unpacklo2(a,b); }

		static INLINE simd_vect unpacklo_twice4 (const simd_vect& a, const simd_vect& b) { return Simd256<uint64_t>::unpacklo_twice(a, b); }

		static INLINE simd_vect unpacklo_twice8 (const simd_vect& a, const simd_vect& b) { return Simd256<uint32_t>::unpacklo_twice(a, b); }

		static INLINE simd_vect unpacklo_twice16 (const simd_vect& a, const simd_vect& b) { return Simd256<uint16_t>::unpacklo_twice(a, b); }

		/********************/
		/* unpackhi_twice   */
		/********************/
		static INLINE simd_vect unpackhi_twice2 (const simd_vect& a, const simd_vect& b) { return unpackhi2(a,b); }

		static INLINE simd_vect unpackhi_twice4 (const simd_vect& a, const simd_vect& b) { return Simd256<uint64_t>::unpackhi_twice(a, b); }

		static INLINE simd_vect unpackhi_twice8 (const simd_vect& a, const simd_vect& b) { return Simd256<uint32_t>::unpackhi_twice(a, b); }

		static INLINE simd_vect unpackhi_twice16 (const simd_vect& a, const simd_vect& b) { return Simd256<uint16_t>::unpackhi_twice(a, b); }

		/********************/
		/* unpacklohi_twice */
		/********************/
		static INLINE void unpacklohi_twice2 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			unpacklohi2(s1, s2, a, b);
		}

		template<typename V = simd_vect, typename std::enable_if<std::is_same<V, __m256d>::value>::type* = nullptr>
		static INLINE void unpacklohi_twice4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<double>;
			s1 = simd256_64::unpacklo_twice(a, b);
			s2 = simd256_64::unpackhi_twice(a, b);
		}

		template<typename V = simd_vect, typename std::enable_if<std::is_same<V, __m256i>::value>::type* = nullptr>
		static INLINE void unpacklohi_twice4 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_64 = Simd256<uint64_t>;
			s1 = simd256_64::unpacklo_twice(a, b);
			s2 = simd256_64::unpackhi_twice(a, b);
		}

		static INLINE void unpacklohi_twice8 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_32 = Simd256<uint32_t>;
			s1 = simd256_32::unpacklo_twice(a, b);
			s2 = simd256_32::unpackhi_twice(a, b);
		}

		static INLINE void unpacklohi_twice16 (simd_vect& s1, simd_vect& s2, const simd_vect& a, const simd_vect& b) {
			using simd256_16 = Simd256<uint16_t>;
			s1 = simd256_16::unpacklo_twice(a, b);
			s2 = simd256_16::unpackhi_twice(a, b);
		}

	};// MemoryOp<T, Simd256<T>>
#endif


#define IS_INTEGRAL \
    typename std::enable_if<std::is_integral<typename Simd::scalar_t>::value>::type* = nullptr
#define IS_FLOATING \
    typename std::enable_if<std::is_floating_point<typename Simd::scalar_t>::value>::type* = nullptr


#define Simd_vect typename Simd::vect_t

	/*
	 * Generic arithmetic operation
	 */
	template <class Simd>
	INLINE Simd_vect reduce (const Simd_vect& a, const Simd_vect& p) {
		Simd_vect t = Simd::greater(p,a);
		return Simd::sub(a, Simd::vandnot(p,t));
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

	template <class Simd>
	INLINE Simd_vect sub_mod (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p) {
		Simd_vect c = Simd::sub(p,b);
		c = Simd::add(a,c);
		return reduce<Simd>(c, p);
	}

	template <class Simd, IS_INTEGRAL>
	INLINE Simd_vect mul_mod (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p, const Simd_vect& bp) {
		//		std::cout << "Inputs of mul_mod : a, b, p, bp, q, c, t, c - t\n";
		Simd_vect q = Simd::mulhi(a,bp);
		Simd_vect c = Simd::mullo(a,b);
		Simd_vect t = Simd::mullo(q,p);
		//		FFLAS::print<Simd>(std::cout, a); std::cout << "\n";
		//		FFLAS::print<Simd>(std::cout, b); std::cout << "\n";
		//		FFLAS::print<Simd>(std::cout, p); std::cout << "\n";
		//		FFLAS::print<Simd>(std::cout, bp); std::cout << "\n";
		//		FFLAS::print<Simd>(std::cout, q); std::cout << "\n";
		//		FFLAS::print<Simd>(std::cout, c); std::cout << "\n";
		//		FFLAS::print<Simd>(std::cout, t); std::cout << "\n";
		//		FFLAS::print<Simd>(std::cout, Simd::sub(c,t)); std::cout << "\n\n";
		return Simd::sub(c,t);
	}

	template <class Simd, IS_FLOATING>
	INLINE Simd_vect mul_mod (const Simd_vect& x, const Simd_vect& y, const Simd_vect& p, const Simd_vect& u) {
		// u = 1/p
		// std::cout << "Inputs of mul_mod : a, b, p, q, c, t, c - t\n";
		Simd_vect h = Simd::mul(x,y);
		Simd_vect l = Simd::fmsub(x,y,h);
		Simd_vect b = Simd::mul(h,u);
		Simd_vect c = Simd::floor(b);
		Simd_vect d = Simd::fnmadd(c,p,h);
		Simd_vect g = Simd::add(d,l);
		Simd_vect t = Simd::sub(g,p);
		g = Simd::blendv(t,g,t);
		t = Simd::add(g,p);
		return Simd::blendv(g,t,g);
	}

	/*
	 * a = [a0, a0, a2, a2, ...]
	* b = [?, b0, ?, b2, ...] with bp its shoup mul_mod precomputation [b0p ? b2p ?, ... ]
	* Return [?, (a0*b0) mod p, ?, (a2*b2) mod p, ... ]
	*/
	template <class Simd, class SimdCompute_t>
	INLINE Simd_vect mul_mod_half (const Simd_vect& a, const Simd_vect& b, const Simd_vect& p, const Simd_vect& bp) {
#if 1
		return mul_mod<Simd>(a, b , p, bp);
#else
		// TODO : DO SOMETHING IF Modular<uint64, uint128> and no mulx exits

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

#undef Simd_vect


#undef IS_INTEGRAL
#undef IS_FLOATING

}

#endif // __LINBOX_simd_additional_functions_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
