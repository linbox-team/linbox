/*
 * Copyright (C) 2019 Cyril Bouvier, Pascal Giorgi, Romain Lebreton
 *
 * Written by Cyril Bouvier <cyril.bouvier@lirmm.fr>
 *            Pascal Giorgi <pascal.giorgi@lirmm.fr>
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

#ifndef __LINBOX_fft_simd_H
#define __LINBOX_fft_simd_H

#include <iostream>
#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include "fflas-ffpack/fflas/fflas_simd.h"

namespace LinBox {

    template<class Field, class Simd,
             typename std::enable_if<Simd::template is_same_element<Field>::value>::type* = nullptr> 
    struct SimdFFT {
        using Element = typename Field::Element;
        using vect_t = typename Simd::vect_t;
        using Compute_t = typename Field::Compute_t;

        template <bool B, class T = void>
        using enable_if_t = typename std::enable_if<B, T>::type;

        template <class T = void>
        using enable_if_integral_t = enable_if_t<Field::is_elt_integral_v, T>;

        template <class T = void>
        using enable_if_floating_point_t = enable_if_t<Field::is_elt_floating_point_v, T>;

        template <class T = void>
        using enable_if_compute_same_size_t = enable_if_t<sizeof(Element) == sizeof(Compute_t), T>;

        template <class T = void>
        using enable_if_compute_twice_size_t = enable_if_t<sizeof(Element) != sizeof(Compute_t), T>;

        template <class T1, class T2, class T = void>
        using enable_if_same_t = enable_if_t<std::is_same<T1, T2>::value, T>;


        /**********************************************************************/
        /* reduce *************************************************************/
        /**********************************************************************/
        /* Reduce from [0..2p[ to [0..p[ */
        static inline vect_t
        reduce (const vect_t& a, const vect_t& p) {
            vect_t t = Simd::greater(p,a);
            return Simd::sub(a, Simd::vandnot(t,p));
        }

        /**********************************************************************/
        /* add_mod ************************************************************/
        /**********************************************************************/
        /* Add two operands in [0..p[ and return the result in [0..p[ */
        static inline vect_t
        add_mod (const vect_t& a, const vect_t& b, const vect_t& p) {
            vect_t c = Simd::add(a,b);
            return reduce(c, p);
        }

        /**********************************************************************/
        /* sub_mod ************************************************************/
        /**********************************************************************/
        /* Substract two operands in [0..p[ and return the result in [0..p[ */
        static inline vect_t
        sub_mod (const vect_t& a, const vect_t& b, const vect_t& p) {
            vect_t c = Simd::sub(p,b);
            c = Simd::add(a,c);
            return reduce(c, p);
        }

        /**********************************************************************/
        /* mul_mod ************************************************************/
        /**********************************************************************/
        /* input and output ranges depend on type of Element */

        /* mul mod for integral type */
        /* Element and Compute has same size, Element can store a full mul */
        template <class T=vect_t, enable_if_integral_t<T>* = nullptr,
                                  enable_if_compute_same_size_t<T>* = nullptr>
        static inline T
        mul_mod (const vect_t& a, const vect_t& b, const vect_t& p,
                 const vect_t& bp) {
#if 1
            vect_t q = Simd::mul (a, bp);
            q = Simd::template sra<4*sizeof (Element)> (q);
            vect_t c = Simd::mul (a, b);
            vect_t t = Simd::mul (q, p);
            return Simd::sub(c, t);
#else
            vect_t q = Simd::mulhi(a,bp);
            vect_t c = Simd::mullo(a,b);
            vect_t t = Simd::mullo(q,p);
            return Simd::sub(c,t);
#endif
        }

        /* Compute is twice as big as Compute */
        template <class T=vect_t, enable_if_integral_t<T>* = nullptr,
                                  enable_if_compute_twice_size_t<T>* = nullptr>
        static inline T
        mul_mod (const vect_t& a, const vect_t& b, const vect_t& p,
                 const vect_t& bp) {
            vect_t q = Simd::mulhi(a,bp);
            vect_t c = Simd::mullo(a,b);
            vect_t t = Simd::mullo(q,p);
            return Simd::sub(c,t);
        }

        /* mul mod for floating type */
        template <class T=vect_t, enable_if_floating_point_t<T>* = nullptr>
        static inline T
        mul_mod (const vect_t& x, const vect_t& y, const vect_t& p,
                 const vect_t& u) {
            // u = 1/p
            vect_t h = Simd::mul(x,y);
            vect_t l = Simd::fmsub(h,x,y); // Beware of the order!
            vect_t b = Simd::mul(h,u);
            vect_t c = Simd::floor(b);
            vect_t d = Simd::fnmadd(h,c,p); // Beware of the order!
            vect_t g = Simd::add(d,l);
            vect_t t = Simd::sub(g,p);
            g = Simd::blendv(t,g,t);
            t = Simd::add(g,p);
            return Simd::blendv(g,t,g);
        }

        /**********************************************************************/
        /* mul_mod_half *******************************************************/
        /**********************************************************************/
        /*
         * a = [a0, a0, a2, a2, ...]
         * b = [?, b0, ?, b2, ...]
         * bp = [b0p ? b2p ?, ... ], the Shoup mul_mod precomputation of b
         * Return [?, (a0*b0) mod p, ?, (a2*b2) mod p, ... ]
         */
        static inline vect_t
        mul_mod_half (const vect_t& a, const vect_t& b, const vect_t& p,
                      const vect_t& bp) {
            return mul_mod (a, b , p, bp);
#if 0
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

        /*
        //do not forget using Givaro::Modular;
        template<class T1, class T2, class Simd>
        static inline typename std::enable_if<std::is_floating_point<T1>::value, typename Simd::vect_t>::type
        SimdFFT<Modular<T1, T2>, Simd>::mul_mod (const typename Simd::vect_t& x, const typename Simd::vect_t& y, const typename Simd::vect_t& p,
                    const typename Simd::vect_t& u) {
                // u = 1/p
                // TODO If fixed argument y, we can save a mul
                typename Simd::vect_t xy = Simd::mul (x, y);
                typename Simd::vect_t xyu = Simd::mul (xy, u);
                typename Simd::vect_t q = Simd::floor (xyu);
                return Simd::fnmadd (xy, q, p);
        }*/

        /******************************************************************/
        /******************************************************************/
        /******************************************************************/
        /* unpacklohi:
            * Input (with n = Simd::vect_size):
            *  a = [ a0, a1, ..., an-1 ]
            *  b = [ b0, b1, ..., bn-1 ]
            * Ouput (with l = n/2 = Simd::vect_size/2):
            *  r1 = [ a0, b0, a1, b1, ..., al-1, bl-1 ]
            *  r2 = [ al, bl, al+1, bl+1, ..., an-1, bn-1 ]
            */
        /* Simd128<float> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<float>>* = nullptr>
        static inline void unpacklohi (vect_t& r1, vect_t& r2, const vect_t a,
                                                        const vect_t b) {
            r1 = _mm_unpacklo_ps (a, b);
            r2 = _mm_unpackhi_ps (a, b);
        }
        /* Simd128<double> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<double>>* = nullptr>
        static inline void unpacklohi (vect_t& r1, vect_t& r2, const vect_t a,
                                                        const vect_t b) {
            r1 = _mm_unpacklo_pd (a, b);
            r2 = _mm_unpackhi_pd (a, b);
        }
        /* Simd128<uint16_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<uint16_t>>* = nullptr>
        static inline void unpacklohi (vect_t& r1, vect_t& r2, const vect_t a,
                                                        const vect_t b) {
            r1 = _mm_unpacklo_epi16 (a, b);
            r2 = _mm_unpackhi_epi16 (a, b);
        }
        /* Simd128<uint32_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<uint32_t>>* = nullptr>
        static inline void unpacklohi (vect_t& r1, vect_t& r2, const vect_t a,
                                                        const vect_t b) {
            r1 = _mm_unpacklo_epi32 (a, b);
            r2 = _mm_unpackhi_epi32 (a, b);
        }
        /* Simd128<uint64_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<uint64_t>>* = nullptr>
        static inline void unpacklohi (vect_t& r1, vect_t& r2, const vect_t a,
                                                        const vect_t b) {
            r1 = _mm_unpacklo_epi64 (a, b);
            r2 = _mm_unpackhi_epi64 (a, b);
        }
        /* Simd256<float> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<float>>* = nullptr>
        static inline void
        unpacklohi (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
/* _mm256_permute4x64_pd needs AVX2 but Simd256<float> only needs AVX */
#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
            vect_t t1, t2;
            /* 0xd8 = 3120 base_4 */
            t1 = _mm256_castpd_ps (_mm256_permute4x64_pd
                                                (_mm256_castps_pd (a), 0xd8));
            t2 = _mm256_castpd_ps (_mm256_permute4x64_pd
                                                (_mm256_castps_pd (b), 0xd8));
            r1 = _mm256_unpacklo_ps (t1, t2);
            r2 = _mm256_unpackhi_ps (t1, t2);
#else /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS not defined */
            vect_t t1, t2;
            t1 = _mm256_unpacklo_ps (a, b);
            t2 = _mm256_unpackhi_ps (a, b);
            r1 = _mm256_permute2f128_ps (t1, t2, 0x20);
            r2 = _mm256_permute2f128_ps (t1, t2, 0x31);
#endif /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS */
        }
        /* Simd256<double> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<double>>* = nullptr>
        static inline void
        unpacklohi (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
/* _mm256_permute4x64_pd needs AVX2 but Simd256<double> only needs AVX */
#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
            vect_t t1, t2;
            /* 0xd8 = 3120 base_4 */
            t1 = _mm256_permute4x64_pd (a, 0xd8);
            t2 = _mm256_permute4x64_pd (b, 0xd8);
            r1 = _mm256_unpacklo_pd (t1, t2);
            r2 = _mm256_unpackhi_pd (t1, t2);
#else /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS not defined */
            vect_t t1, t2;
            t1 = _mm256_unpacklo_pd (a, b);
            t2 = _mm256_unpackhi_pd (a, b);
            r1 = _mm256_permute2f128_pd (t1, t2, 0x20);
            r2 = _mm256_permute2f128_pd (t1, t2, 0x31);
#endif /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS */
        }
        /* Simd256<uint16_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<uint16_t>>* = nullptr>
        static inline void
        unpacklohi (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            vect_t t1, t2;
            /* 0xd8 = 3120 base_4 */
            t1 = _mm256_permute4x64_epi64 (a, 0xd8);
            t2 = _mm256_permute4x64_epi64 (b, 0xd8);
            r1 = _mm256_unpacklo_epi16 (t1, t2);
            r2 = _mm256_unpackhi_epi16 (t1, t2);
        }
        /* Simd256<uint32_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<uint32_t>>* = nullptr>
        static inline void
        unpacklohi (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            vect_t t1, t2;
            /* 0xd8 = 3120 base_4 */
            t1 = _mm256_permute4x64_epi64 (a, 0xd8);
            t2 = _mm256_permute4x64_epi64 (b, 0xd8);
            r1 = _mm256_unpacklo_epi32 (t1, t2);
            r2 = _mm256_unpackhi_epi32 (t1, t2);
        }
        /* Simd256<uint64_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<uint64_t>>* = nullptr>
        static inline void
        unpacklohi (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            vect_t t1, t2;
            /* 0xd8 = 3120 base_4 */
            t1 = _mm256_permute4x64_epi64 (a, 0xd8);
            t2 = _mm256_permute4x64_epi64 (b, 0xd8);
            r1 = _mm256_unpacklo_epi64 (t1, t2);
            r2 = _mm256_unpackhi_epi64 (t1, t2);
        }

        /******************************************************************/
        /******************************************************************/
        /******************************************************************/
        /* pack:
            * Input (with n = Simd::vect_size):
            *  a = [ a0, a1, ..., an-1 ]
            *  b = [ b0, b1, ..., bn-1 ]
            * Ouput (with l = n/2 = Simd::vect_size/2):
            *  r1 = [ a0, a2, ..., an-2, b0, b2, ..., bn-2 ]
            *  r2 = [ a1, a3, ..., an-1, b1, b3, ..., bn-1 ]
            */
        /* Simd128<float> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<float>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            /* 0xd8 = 3120 base_4 */
            __m128d t1 = _mm_castps_pd (_mm_permute_ps (a, 0xd8));
            __m128d t2 = _mm_castps_pd (_mm_permute_ps (b, 0xd8));
            r1 = _mm_castpd_ps (_mm_unpacklo_pd (t1, t2));
            r2 = _mm_castpd_ps (_mm_unpackhi_pd (t1, t2));
        }
        /* Simd128<double> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<double>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            r1 = _mm_unpacklo_pd (a, b);
            r2 = _mm_unpackhi_pd (a, b);
        }
        /* Simd128<uint16_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<uint16_t>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            vect_t t1, t2, idx;
            /* 0x =  base_16 */
            idx = _mm_set_epi8 (15,14,11,10,7,6,3,2,13,12,9,8,5,4,1,0);
            t1 = _mm_shuffle_epi8 (a, idx);
            t2 = _mm_shuffle_epi8 (b, idx);
            r1 = _mm_unpacklo_epi64 (t1, t2);
            r2 = _mm_unpackhi_epi64 (t1, t2);
        }
        /* Simd128<uint32_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<uint32_t>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            vect_t t1, t2;
            /* 0xd8 = 3120 base_4 */
            t1 = _mm_shuffle_epi32 (a, 0xd8);
            t2 = _mm_shuffle_epi32 (b, 0xd8);
            r1 = _mm_unpacklo_epi64 (t1, t2);
            r2 = _mm_unpackhi_epi64 (t1, t2);
        }
        /* Simd128<uint64_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd128<uint64_t>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            r1 = _mm_unpacklo_epi64 (a, b);
            r2 = _mm_unpackhi_epi64 (a, b);
        }
        /* Simd256<float> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<float>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
/* _mm256_permute4x64_pd needs AVX2 but Simd256<float> only needs AVX */
#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
            /* 0xd8 = 3120 base_4 */
            __m256d t1 = _mm256_castps_pd (_mm256_permute_ps (a, 0xd8));
            __m256d t2 = _mm256_castps_pd (_mm256_permute_ps (b, 0xd8));
            __m256d p1 = _mm256_unpacklo_pd (t1, t2);
            __m256d p2 = _mm256_unpackhi_pd (t1, t2);
            /* 0xd8 = 3120 base_4 */
            r1 = _mm256_castpd_ps (_mm256_permute4x64_pd (p1, 0xd8));
            r2 = _mm256_castpd_ps (_mm256_permute4x64_pd (p2, 0xd8));
#else /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS not defined */
            /* 0xd8 = 3120 base_4 */
            __m256d pa = _mm256_castps_pd (_mm256_permute_ps (a, 0xd8));
            __m256d pb = _mm256_castps_pd (_mm256_permute_ps (b, 0xd8));
            __m256d t1 = _mm256_permute2f128_pd (pa, pb, 0x20);
            __m256d t2 = _mm256_permute2f128_pd (pa, pb, 0x31);
            r1 = _mm256_castpd_ps (_mm256_unpacklo_pd (t1, t2));
            r2 = _mm256_castpd_ps (_mm256_unpackhi_pd (t1, t2));

#endif /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS */
        }
        /* Simd256<double> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<double>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
/* _mm256_permute4x64_pd needs AVX2 but Simd256<double> only needs AVX */
#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
            r1 = _mm256_unpacklo_pd (a, b);
            r2 = _mm256_unpackhi_pd (a, b);
            /* 0xd8 = 3120 base_4 */
            r1 = _mm256_permute4x64_pd (r1, 0xd8);
            r2 = _mm256_permute4x64_pd (r2, 0xd8);
#else /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS not defined */
            vect_t t1, t2;
            t1 = _mm256_permute2f128_pd (a, b, 0x20);
            t2 = _mm256_permute2f128_pd (a, b, 0x31);
            r1 = _mm256_unpacklo_pd (t1, t2);
            r2 = _mm256_unpackhi_pd (t1, t2);
#endif /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS */
        }
        /* Simd256<uint32_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<uint32_t>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            vect_t t1, t2;
            /* 0xd8 = 3120 base_4 */
            r1 = _mm256_shuffle_epi32 (a, 0xd8);
            r2 = _mm256_shuffle_epi32 (b, 0xd8);
            t1 = _mm256_unpacklo_epi64 (r1, r2);
            t2 = _mm256_unpackhi_epi64 (r1, r2);
            /* 0xd8 = 3120 base_4 */
            r1 = _mm256_permute4x64_epi64 (t1, 0xd8);
            r2 = _mm256_permute4x64_epi64 (t2, 0xd8);
        }
        /* Simd256<uint64_t> */
        template <typename S = Simd,
                    enable_if_same_t<S, Simd256<uint64_t>>* = nullptr>
        static inline void
        pack (vect_t& r1, vect_t& r2, const vect_t a, const vect_t b) {
            r1 = _mm256_unpacklo_epi64 (a, b);
            r2 = _mm256_unpackhi_epi64 (a, b);
            /* 0xd8 = 3120 base_4 */
            r1 = _mm256_permute4x64_epi64 (r1, 0xd8);
            r2 = _mm256_permute4x64_epi64 (r2, 0xd8);
        }
    };
}

#endif /* __LINBOX_fft_simd_H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
