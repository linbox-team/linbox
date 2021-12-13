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


#ifndef __LINBOX_fft_floating_INL
#define __LINBOX_fft_floating_INL

#include <type_traits>
#include <givaro/givpower.h>
#include "linbox/linbox-config.h"
#include "linbox/algorithms/polynomial-matrix/fft-simd.h"
#include "linbox/algorithms/polynomial-matrix/fft-utils.h"

//#define _FFT_NO_INLINED_BUTTERFLY 1

namespace LinBox {

    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
    template <typename Field, typename Simd>
    class FFT_base<Field, Simd,
                typename std::enable_if<Field::is_elt_floating_point_v
                    && Simd::template is_same_element<Field>::value>::type> {
        protected:
            /******************************************************************/
            /* Types **********************************************************/
            /******************************************************************/
            using Element = typename Field::Element;
            using simd_vect_t = typename Simd::vect_t;
            using SimdExtra = SimdFFT<Field, Simd>;

            /******************************************************************/
            /* Attributes *****************************************************/
            /******************************************************************/
            const Field *fld;
            size_t l2n; /* log2 of size */
            size_t n; /* 2^l2n */

            /* pow_w is the table of roots of unity. Its size is n-1.
             * If w = primitive n-th root, then the table is:
             *  1, w, w^2, ..., w^{n/2-1},        [ #elements = n/2 ]
             *  1, w^2, w^4, ..., w^{n/2-2},      [ #elements = n/4 ]
             *  1, w^4, w^8, ..., w^{n/2-4},      [ #elements = n/8 ]
             *  ...
             *  1, w^{n/8}, w^{n/4}, w^{3n/8},    [ #elements = 4 ]
             *  1, w^{n/4},                       [ #elements = 2 ]
             *  1.                                [ #elements = 1 ]
             *
             * pow_w_br is the same as pow_w with each subarray in bitreverse
             * order.
             */
            typename Simd::aligned_vector pow_w;
            typename Simd::aligned_vector pow_w_br;

            /******************************************************************/
            /* constructor ****************************************************/
            /******************************************************************/
            FFT_base (const Field& F, size_t k, Element w)
                                            : fld(&F), l2n(k), n(1UL << l2n),
                                              pow_w(n-1) , pow_w_br(n-1) {
                init_powers (w);
            }
        public:
            void
            DIF (Element *coeffs) const {
                /* w = n/2, f = 1 */
                DIF_core (coeffs, n >> 1, 1, pow_w.data(),
                                            FFTSimdHelper<Simd::vect_size>());
            }

            void
            DIT_reversed (Element *coeffs) const {
                /* w = n/2, f = 1 */
                DIT_reversed_core (coeffs, n >> 1, 1, pow_w_br.data()+ (n-2),
                                            FFTSimdHelper<Simd::vect_size>());
            }

            void
            DIT (Element *coeffs) const {
                /* w = 1, f = n / 2 */
                DIT_core (coeffs, 1, n >> 1, pow_w.data() + (n-2),
                                            FFTSimdHelper<Simd::vect_size>());
            }

            void
            DIF_reversed (Element *coeffs) const {
                /* w = 1, f = n / 2 */
                DIF_reversed_core (coeffs, 1, n >> 1, pow_w_br.data(),
                                            FFTSimdHelper<Simd::vect_size>());
            }

        protected:
            /******************************************************************/
            /* Core functions *************************************************/
            /******************************************************************/
            /* In the _core functions:
             *   n : length of the array 'coeffs' (always of power of 2)
             *   f : number of families of butterflies
             *   w : width of butterflies
             *   (outmost) loop invariant : 2*f*w == n
             */

            /* DIF ************************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIF_core (Element *coeffs, size_t w, size_t f,
                      const Element *pow, FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                for ( ; w >= Simd::vect_size; pow += w, f <<= 1, w >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w)
                        for (size_t j = 0; j < w; j += Simd::vect_size,
                                                    Aptr += Simd::vect_size,
                                                    Bptr += Simd::vect_size)
                            Butterfly_DIF (Aptr, Bptr, pow+j, P, U);
                }

                DIF_core_laststeps (coeffs, w, f, pow, P, U, h);
            }

            /* NoSimd */
            void
            DIF_core (Element *coeffs, size_t w, size_t f,
                      const Element *pow, FFTSimdHelper<1>) const {
                for ( ; w > 0; pow += w, f <<= 1, w >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w)
                        for (size_t j = 0; j < w; j++, Aptr++, Bptr++)
                            Butterfly_DIF (*Aptr, *Bptr, pow[j]);
                }
            }

            /* DIT reversed ***************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                            const Element *pow,
                                            FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                for ( ; w >= Simd::vect_size; f <<= 1, w >>= 1, pow -= f) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        for (size_t j = 0; j < w; j += Simd::vect_size,
                                                    Aptr += Simd::vect_size,
                                                    Bptr += Simd::vect_size)
                            Butterfly_DIT (Aptr, Bptr, alpha, P, U);
                    }
                }

                DIT_reversed_core_laststeps (coeffs, w, f, pow, P, U, h);
            }

            /* NoSimd */
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                                const Element *pow,
                                                FFTSimdHelper<1>) const {
                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w) {
                        Element alpha = pow[i];
                        for (size_t j = 0; j < w; j++, Aptr++, Bptr++)
                            Butterfly_DIT (*Aptr, *Bptr, alpha);
                    }
                }
            }

            /* DIT ************************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIT_core (Element *coeffs, size_t w, size_t f,
                      const Element *pow, FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                DIT_core_firststeps (coeffs, w, f, pow, P, U, h);

                for ( ; w < n; w <<= 1, f >>= 1, pow -= w) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w)
                        for (size_t j = 0; j < w; j += Simd::vect_size,
                                                    Aptr += Simd::vect_size,
                                                    Bptr += Simd::vect_size)
                            Butterfly_DIT (Aptr, Bptr, pow+j, P, U);
                }
            }

            /* NoSimd */
            void
            DIT_core (Element *coeffs, size_t w, size_t f,
                                        const Element *pow, FFTSimdHelper<1>,
                                        size_t bound = 0) const {
                bound = bound ? std::min (bound, n) : n;
                for ( ; w < bound; w <<= 1, f >>= 1, pow -= w) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w)
                        for (size_t j = 0; j < w; j++, Aptr++, Bptr++)
                            Butterfly_DIT (*Aptr, *Bptr, pow[j]);
                }
            }

            /* DIF reversed ***************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                            const Element *pow,
                                            FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                DIF_reversed_core_firststeps (coeffs, w, f, pow, P, U, h);

                for ( ; w < n; pow += f, w <<= 1, f >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        for (size_t j = 0; j < w; j += Simd::vect_size,
                                                    Aptr += Simd::vect_size,
                                                    Bptr += Simd::vect_size)
                            Butterfly_DIF (Aptr, Bptr, alpha, P, U);
                    }
                }
            }

            /* NoSimd */
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                        const Element *pow, FFTSimdHelper<1>,
                                        size_t bound = 0) const {
                bound = bound ? std::min (bound, n) : n;
                for ( ; w < bound; pow += f, w <<= 1, f >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w) {
                        Element alpha = pow[i];
                        for (size_t j = 0; j < w; j++, Aptr++, Bptr++)
                            Butterfly_DIF (*Aptr, *Bptr, alpha);
                    }
                }
            }

            /******************************************************************/
            /* Butterflies ****************************************************/
            /******************************************************************/
            /* Compute A[i]+B[i], (A[i]-B[i])*alpha[i] using Harvey's algorithm,
             * for 0 <= i < Simd::vect_size.
             * Input must satisfy:
             *  - 0 <= A[i],B[i],alpha[i] < p
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < p
             */
            void
            Butterfly_DIF (Element& A, Element& B, const Element& alpha) const {
                Element tmp;
                this->fld->assign (tmp, A);
                this->fld->addin (A, B);
                this->fld->sub (B, tmp, B);
                this->fld->mulin (B, alpha);
            }

            void
            Butterfly_DIF (simd_vect_t& A, simd_vect_t& B,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
                simd_vect_t T1, T2;

                /* A+B mod p */
                T1 = SimdExtra::add_mod (A, B, P);
                /* A-B mod p */
                T2 = SimdExtra::sub_mod (A, B, P);
                /* multiply A-B by alpha and store it in Bptr */
                B = SimdExtra::mul_mod (T2, alpha, P, U);

                A = T1;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                            const Element *alpha_ptr, const simd_vect_t& P,
                            const simd_vect_t& U) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3, alpha;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);
                alpha = Simd::load (alpha_ptr);

                /* A+B mod p */
                T1 = SimdExtra::add_mod (A, B, P);
                Simd::store (Aptr, T1);
                /* A-B mod p */
                T2 = SimdExtra::sub_mod (A, B, P);
                /* multiply A-B by alpha and store it in Bptr */
                T3 = SimdExtra::mul_mod (T2, alpha, P, U);
                Simd::store (Bptr, T3);
#else
                simd_vect_t alpha = Simd::load (alpha_ptr);
                Butterfly_DIF (Aptr, Bptr, alpha, P, U);
#endif
            }

            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* A+B mod p */
                T1 = SimdExtra::add_mod (A, B, P);
                Simd::store (Aptr, T1);
                /* A-B mod p */
                T2 = SimdExtra::sub_mod (A, B, P);
                /* multiply A-B by alpha and store it in Bptr */
                T3 = SimdExtra::mul_mod (T2, alpha, P, U);
                Simd::store (Bptr, T3);
#else
                simd_vect_t A = Simd::load (Aptr);
                simd_vect_t B = Simd::load (Bptr);
                Butterfly_DIF (A, B, alpha, P, U);
                Simd::store (Aptr, A);
                Simd::store (Bptr, B);
#endif
            }

            /* Compute A[i]+B[i]*alpha[i], A[i]-B[i]*alpha[i] using Harvey's
             * algorithm, for 0 <= i < simd::vect_size.
             * Input must satisfy:
             *  - 0 <= A[i],B[i],alpha[i] < p
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < p
             */
            void
            Butterfly_DIT (Element& A, Element& B, const Element& alpha) const {
                Element tmp;
                this->fld->mul (tmp, alpha, B);
                this->fld->sub (B, A, tmp);
                this->fld->addin (A, tmp);
            }

            void
            Butterfly_DIT (simd_vect_t& A, simd_vect_t& B,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
                simd_vect_t T1, T2;

                /* B*alpha mod P */
                T1 = SimdExtra::mul_mod (B, alpha, P, U);
                /* A+B*alpha */
                T2 = SimdExtra::add_mod (A, T1, P);
                /* A-B*alpha */
                B = SimdExtra::sub_mod (A, T1, P);

                A = T2;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                            const Element *alpha_ptr, const simd_vect_t& P,
                            const simd_vect_t& U) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3, alpha;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);
                alpha = Simd::load (alpha_ptr);

                /* B*alpha mod P */
                T1 = SimdExtra::mul_mod (B, alpha, P, U);
                /* A+B*alpha */
                T2 = SimdExtra::add_mod (A, T1, P);
                Simd::store (Aptr, T2);
                /* A-B*alpha */
                T3 = SimdExtra::sub_mod (A, T1, P);
                Simd::store (Bptr, T3);
#else
                simd_vect_t alpha = Simd::load (alpha_ptr);
                Butterfly_DIT (Aptr, Bptr, alpha, P, U);
#endif
            }

            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* B*alpha mod P */
                T1 = SimdExtra::mul_mod (B, alpha, P, U);
                /* A+B*alpha */
                T2 = SimdExtra::add_mod (A, T1, P);
                Simd::store (Aptr, T2);
                /* A-B*alpha */
                T3 = SimdExtra::sub_mod (A, T1, P);
                Simd::store (Bptr, T3);
#else
                simd_vect_t A = Simd::load (Aptr);
                simd_vect_t B = Simd::load (Bptr);
                Butterfly_DIT (A, B, alpha, P, U);
                Simd::store (Aptr, A);
                Simd::store (Bptr, B);
#endif
            }

            /******************************************************************/
            /* Laststeps for DIF and DIT reversed *****************************/
            /******************************************************************/

            /* Laststeps perform the last log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = Simd::vect_size/2, f = 4*n/Simd::vect_size, and the correct
             * pow pointer.
             */

            /* Fallback code, use no simd version of DIF_core */
            template<size_t VecSize>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const simd_vect_t& P,
                                const simd_vect_t& U,
                                FFTSimdHelper<VecSize>) const {
                /* P and U are unused with no simd fallback */
                DIF_core (coeffs, w, f, pow, FFTSimdHelper<1>());
            }

            /* For vect_size == 2 */
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const simd_vect_t& P,
                                const simd_vect_t& U,
                                FFTSimdHelper<2>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_core_laststeps (coeffs, w, f, pow, P, U,
                                                            FFTSimdHelper<1>());
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform V1 = [A B], V2 = [C D]
                         *      into V1 = [A C], V2 = [B D]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdExtra::add_mod (V1, V2, P);
                        V2 = SimdExtra::sub_mod (V1, V2, P);

                        /* Result in T = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        Simd::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            /* For vect_size == 4 */
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const simd_vect_t& P,
                                const simd_vect_t& U,
                                FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_core_laststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W, P, U);
                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdExtra::add_mod (V1, V2, P);
                        V2 = SimdExtra::sub_mod (V1, V2, P);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        Simd::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            /* For vect_size == 8 */
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const simd_vect_t& P,
                                const simd_vect_t& U,
                                FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_core_laststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1],
                                               pow[2], pow[2], pow[3], pow[3]);
                    simd_vect_t W2 = Simd::set (pow[4], pow[4], pow[4], pow[4],
                                                pow[5], pow[5], pow[5], pow[5]);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** step *********************************************/
                        Butterfly_DIF (V1, V2, W, P, U);
                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W2, P, U);
                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdExtra::add_mod (V1, V2, P);
                        V2 = SimdExtra::sub_mod (V1, V2, P);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        Simd::unpacklohi (V1, V2, T, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            /* Fallback code, use no simd version of DIT_reversed_core */
            template<size_t VecSize>
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<VecSize>) const {
                /* P and U are unused with no_simd fallback */
                DIT_reversed_core (coeffs, w, f, pow, FFTSimdHelper<1>());
            }

            /* For vect_size == 2 */
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const simd_vect_t& P,
                                const simd_vect_t& U,
                                FFTSimdHelper<2>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_reversed_core_laststeps (coeffs, w, f, pow, P, U,
                                                            FFTSimdHelper<1>());
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, W;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::load (pow+i);

                        /* transform V1 = [A B], V2 = [C D]
                         *      into V1 = [A C], V2 = [B D]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step ****************************************/
                        Butterfly_DIT (V1, V2, W, P, U);

                        /* Result in T = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            /* For vect_size == 4 */
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const simd_vect_t& P,
                                const simd_vect_t& U,
                                FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_reversed_core_laststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    const Element *pow0 = pow - (f << 1);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr,
                                                    pow0 += Simd::vect_size) {
                        simd_vect_t V1, V2, W;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::set (pow[i], pow[i+1], pow[i], pow[i+1]);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIT (V1, V2, W, P, U);
                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step ****************************************/
                        W = Simd::load (pow0);
                        Butterfly_DIT (V1, V2, W, P, U);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            /* For vect_size == 8 */
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const simd_vect_t& P,
                                const simd_vect_t& U,
                                FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_reversed_core_laststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    const Element *pow1 = pow - (f << 1);
                    const Element *pow0 = pow1 - (f << 2);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr,
                                                    pow0 += Simd::vect_size) {
                        simd_vect_t V1, V2, W;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::set (pow[i], pow[i+1], pow[i], pow[i+1],
                                            pow[i], pow[i+1], pow[i], pow[i+1]);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** step *********************************************/
                        Butterfly_DIT (V1, V2, W, P, U);
                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        W = Simd::set (pow1[2*i], pow1[2*i+1], pow1[2*i+2],
                                                                    pow1[2*i+3],
                                       pow1[2*i], pow1[2*i+1], pow1[2*i+2],
                                                                pow1[2*i+3]);
                        Butterfly_DIT (V1, V2, W, P, U);
                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        W = Simd::load (pow0);
                        Butterfly_DIT (V1, V2, W, P, U);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            /******************************************************************/
            /* Firststeps for DIT and DIF reversed ****************************/
            /******************************************************************/

            /* Firststeps perform the first log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = 1, f = n/2, and the correct pow pointer.
             */

            /* Fallback code, use no simd version of DIT_core */
            template<size_t VecSize>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<VecSize>) const {
                /* P and U are unused with no_simd fallback */
                DIT_core (coeffs, w, f, pow, FFTSimdHelper<1>(),
                                                            Simd::vect_size);
                for ( ; w < Simd::vect_size; w <<= 1, f >>= 1, pow -= w) ;
            }

            /* For vect_size == 2 */
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<2>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_core_firststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform V1 = [A B], V2 = [C D]
                         *      into V1 = [A C], V2 = [B D]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdExtra::add_mod (V1, V2, P);
                        V2 = SimdExtra::sub_mod (V1, V2, P);

                        /* Result in T = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        Simd::pack (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    w <<= 1;
                    f >>= 1;
                    pow -= w;
                }
            }

            /* For vect_size == 4 */
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_core_firststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    f >>= 2;
                    w <<= 2;
                    pow -= 2;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        Simd::pack (V1, V2, V1, V2);


                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdExtra::add_mod (V1, V2, P);
                        V2 = SimdExtra::sub_mod (V1, V2, P);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, P, U);

                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        Simd::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                }
            }

            /* For vect_size == 8 */
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_core_firststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    f >>= 3;
                    w <<= 3;
                    pow -= 2;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[0], pow[0],
                                               pow[1], pow[1], pow[1], pow[1]);
                    pow -= 4;
                    simd_vect_t W2 = Simd::set (pow[0], pow[0], pow[1], pow[1],
                                                pow[2], pow[2], pow[3], pow[3]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdExtra::add_mod (V1, V2, P);
                        V2 = SimdExtra::sub_mod (V1, V2, P);

                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        Simd::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, P, U);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** third step ***************************************/
                        Butterfly_DIT (V1, V2, W2, P, U);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                }
            }


            /* Fallback code, use no simd version of DIF_reversed_core */
            template<size_t VecSize>
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<VecSize>) const {
                /* P and U are unused with no_simd fallback */
                DIF_reversed_core (coeffs, w, f, pow, FFTSimdHelper<1>(),
                                                            Simd::vect_size);
                for ( ; w < Simd::vect_size; pow+=f, w <<= 1, f >>= 1) ;
            }

            /* For vect_size == 2 */
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<2>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_reversed_core_firststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, W;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::load (pow+i);

                        /* transform V1 = [A B], V2 = [C D]
                         *      into V1 = [A C], V2 = [B D]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** first step ***************************************/
                        Butterfly_DIF (V1, V2, W, P, U);

                        /* Result in V1 = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        Simd::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow += f;
                    w <<= 1;
                    f >>= 1;
                }
            }

            /* For vect_size == 4 */
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_reversed_core_firststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    const Element *pow0 = pow;
                    pow += f;
                    f >>= 1;
                    w <<= 2;
                    for (size_t i = 0; i < f; i += 2, coeffs += incr,
                                                    pow0 += Simd::vect_size) {
                        simd_vect_t V1, V2, W;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::load (pow0);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        Simd::pack (V1, V2, V1, V2);


                        /*** first step ***************************************/
                        Butterfly_DIF (V1, V2, W, P, U);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** second step **************************************/
                        W = Simd::set (pow[i], pow[i+1], pow[i], pow[i+1]);
                        Butterfly_DIF (V1, V2, W, P, U);

                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        Simd::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow += f;
                    f >>= 1;
                }
            }

            /* For vect_size == 8 */
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U,
                                    FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_reversed_core_firststeps (coeffs, w, f, pow, P, U,
                                                        FFTSimdHelper<1>());
                } else {
                    const Element *pow0 = pow;
                    pow += f;
                    f >>= 1;
                    const Element *pow1 = pow;
                    pow += f;
                    f >>= 1;
                    w <<= 3;
                    for (size_t i = 0; i < f; i += 2, coeffs += incr,
                                                    pow0 += Simd::vect_size) {
                        simd_vect_t V1, V2, W;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::load (pow0);

                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** first step ***************************************/
                        Butterfly_DIF (V1, V2, W, P, U);

                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** second step **************************************/
                        W = Simd::set (pow1[2*i], pow1[2*i+1], pow1[2*i+2],
                                                                    pow1[2*i+3],
                                       pow1[2*i], pow1[2*i+1], pow1[2*i+2],
                                                                pow1[2*i+3]);
                        Butterfly_DIF (V1, V2, W, P, U);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** third step ***************************************/
                        W = Simd::set (pow[i], pow[i+1], pow[i], pow[i+1],
                                            pow[i], pow[i+1], pow[i], pow[i+1]);
                        Butterfly_DIF (V1, V2, W, P, U);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow += f;
                    f >>= 1;
                }
            }

            /******************************************************************/
            /* Utils **********************************************************/
            /******************************************************************/
            void
            init_powers (const Element & w) {
                /* compute w^i and set first subarray */
                this->fld->assign (pow_w[0], this->fld->one);
                for (size_t i = 1; i < n/2; i++)
                    this->fld->mul (pow_w[i], pow_w[i-1], w);

                /* Other elements can be set from previously computed values */
                size_t idx = n/2; /* index for next value to be written in */
                for (size_t k=2; k <= n/2; k<<=1)
                    for(size_t i = 0; i < n/2; i+=k, idx++)
                        pow_w[idx]  = pow_w[i];

                /* init powers in bitreverse order */
                size_t l = l2n, len = n >> 1, base_idx = 0;
                do
                {
                    l--;
                    for (size_t i = 0; i < len; i++)
                    {
                        size_t i_br = FFT_utils::bitreverse (i, l);
                        pow_w_br[base_idx + i] = pow_w[base_idx + i_br];
                    }
                    base_idx += len;
                    len >>= 1;
                }
                while (l) ;
            }
    };

    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
    template <typename Field, typename Simd>
    class FFT_multi_base<Field, Simd,
                typename std::enable_if<Field::is_elt_floating_point_v
                    && Simd::template is_same_element<Field>::value>::type> {
        protected:
            /******************************************************************/
            /* Types **********************************************************/
            /******************************************************************/
            using Element = typename Field::Element;
            using Residu_t = typename Field::Residu_t;
            using simd_vect_t = typename Simd::vect_t;
            using SimdExtra = SimdFFT<Field, Simd>;

            /******************************************************************/
            /* Attributes *****************************************************/
            /******************************************************************/
            const Field *fld;
            size_t l2n; /* log2 of size */
            size_t n; /* 2^l2n */

            /* pow_w is the table of roots of unity. Its size is n-1.
             * If w = primitive n-th root, then the table is:
             *  1, w, w^2, ..., w^{n/2-1},        [ #elements = n/2 ]
             *  1, w^2, w^4, ..., w^{n/2-2},      [ #elements = n/4 ]
             *  1, w^4, w^8, ..., w^{n/2-4},      [ #elements = n/8 ]
             *  ...
             *  1, w^{n/8}, w^{n/4}, w^{3n/8},    [ #elements = 4 ]
             *  1, w^{n/4},                       [ #elements = 2 ]
             *  1.                                [ #elements = 1 ]
             *
             * pow_w_br is the same as pow_w with each subarray in bitreverse
             * order.
             *
             * pow_wp[i] := precomp_b (pow_w[i])
             * pow_wp_br[i] := precomp_b (pow_w_br[i])
             */
            typename Simd::aligned_vector pow_w;
            typename Simd::aligned_vector pow_w_br;

            /******************************************************************/
            /* constructor ****************************************************/
            /******************************************************************/
            FFT_multi_base (const Field& F, size_t k, Element w)
                                            : fld(&F), l2n(k), n(1UL << l2n),
                                              pow_w(n-1), pow_w_br(n-1) {
                init_powers (w);
            }

        public:
            void
            DIF (Element *coeffs, size_t stride) const {
                /* w = n/2, f = 1 */
                DIF_core (coeffs, n >> 1, 1, stride, pow_w.data(),
                                            FFTSimdHelper<Simd::vect_size>());
            }

            void
            DIT_reversed (Element *coeffs, size_t stride) const {
                /* w = n/2, f = 1 */
                DIT_reversed_core (coeffs, n >> 1, 1, stride,
                                            pow_w_br.data()+ (n-2),
                                            FFTSimdHelper<Simd::vect_size>());
            }

            void
            DIT (Element *coeffs, size_t stride) const {
                /* w = 1, f = n / 2 */
                DIT_core (coeffs, 1, n >> 1, stride, pow_w.data() + (n-2),
                                            FFTSimdHelper<Simd::vect_size>());
            }

            void
            DIF_reversed (Element *coeffs, size_t stride) const {
                /* w = 1, f = n / 2 */
                DIF_reversed_core (coeffs, 1, n >> 1, stride, pow_w_br.data(),
                                            FFTSimdHelper<Simd::vect_size>());
            }

       protected:
            /******************************************************************/
            /* Core functions *************************************************/
            /******************************************************************/
            /* In the _core functions:
             *   n : length of the array 'coeffs' (always of power of 2)
             *   f : number of families of butterflies
             *   w : width of butterflies
             *   (outmost) loop invariant : 2*f*w == n
             */

            /* DIF ************************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIF_core (Element *coeffs, size_t w, size_t f, size_t stride,
                                            const Element *pow,
                                            FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                for ( ; w > 0; pow += w, f <<= 1, w >>= 1) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws)
                        for (size_t j = 0; j < w; j += 1) {
                            simd_vect_t alpha = Simd::set1 (pow[j]);
                            size_t l = 0;
                            for ( ; l + Simd::vect_size <= stride;
                                                        l += Simd::vect_size,
                                                        Aptr += Simd::vect_size,
                                                        Bptr += Simd::vect_size)
                                Butterfly_DIF (Aptr, Bptr, alpha, P, U);
                            for ( ; l < stride; l++, Aptr++, Bptr++)
                                Butterfly_DIF (*Aptr, *Bptr, pow[j]);
                        }
                }
            }

            /* NoSimd */
            void
            DIF_core (Element *coeffs, size_t w, size_t f, size_t stride,
                                                    const Element *pow,
                                                    FFTSimdHelper<1>) const {
                for ( ; w > 0; pow += w, f <<= 1, w >>= 1) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws)
                        for (size_t j = 0; j < w; j += 1)
                            for (size_t l = 0; l < stride; l++, Aptr++, Bptr++)
                                Butterfly_DIF (*Aptr, *Bptr, pow[j]);
                }
            }

            /* DIT reversed ***************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                            size_t stride, const Element *pow,
                                            FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        size_t j = 0;
                        for ( ; j + Simd::vect_size <= ws; j += Simd::vect_size,
                                                        Aptr += Simd::vect_size,
                                                        Bptr += Simd::vect_size)
                            Butterfly_DIT (Aptr, Bptr, alpha, P, U);
                        for ( ; j < ws; j++, Aptr++, Bptr++)
                            Butterfly_DIT (*Aptr, *Bptr, pow[i]);
                    }
                }
            }

            /* NoSimd */
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                            size_t stride, const Element *pow,
                                            FFTSimdHelper<1>) const {
                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws) {
                        Element alpha = pow[i];
                        for (size_t j = 0; j < ws; j += 1, Aptr++, Bptr++)
                            Butterfly_DIT (*Aptr, *Bptr, alpha);
                    }
                }
            }

            /* DIT ************************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIT_core (Element *coeffs, size_t w, size_t f, size_t stride,
                        const Element *pow, FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                for ( ; w < n; w <<= 1, f >>= 1, pow -= w) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws)
                        for (size_t j = 0; j < w; j += 1) {
                            simd_vect_t alpha = Simd::set1 (pow[j]);
                            size_t l = 0;
                            for ( ; l + Simd::vect_size <= stride;
                                                        l += Simd::vect_size,
                                                        Aptr += Simd::vect_size,
                                                        Bptr += Simd::vect_size)
                                Butterfly_DIT (Aptr, Bptr, alpha, P, U);
                            for ( ; l < stride; l++, Aptr++, Bptr++)
                                Butterfly_DIT (*Aptr, *Bptr, pow[j]);
                        }
                }
            }

            /* NoSimd */
            void
            DIT_core (Element *coeffs, size_t w, size_t f, size_t stride,
                                const Element *pow, FFTSimdHelper<1>) const {
                for ( ; w < n; w <<= 1, f >>= 1, pow -= w) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws)
                        for (size_t j = 0; j < w; j += 1)
                            for (size_t l = 0; l < stride; l++, Aptr++, Bptr++)
                                Butterfly_DIT (*Aptr, *Bptr, pow[j]);
                }
            }

            /* DIF reversed ***************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                            size_t stride, const Element *pow,
                                            FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                for ( ; w < n; pow += f, w <<= 1, f >>= 1) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        size_t j = 0;
                        for ( ; j + Simd::vect_size <= ws; j += Simd::vect_size,
                                                        Aptr += Simd::vect_size,
                                                        Bptr += Simd::vect_size)
                            Butterfly_DIF (Aptr, Bptr, alpha, P, U);
                        for ( ; j < ws; j++, Aptr++, Bptr++)
                            Butterfly_DIF (*Aptr, *Bptr, pow[i]);
                    }
                }
            }

            /* NoSimd */
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                            size_t stride, const Element *pow,
                                            FFTSimdHelper<1>) const {
                for ( ; w < n; pow += f, w <<= 1, f >>= 1) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws) {
                        Element alpha = pow[i];
                        for (size_t j = 0; j < ws; j += 1, Aptr++, Bptr++)
                                Butterfly_DIF (*Aptr, *Bptr, alpha);
                    }
                }
            }

            /******************************************************************/
            /* Butterflies ****************************************************/
            /******************************************************************/
            /* Compute A[i]+B[i], (A[i]-B[i])*alpha[i] using Harvey's algorithm,
             * for 0 <= i < Simd::vect_size.
             * Input must satisfy:
             *  - 0 <= A[i],B[i],alpha[i] < p
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < p
             */
            void
            Butterfly_DIF (Element& A, Element& B, const Element& alpha) const {
                Element tmp;
                this->fld->assign (tmp, A);
                this->fld->addin (A, B);
                this->fld->sub (B, tmp, B);
                this->fld->mulin (B, alpha);
            }

            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::loadu (Aptr);
                B = Simd::loadu (Bptr);

                /* A+B mod p */
                T1 = SimdExtra::add_mod (A, B, P);
                Simd::storeu (Aptr, T1);
                /* A-B mod p */
                T2 = SimdExtra::sub_mod (A, B, P);
                /* multiply A-B by alpha and store it in Bptr */
                T3 = SimdExtra::mul_mod (T2, alpha, P, U);
                Simd::storeu (Bptr, T3);
            }

            /* Compute A[i]+B[i]*alpha[i], A[i]-B[i]*alpha[i] using Harvey's
             * algorithm, for 0 <= i < simd::vect_size.
             * Input must satisfy:
             *  - 0 <= A[i],B[i],alpha[i] < p
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < p
             */
            void
            Butterfly_DIT (Element& A, Element& B, const Element& alpha) const {
                Element tmp;
                this->fld->mul (tmp, alpha, B);
                this->fld->sub (B, A, tmp);
                this->fld->addin (A, tmp);
            }

            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::loadu (Aptr);
                B = Simd::loadu (Bptr);

                /* B*alpha mod P */
                T1 = SimdExtra::mul_mod (B, alpha, P, U);
                /* A+B*alpha */
                T2 = SimdExtra::add_mod (A, T1, P);
                Simd::storeu (Aptr, T2);
                /* A-B*alpha */
                T3 = SimdExtra::sub_mod (A, T1, P);
                Simd::storeu (Bptr, T3);
            }

            /******************************************************************/
            /* Utils **********************************************************/
            /******************************************************************/
            void
            init_powers (const Element & w) {
                /* compute w^i and set first subarray */
                this->fld->assign (pow_w[0], this->fld->one);
                for (size_t i = 1; i < n/2; i++)
                    this->fld->mul (pow_w[i], pow_w[i-1], w);

                /* Other elements can be set from previously computed values */
                size_t idx = n/2; /* index for next value to be written in */
                for (size_t k=2; k <= n/2; k<<=1)
                    for(size_t i = 0; i < n/2; i+=k, idx++)
                        pow_w[idx]  = pow_w[i];

                /* init powers in bitreverse order */
                size_t l = l2n, len = n >> 1, base_idx = 0;
                do
                {
                    l--;
                    for (size_t i = 0; i < len; i++)
                    {
                        size_t i_br = FFT_utils::bitreverse (i, l);
                        pow_w_br[base_idx + i] = pow_w[base_idx + i_br];
                    }
                    base_idx += len;
                    len >>= 1;
                }
                while (l) ;
            }
    };
}
#endif /* __LINBOX_fft_floating_INL */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
