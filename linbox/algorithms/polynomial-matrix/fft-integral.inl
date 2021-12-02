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


#ifndef __LINBOX_fft_integral_INL
#define __LINBOX_fft_integral_INL

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
                typename std::enable_if<Field::is_elt_integral_v
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
            const Residu_t p; /* p = field characteristic */
            const Residu_t p2; /* p2 = 2*p */

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
            typename Simd::aligned_vector pow_wp;
            typename Simd::aligned_vector pow_wp_br;

            /******************************************************************/
            /* constructor ****************************************************/
            /******************************************************************/
            FFT_base (const Field& F, size_t k, Element w)
                                            : fld(&F), l2n(k), n(1UL << l2n),
                                              p(F.characteristic()), p2(p << 1),
                                              pow_w(n-1), pow_w_br(n-1),
                                              pow_wp(n-1), pow_wp_br(n-1) {
                // TODO check size of p for integral
                init_powers (w);
            }

        public:
            void
            DIF (Element *coeffs) const {
                FFTSimdHelper<Simd::vect_size> h;
                /* w = n/2, f = 1 */
                DIF_core (coeffs, n >> 1, 1, pow_w.data(), pow_wp.data(), h);
                reduce_coeffs_2p (coeffs, h);
            }

            void
            DIT_reversed (Element *coeffs) const {
                FFTSimdHelper<Simd::vect_size> h;
                /* w = n/2, f = 1 */
                DIT_reversed_core (coeffs, n >> 1, 1, pow_w_br.data()+ (n-2),
                                                  pow_wp_br.data()+ (n-2), h);
                reduce_coeffs_4p (coeffs, h);
            }

            void
            DIT (Element *coeffs) const {
                FFTSimdHelper<Simd::vect_size> h;
                /* w = 1, f = n / 2 */
                DIT_core (coeffs, 1, n >> 1, pow_w.data() + (n-2),
                                             pow_wp.data() + (n-2), h);
                reduce_coeffs_4p (coeffs, h);
            }

            void
            DIF_reversed (Element *coeffs) const {
                FFTSimdHelper<Simd::vect_size> h;
                /* w = 1, f = n / 2 */
                DIF_reversed_core (coeffs, 1, n >> 1, pow_w_br.data(),
                                                      pow_wp_br.data(), h);
                reduce_coeffs_2p (coeffs, h);
            }

        protected:
            /******************************************************************/
            /* reduce *********************************************************/
            /******************************************************************/
            /* NoSimd */
            void
            reduce (Element &v, const Residu_t m) const {
                v -= (v >= m) ? m : 0;
            }

            void
            reduce_coeffs_2p (Element *coeffs, FFTSimdHelper<1>) const {
                for (size_t i = 0; i < n; i++)
                    reduce (coeffs[i], p);
            }

            void
            reduce_coeffs_4p (Element *coeffs, FFTSimdHelper<1>) const {
                for (size_t i = 0; i < n; i++) {
                    reduce (coeffs[i], p2);
                    reduce (coeffs[i], p);
                }
            }

            /* Simd */
            template <size_t VecSize>
            void
            reduce_coeffs_2p (Element *coeffs, FFTSimdHelper<VecSize>) const {
                if (n < Simd::vect_size) {
                    reduce_coeffs_2p (coeffs, FFTSimdHelper<1>());
                } else {
                    simd_vect_t P = Simd::set1 (p);
                    for (size_t i = 0; i < n; i += Simd::vect_size) {
                        simd_vect_t T = Simd::load (coeffs+i);
                        T = SimdExtra::reduce (T, P);
                        Simd::store (coeffs+i, T);
                    }
                }
            }

            template <size_t VecSize>
            void
            reduce_coeffs_4p (Element *coeffs, FFTSimdHelper<VecSize>) const {
                if (n < Simd::vect_size) {
                    reduce_coeffs_4p (coeffs, FFTSimdHelper<1>());
                } else {
                    simd_vect_t P = Simd::set1 (p);
                    simd_vect_t P2 = Simd::set1 (p2);
                    for (size_t i = 0; i < n; i += Simd::vect_size) {
                        simd_vect_t T = Simd::load (coeffs+i);
                        T = SimdExtra::reduce (T, P2);
                        T = SimdExtra::reduce (T, P);
                        Simd::store (coeffs+i, T);
                    }
                }
            }

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
            DIF_core (Element *coeffs, size_t w, size_t f, const Element *pow,
                        const Element *powp, FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                for ( ; w >= Simd::vect_size; pow += w, powp += w, f <<= 1,
                                                                   w >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w)
                        for (size_t j = 0; j < w; j += Simd::vect_size,
                                                    Aptr += Simd::vect_size,
                                                    Bptr += Simd::vect_size)
                            Butterfly_DIF (Aptr, Bptr, pow+j, powp+j, P, P2);
                }

                DIF_core_laststeps (coeffs, w, f, pow, powp, P, P2, h);
            }

            /* NoSimd */
            void
            DIF_core (Element *coeffs, size_t w, size_t f, const Element *pow,
                        const Element *powp, FFTSimdHelper<1>) const {
                for ( ; w > 0; pow += w, powp += w, f <<= 1, w >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w)
                        for (size_t j = 0; j < w; j++, Aptr++, Bptr++)
                            Butterfly_DIF (*Aptr, *Bptr, pow[j], powp[j], p2);
                }
            }

            /* DIT reversed ***************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                for ( ; w >= Simd::vect_size; f <<= 1, w >>= 1, pow -= f,
                                                                powp -= f) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        simd_vect_t alphap = Simd::set1 (powp[i]);
                        for (size_t j = 0; j < w; j += Simd::vect_size,
                                                    Aptr += Simd::vect_size,
                                                    Bptr += Simd::vect_size)
                            Butterfly_DIT (Aptr, Bptr, alpha, alphap, P, P2);
                    }
                }

                DIT_reversed_core_laststeps (coeffs, w, f, pow, powp, P, P2, h);
            }

            /* NoSimd */
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                FFTSimdHelper<1>) const {
                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f, powp -= f) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w) {
                        Element alpha = pow[i];
                        Element alphap = powp[i];
                        for (size_t j = 0; j < w; j++, Aptr++, Bptr++)
                            Butterfly_DIT (*Aptr, *Bptr, alpha, alphap, p2);
                    }
                }
            }

            /* DIT ************************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIT_core (Element *coeffs, size_t w, size_t f, const Element *pow,
                        const Element *powp, FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                DIT_core_firststeps (coeffs, w, f, pow, powp, P, P2, h);

                for ( ; w < n; w <<= 1, f >>= 1, pow -= w, powp -= w) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w)
                        for (size_t j = 0; j < w; j += Simd::vect_size,
                                                    Aptr += Simd::vect_size,
                                                    Bptr += Simd::vect_size)
                            Butterfly_DIT (Aptr, Bptr, pow+j, powp+j, P, P2);
                }
            }

            /* NoSimd */
            void
            DIT_core (Element *coeffs, size_t w, size_t f, const Element *pow,
                        const Element *powp, FFTSimdHelper<1>,
                        size_t bound = 0) const {
                bound = bound ? std::min (bound, n) : n;
                for ( ; w < bound; w <<= 1, f >>= 1, pow -= w, powp -= w) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w)
                        for (size_t j = 0; j < w; j++, Aptr++, Bptr++)
                            Butterfly_DIT (*Aptr, *Bptr, pow[j], powp[j], p2);
                }
            }

            /* DIF reversed ***************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                DIF_reversed_core_firststeps (coeffs, w, f, pow, powp, P, P2,h);

                for ( ; w < n; pow += f, powp += f, w <<= 1, f >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        simd_vect_t alphap = Simd::set1 (powp[i]);
                        for (size_t j = 0; j < w; j += Simd::vect_size,
                                                    Aptr += Simd::vect_size,
                                                    Bptr += Simd::vect_size)
                            Butterfly_DIF (Aptr, Bptr, alpha, alphap, P, P2);
                    }
                }
            }

            /* NoSimd */
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                FFTSimdHelper<1>, size_t bound = 0) const {
                bound = bound ? std::min (bound, n) : n;
                for ( ; w < bound; pow += f, powp += f, w <<= 1, f >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++, Aptr += w, Bptr += w) {
                        Element alpha = pow[i];
                        Element alphap = powp[i];
                        for (size_t j = 0; j < w; j++, Aptr++, Bptr++)
                            Butterfly_DIF (*Aptr, *Bptr, alpha, alphap, p2);
                    }
                }
            }

            /******************************************************************/
            /* Butterflies ****************************************************/
            /******************************************************************/
            /* Compute A[i]+B[i], (A[i]-B[i])*alpha[i] using Harvey's algorithm,
             * for 0 <= i < Simd::vect_size
             * Input must satisfy:
             *  - 0 <= A[i],B[i] < 2*p
             *  - 0 <= alpha[i] < p
             *  - p < 2^#nbits(Element) / 4
             *  - alphap[i] = Floor(alpha[i] * 2^#nbits(Element) / p)
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < 2*p
             *
             * Note: maybe 2^#nbits(Element) should be maxCardinality ? (in p<)
             */
            void
            Butterfly_DIF (Element& A, Element& B, const Element& alpha,
                                                const Element& alphap,
                                                const Element& p2) const {
                Element tmp = A;
                A += B;
                reduce (A, p2); /* A -= 2p if A >= 2p */
                B = tmp + (p2 - B);
                this->fld->mul_precomp_b_without_reduction (B, B, alpha, alphap);
            }

            void
            Butterfly_DIF (simd_vect_t& A, simd_vect_t& B,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
                simd_vect_t T1, T2, T3;

                /* A+B mod 2p */
                T1 = SimdExtra::add_mod (A, B, P2);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha */
                B = SimdExtra::mul_mod (T3, alpha, P, alphap);

                A = T1;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                            const Element *alpha_ptr, const Element *alphap_ptr,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3, T4, alpha, alphap;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);
                alpha = Simd::load (alpha_ptr);
                alphap = Simd::load (alphap_ptr);

                /* A+B mod 2p and store it in Aptr */
                T1 = SimdExtra::add_mod (A, B, P2);
                Simd::store (Aptr, T1);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha and store it in Bptr */
                T4 = SimdExtra::mul_mod (T3, alpha, P, alphap);
                Simd::store (Bptr, T4);
#else
                simd_vect_t alpha = Simd::load (alpha_ptr);
                simd_vect_t alphap = Simd::load (alphap_ptr);
                Butterfly_DIF (Aptr, Bptr, alpha, alphap, P, P2);
#endif
            }

            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3, T4;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* A+B mod 2p and store it in Aptr */
                T1 = SimdExtra::add_mod (A, B, P2);
                Simd::store (Aptr, T1);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha and store it in Bptr */
                T4 = SimdExtra::mul_mod (T3, alpha, P, alphap);
                Simd::store (Bptr, T4);
#else
                simd_vect_t A = Simd::load (Aptr);
                simd_vect_t B = Simd::load (Bptr);
                Butterfly_DIF (A, B, alpha, alphap, P, P2);
                Simd::store (Aptr, A);
                Simd::store (Bptr, B);
#endif
            }

            /* Compute A[i]+B[i]*alpha[i], A[i]-B[i]*alpha[i] using Harvey's
             * algorithm, for 0 <= i < simd::vect_size.
             * Input must satisfy:
             *  - 0 <= A[i],B[i] < 4*p
             *  - 0 <= alpha[i] < p
             *  - p < 2^#nbits(Element) / 4
             *  - alphap[i] = Floor(alpha[i] * 2^#nbits(Element) / p)
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < 4*p
             *
             * Note: maybe 2^#nbits(Element) should be maxCardinality ? (in p<)
             */
            void
            Butterfly_DIT (Element& A, Element& B, const Element& alpha,
                                                const Element& alphap,
                                                const Element& p2) const {
                reduce (A, p2); /* A -= 2p if A >= 2p */
                Element tmp;
                this->fld->mul_precomp_b_without_reduction (tmp, B, alpha, alphap);
                B = A + (p2 - tmp);
                A += tmp;
            }

            void
            Butterfly_DIT (simd_vect_t& A, simd_vect_t& B,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
                simd_vect_t T1, T2, T3;
                T1 = SimdExtra::reduce (A, P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdExtra::mul_mod (B, alpha, P, alphap);
                /* A+B*alpha */
                A = Simd::add (T1, T2);
                /* A-B*alpha (computed as A+(2p-B*alpha)) */
                T3 = Simd::sub (P2, T2);
                B = Simd::add (T1, T3);
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                            const Element *alpha_ptr, const Element *alphap_ptr,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3, alpha, alphap;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);
                alpha = Simd::load (alpha_ptr);
                alphap = Simd::load (alphap_ptr);

                T1 = SimdExtra::reduce (A, P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdExtra::mul_mod (B, alpha, P, alphap);
                /* A+B*alpha */
                A = Simd::add (T1, T2);
                Simd::store (Aptr, A);
                /* A-B*alpha (computed as A+(2p-B*alpha)) */
                T3 = Simd::sub (P2, T2);
                B = Simd::add (T1, T3);
                Simd::store (Bptr, B);
#else
                simd_vect_t alpha = Simd::load (alpha_ptr);
                simd_vect_t alphap = Simd::load (alphap_ptr);
                Butterfly_DIT (Aptr, Bptr, alpha, alphap, P, P2);
#endif
            }

            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                T1 = SimdExtra::reduce (A, P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdExtra::mul_mod (B, alpha, P, alphap);
                /* A+B*alpha */
                A = Simd::add (T1, T2);
                Simd::store (Aptr, A);
                /* A-B*alpha (computed as A+(2p-B*alpha)) */
                T3 = Simd::sub (P2, T2);
                B = Simd::add (T1, T3);
                Simd::store (Bptr, B);
#else
                simd_vect_t A = Simd::load (Aptr);
                simd_vect_t B = Simd::load (Bptr);
                Butterfly_DIT (A, B, alpha, alphap, P, P2);
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
                                const Element *pow, const Element *powp,
                                const simd_vect_t& P, const simd_vect_t& P2,
                                FFTSimdHelper<VecSize>) const {
                /* P and P2 are unused with no_simd fallback */
                DIF_core (coeffs, w, f, pow, powp, FFTSimdHelper<1>());
            }

            /* For vect_size == 2 */
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                const simd_vect_t& P, const simd_vect_t& P2,
                                FFTSimdHelper<2>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_core_laststeps (coeffs, w, f, pow, powp, P, P2,
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
                        T = SimdExtra::add_mod (V1, V2, P2);
                        V2 = SimdExtra::sub_mod (V1, V2, P2);

                        /* Result in T = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        Simd::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            /* For vect_size == 4 and Field != Modular<uint32_t, uint64_t> */
            template <typename E=Element, typename C=typename Field::Compute_t,
                      typename std::enable_if<!std::is_same<E, uint32_t>::value || !std::is_same<C, uint64_t>::value>::type* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                const simd_vect_t& P, const simd_vect_t& P2,
                                FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_core_laststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    simd_vect_t Wp = Simd::set(powp[0],powp[0],powp[1],powp[1]);

                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);
                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdExtra::add_mod (V1, V2, P2);
                        V2 = SimdExtra::sub_mod (V1, V2, P2);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        Simd::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
            /* For vect_size == 4 and Field == Modular<uint32_t, uint64_t> */
            template <typename E=Element, typename C=typename Field::Compute_t,
                      typename std::enable_if<std::is_same<E, uint32_t>::value && std::is_same<C, uint64_t>::value>::type* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                const simd_vect_t& P, const simd_vect_t& P2,
                                FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_core_laststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    simd_vect_t W = Simd::set1 (pow[1]);
                    simd_vect_t Wp = Simd::set1 (powp[1]);

                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, V3, V4, V5, V6, V7, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        T = SimdExtra::add_mod (V1, V2, P2);
                        V2 = SimdExtra::sub_mod (V1, V2, P2);

                        /* V4 = [D D H H] */
                        V4 = Simd::unpackhi (V2, V2);
                        /* Using extended mul (mulx) to compute V4*Wp as we only
                         * need to compute the product for half of the entries.
                         */
                        V7 = Simd128<uint64_t>::mulx (V4, Wp);
                        V5 = Simd::mullo (V7, P);
                        V6 = Simd::mullo (V4, W);
                        V4 = Simd::sub (V6, V5);
                        /* V3 = [ ? ? D H ] */
                        V3 = Simd::template shuffle<0xDD> (V4);
                        /* We need V1 = [A C E G], V2 = [B D F H] */
                        V1 = Simd::unpacklo (T, V2);
                        V2 = Simd::unpackhi (T, V3);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdExtra::add_mod (V1, V2, P2);
                        V2 = SimdExtra::sub_mod (V1, V2, P2);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        Simd::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }
#endif

            /* For vect_size == 8 and Field != Modular<uint32_t, uint64_t> */
            template <typename E=Element, typename C=typename Field::Compute_t,
                      typename std::enable_if<!std::is_same<E, uint32_t>::value || !std::is_same<C, uint64_t>::value>::type* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                const simd_vect_t& P, const simd_vect_t& P2,
                                FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_core_laststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1],
                                               pow[2], pow[2], pow[3], pow[3]);
                    simd_vect_t Wp = Simd::set (powp[0],powp[0],powp[1],powp[1],
                                               powp[2],powp[2],powp[3],powp[3]);
                    simd_vect_t W2 = Simd::set (pow[4], pow[4], pow[4], pow[4],
                                                pow[5], pow[5], pow[5], pow[5]);
                    simd_vect_t W2p = Simd::set(powp[4],powp[4],powp[4],powp[4],
                                               powp[5],powp[5],powp[5],powp[5]);

                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** step *********************************************/
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);
                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W2, W2p, P, P2);
                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdExtra::add_mod (V1, V2, P2);
                        V2 = SimdExtra::sub_mod (V1, V2, P2);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        Simd::unpacklohi (V1, V2, T, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
            /* For vect_size == 8 and Field == Modular<uint32_t, uint64_t> */
            template <typename E=Element, typename C=typename Field::Compute_t,
                      typename std::enable_if<std::is_same<E, uint32_t>::value && std::is_same<C, uint64_t>::value>::type* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                const Element *pow, const Element *powp,
                                const simd_vect_t& P, const simd_vect_t& P2,
                                FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_core_laststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    simd_vect_t W = Simd::set (pow[0], pow[1], pow[2], pow[3],
                                               pow[0], pow[1], pow[2], pow[3]);
                    simd_vect_t Wp = Simd::set (powp[0],powp[1],powp[2],powp[3],
                                               powp[0],powp[1],powp[2],powp[3]);
                    simd_vect_t W2 = Simd::set1 (pow[5]);
                    simd_vect_t W2p = Simd::set1 (powp[5]);

                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;
                        simd_vect_t V3, V4, V5, V6, V7, Q;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform into
                         *      V3 = [A B C D I J K L], V4 = [E F G H M N O P]
                         */
                        V3 = Simd256<uint64_t>::permute128<0x20> (V1, V2);
                        V4 = Simd256<uint64_t>::permute128<0x31> (V1, V2);

                        /*** step *********************************************/
                        Butterfly_DIF (V3, V4, W, Wp, P, P2);

                        /* transform into
                         *      V1 = [A E B F I M J N], V2 = [C G D H K O L P]
                         */
                        V1 = Simd::unpacklo_intrinsic (V3, V4);
                        V2 = Simd::unpackhi_intrinsic (V3, V4);

                        /*** last but one step ********************************/
                        T = SimdExtra::add_mod (V1, V2, P2);
                        V7 = SimdExtra::sub_mod (V1, V2, P2);

                        /* V4 = [D D H H L L P P ] */
                        V4 = Simd::unpackhi_intrinsic (V7, V7);
                        /* Using extended mul (mulx) to compute V4*Wp as we only
                         * need to compute the product for half of the entries.
                         */
                        Q = Simd256<uint64_t>::mulx (V4, W2p);
                        V5 = Simd::mullo (Q, P);
                        V6 = Simd::mullo (V4, W2);
                        V3 = Simd::sub (V6, V5);
                        /* V2 = [* * D H * * L P] */
                        V2 = Simd::template shuffle_twice<0xDD> (V3);
                        /* We need
                         *      V3 = [A C E G I K M O], V4 = [B D F H J L N P]
                         */
                        V1 = Simd::unpacklo_intrinsic (T, V7);
                        V2 = Simd::unpackhi_intrinsic (T, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdExtra::add_mod (V1, V2, P2);
                        V2 = SimdExtra::sub_mod (V1, V2, P2);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        Simd::unpacklohi (V1, V2, T, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }
#endif /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS */

            /* Fallback code, use no simd version of DIT_reversed_core */
            template<size_t VecSize>
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *pow, const Element *powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<VecSize>) const {
                /* P and P2 are unused with no_simd fallback */
                DIT_reversed_core (coeffs, w, f, pow, powp, FFTSimdHelper<1>());
            }

            /* For vect_size == 2 */
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *pow, const Element *powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<2>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_reversed_core_laststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, W, Wp;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::load (pow+i);
                        Wp = Simd::load (powp+i);

                        /* transform V1 = [A B], V2 = [C D]
                         *      into V1 = [A C], V2 = [B D]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step ****************************************/
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);

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
                                    const Element *pow, const Element *powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_reversed_core_laststeps (coeffs, w, f, pow, powp, P, P2,
                                                        FFTSimdHelper<1>());
                } else {
                    const Element *pow0 = pow - (f << 1);
                    const Element *powp0 = powp - (f << 1);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr,
                                                    pow0 += Simd::vect_size,
                                                    powp0 += Simd::vect_size) {
                        simd_vect_t V1, V2, W, Wp;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::set (pow[i], pow[i+1], pow[i], pow[i+1]);
                        Wp = Simd::set (powp[i], powp[i+1], powp[i], powp[i+1]);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);
                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step ****************************************/
                        W = Simd::load (pow0);
                        Wp = Simd::load (powp0);
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);

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
                                    const Element *pow, const Element *powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_reversed_core_laststeps (coeffs, w, f, pow, powp, P, P2,
                                                        FFTSimdHelper<1>());
                } else {
                    const Element *pow1 = pow - (f << 1);
                    const Element *powp1 = powp - (f << 1);
                    const Element *pow0 = pow1 - (f << 2);
                    const Element *powp0 = powp1 - (f << 2);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr,
                                                    pow0 += Simd::vect_size,
                                                    powp0 += Simd::vect_size) {
                        simd_vect_t V1, V2, W, Wp;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::set (pow[i], pow[i+1], pow[i], pow[i+1],
                                            pow[i], pow[i+1], pow[i], pow[i+1]);
                        Wp = Simd::set (powp[i], powp[i+1], powp[i], powp[i+1],
                                        powp[i], powp[i+1], powp[i], powp[i+1]);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** step *********************************************/
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);
                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        W = Simd::set (pow1[2*i], pow1[2*i+1], pow1[2*i+2],
                                                                    pow1[2*i+3],
                                       pow1[2*i], pow1[2*i+1], pow1[2*i+2],
                                                                pow1[2*i+3]);
                        Wp = Simd::set (powp1[2*i], powp1[2*i+1], powp1[2*i+2],
                                                                powp1[2*i+3],
                                        powp1[2*i], powp1[2*i+1], powp1[2*i+2],
                                                                powp1[2*i+3]);
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);
                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        Simd::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        W = Simd::load (pow0);
                        Wp = Simd::load (powp0);
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);

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
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<VecSize>) const {
                /* P and P2 are unused with no_simd fallback */
                DIT_core (coeffs, w, f, pow, powp, FFTSimdHelper<1>(),
                                                            Simd::vect_size);
                for ( ; w < Simd::vect_size; w <<= 1, f >>= 1, pow-=w, powp-=w);
            }

            /* For vect_size == 2 */
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<2>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_core_firststeps (coeffs, w, f, pow, powp, P, P2,
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
                        /* We know that entries of V1 and V2 are < P (because
                         * this is the first step), so the addition and the
                         * substration can be done with P, not P2.
                         */
                        T = Simd::add (V1, V2);
                        V2 = Simd::sub (P, V2);
                        V2 = Simd::add (V1, V2);

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
                    powp -= w;
                }
            }

            /* For vect_size == 4 and Field != Modular<uint32_t, uint64_t> */
            template <typename E=Element, typename C=typename Field::Compute_t,
                      typename std::enable_if<!std::is_same<E, uint32_t>::value || !std::is_same<C, uint64_t>::value>::type* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_core_firststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    f >>= 2;
                    w <<= 2;
                    pow -= 2;
                    powp -= 2;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    simd_vect_t Wp = Simd::set (powp[0],powp[0],powp[1],powp[1]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** first step (special butterfly with mul by 1) *****/
                        /* We know that entries of V1 and V2 are < P (because
                         * this is the first step), so the addition and the
                         * substration can be done with P, not P2.
                         */
                        T = Simd::add (V1, V2);
                        V2 = Simd::sub (P, V2);
                        V2 = Simd::add (V1, V2);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);

                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        Simd::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                    powp -= w;
                }
            }

#ifdef __FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS
            /* For vect_size == 4 and Field == Modular<uint32_t, uint64_t> */
            template <typename E=Element, typename C=typename Field::Compute_t,
                      typename std::enable_if<std::is_same<E, uint32_t>::value && std::is_same<C, uint64_t>::value>::type* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_core_firststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    f >>= 2;
                    w <<= 2;
                    pow -= 2;
                    powp -= 2;
                    simd_vect_t W = Simd::set1 (pow[1]);
                    simd_vect_t Wp = Simd::set1 (powp[1]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, V3, V4, T1, T2, T3, T4;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A E C G], V2 = [B F D H]
                         */
                        T1 = Simd::template shuffle<0xD8> (V1);
                        T2 = Simd::template shuffle<0xD8> (V2);
                        V1 = Simd::unpacklo (T1, T2);
                        V2 = Simd::unpackhi (T1, T2);

                        /*** first step (special butterfly with mul by 1) *****/
                        /* We know that entries of V1 and V2 are < P (because
                         * this is the first step), so the addition and the
                         * substration can be done with P, not P2.
                         */
                        V3 = Simd::add (V1, V2);
                        V4 = SimdExtra::sub_mod (V1, V2, P);

                        /*** second step **************************************/
                        /* T1 = [D D H H] */
                        T1 = Simd::unpackhi (V4, V4);
                        /* Using extended mul (mulx) to compute T1*Wp as we only
                         * need to compute the product for half of the entries.
                         */
                        T2 = Simd128<uint64_t>::mulx (T1, Wp);
                        T3 = Simd::mullo (T2, P);
                        T4 = Simd::mullo (T1, W);
                        T1 = Simd::sub (T4, T3);
                        /* T2 = [ ? ? D H ] */
                        T2 = Simd::template shuffle<0xDD> (T1);
                        /* We need V1 = [A B E F], V2 = [C D G H] */
                        V1 = Simd::unpacklo (V3, V4);
                        V2 = Simd::unpackhi (V3, T2);

                        T1 = Simd::add (V1, V2);
                        T2 = SimdExtra::sub_mod (V1, V2, P2);

                        /* transform T1 = [A B E F], T2 = [C D G H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        V1 = Simd::unpacklo (T1, T2);
                        V2 = Simd::unpackhi (T1, T2);
                        T1 = Simd::template shuffle<0xD8> (V1);
                        T2 = Simd::template shuffle<0xD8> (V2);

                        Simd::store (coeffs, T1);
                        Simd::store (coeffs + Simd::vect_size, T2);
                    }
                    pow -= w;
                    powp -= w;
                }
            }
#endif

            /* For vect_size == 8 and Field != Modular<uint32_t, uint64_t> */
            template <typename E=Element, typename C=typename Field::Compute_t,
                      typename std::enable_if<!std::is_same<E, uint32_t>::value || !std::is_same<C, uint64_t>::value>::type* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_core_firststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    f >>= 3;
                    w <<= 3;
                    pow -= 2;
                    powp -= 2;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[0], pow[0],
                                               pow[1], pow[1], pow[1], pow[1]);
                    simd_vect_t Wp = Simd::set (powp[0],powp[0],powp[0],powp[0],
                                                powp[1],powp[1],powp[1],powp[1]);
                    pow -= 4;
                    powp -= 4;
                    simd_vect_t W2 = Simd::set (pow[0], pow[0], pow[1], pow[1],
                                                pow[2], pow[2], pow[3], pow[3]);
                    simd_vect_t W2p = Simd::set(powp[0],powp[0],powp[1],powp[1],
                                               powp[2],powp[2],powp[3],powp[3]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** first step (special butterfly with mul by 1) *****/
                        /* We know that entries of V1 and V2 are < P (because
                         * this is the first step), so the addition and the
                         * substration can be done with P, not P2.
                         */
                        T = Simd::add (V1, V2);
                        V2 = Simd::sub (P, V2);
                        V2 = Simd::add (V1, V2);

                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        Simd::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** third step ***************************************/
                        Butterfly_DIT (V1, V2, W2, W2p, P, P2);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                    powp -= w;
                }
            }

#ifdef __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS
            /* For vect_size == 8 and Field == Modular<uint32_t, uint64_t> */
            template <typename E=Element, typename C=typename Field::Compute_t,
                      typename std::enable_if<std::is_same<E, uint32_t>::value && std::is_same<C, uint64_t>::value>::type* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIT_core_firststeps (coeffs, w, f, pow, powp, P, P2,
                                                            FFTSimdHelper<1>());
                } else {
                    f >>= 3;
                    w <<= 3;
                    pow -= 2;
                    powp -= 2;
                    simd_vect_t W = Simd::set1 (pow[1]);
                    simd_vect_t Wp = Simd::set1 (powp[1]);
                    pow -= 4;
                    powp -= 4;
                    simd_vect_t W2 = Simd::set (pow[0], pow[1], pow[2], pow[3],
                                                pow[0], pow[1], pow[2], pow[3]);
                    simd_vect_t W2p = Simd::set(powp[0],powp[1],powp[2],powp[3],
                                               powp[0],powp[1],powp[2],powp[3]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, V3, V4, V5, V6, V7, Q;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);

                        /* transform into
                         *      V3 = [A I C K E M G O], V4 = [B J D L F N H P]
                         */
                        V6 = Simd::unpacklo_intrinsic(V1,V2);
                        V7 = Simd::unpackhi_intrinsic(V1,V2);
                        V3 = Simd256<uint64_t>::unpacklo_intrinsic(V6,V7);
                        V4 = Simd256<uint64_t>::unpackhi_intrinsic(V6,V7);

                        /*** first step (special butterfly with mul by 1) *****/
                        /* We know that entries of V1 and V2 are < P (because
                         * this is the first step), so the addition and the
                         * substration can be done with P, not P2.
                         */
                        V1 = Simd::add(V3,V4);
                        V2 = SimdExtra::sub_mod (V3, V4, P);

                        /*** second step **************************************/
                        /* V5 = [D D L L H H P P] */
                        V5 = Simd::unpackhi_intrinsic (V2, V2);
                        /* Using extended mul (mulx) to compute V5*Wp as we only
                         * need to compute the product for half of the entries.
                         */
                        Q = Simd256<uint64_t>::mulx (V5, Wp);
                        V6 = Simd::mullo (Q, P);
                        V7 = Simd::mullo (V5, W);
                        V3 = Simd::sub (V7, V6);
                        /* V7 = [D L * * H P * *] */
                        V7 = Simd::template shuffle_twice<0xFD> (V3);
                        /* We need
                         *      V3 = [A B I J E F M N], V4 = [C D K L G H O P]
                         */
                        V6 = Simd256<uint64_t>::unpacklo_intrinsic (V2, V7);
                        V3 = Simd::unpacklo_intrinsic (V1, V6);
                        V4 = Simd::unpackhi_intrinsic (V1, V6);

                        V1 = Simd::add (V3, V4);
                        V2 = SimdExtra::sub_mod (V3, V4, P2);

                        /* transform into
                         *      V3 = [A B C D I J K L], V4 = [E F G H M N O P]
                         */
                        V6 = Simd256<uint64_t>::unpacklo_intrinsic (V1, V2);
                        V7 = Simd256<uint64_t>::unpackhi_intrinsic (V1, V2);
                        V3 = Simd256<uint64_t>::permute128<0x20> (V6, V7);
                        V4 = Simd256<uint64_t>::permute128<0x31> (V6, V7);

                        /*** third step ***************************************/
                        Butterfly_DIT (V3, V4, W2, W2p, P, P2);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        V1 = Simd256<uint64_t>::permute128<0x20> (V3, V4);
                        V2 = Simd256<uint64_t>::permute128<0x31> (V3, V4);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                    powp -= w;
                }
            }
#endif /* __FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS */

            /* Fallback code, use no simd version of DIF_reversed_core */
            template<size_t VecSize>
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<VecSize>) const {
                /* P and P2 are unused with no_simd fallback */
                DIF_reversed_core (coeffs, w, f, pow, powp, FFTSimdHelper<1>(),
                                                            Simd::vect_size);
                for ( ; w < Simd::vect_size; pow+=f, powp+=f, w <<= 1, f >>= 1);
            }

            /* For vect_size == 2 */
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<2>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_reversed_core_firststeps (coeffs, w, f, pow, powp, P,
                                                        P2, FFTSimdHelper<1>());
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, W, Wp;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::load (pow+i);
                        Wp = Simd::load (powp+i);

                        /* transform V1 = [A B], V2 = [C D]
                         *      into V1 = [A C], V2 = [B D]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** first step ***************************************/
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);

                        /* Result in V1 = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        Simd::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow += f;
                    powp += f;
                    w <<= 1;
                    f >>= 1;
                }
            }

            /* For vect_size == 4 */
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<4>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_reversed_core_firststeps (coeffs, w, f, pow, powp, P,
                                                        P2, FFTSimdHelper<1>());
                } else {
                    const Element *pow0 = pow;
                    const Element *powp0 = powp;
                    pow += f;
                    powp += f;
                    f >>= 1;
                    w <<= 2;
                    for (size_t i = 0; i < f; i += 2, coeffs += incr,
                                                    pow0 += Simd::vect_size,
                                                    powp0 += Simd::vect_size) {
                        simd_vect_t V1, V2, W, Wp;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::load (pow0);
                        Wp = Simd::load (powp0);

                        /* transform V1 = [A B C D], V2 = [E F G H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        Simd::pack (V1, V2, V1, V2);


                        /*** first step ***************************************/
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** second step **************************************/
                        W = Simd::set (pow[i], pow[i+1], pow[i], pow[i+1]);
                        Wp = Simd::set (powp[i], powp[i+1], powp[i], powp[i+1]);
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);

                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        Simd::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow += f;
                    powp += f;
                    f >>= 1;
                }
            }

            /* For vect_size == 8 */
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P, const simd_vect_t& P2,
                                    FFTSimdHelper<8>) const {
                const constexpr size_t incr = Simd::vect_size << 1;
                if (n < incr) {
                    DIF_reversed_core_firststeps (coeffs, w, f, pow, powp, P,
                                                        P2, FFTSimdHelper<1>());
                } else {
                    const Element *pow0 = pow;
                    const Element *powp0 = powp;
                    pow += f;
                    powp += f;
                    f >>= 1;
                    const Element *pow1 = pow;
                    const Element *powp1 = powp;
                    pow += f;
                    powp += f;
                    f >>= 1;
                    w <<= 3;
                    for (size_t i = 0; i < f; i += 2, coeffs += incr,
                                                    pow0 += Simd::vect_size,
                                                    powp0 += Simd::vect_size) {
                        simd_vect_t V1, V2, W, Wp;

                        V1 = Simd::load (coeffs);
                        V2 = Simd::load (coeffs + Simd::vect_size);
                        W = Simd::load (pow0);
                        Wp = Simd::load (powp0);

                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** first step ***************************************/
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);

                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** second step **************************************/
                        W = Simd::set (pow1[2*i], pow1[2*i+1], pow1[2*i+2],
                                                                    pow1[2*i+3],
                                       pow1[2*i], pow1[2*i+1], pow1[2*i+2],
                                                                pow1[2*i+3]);
                        Wp = Simd::set (powp1[2*i], powp1[2*i+1], powp1[2*i+2],
                                                                powp1[2*i+3],
                                       powp1[2*i], powp1[2*i+1], powp1[2*i+2],
                                                                powp1[2*i+3]);
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        /*** third step ***************************************/
                        W = Simd::set (pow[i], pow[i+1], pow[i], pow[i+1],
                                            pow[i], pow[i+1], pow[i], pow[i+1]);
                        Wp = Simd::set (powp[i], powp[i+1], powp[i], powp[i+1],
                                        powp[i], powp[i+1], powp[i], powp[i+1]);
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        Simd::pack (V1, V2, V1, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow += f;
                    powp += f;
                    f >>= 1;
                }
            }


            /******************************************************************/
            /* Utils **********************************************************/
            /******************************************************************/
            void
            init_powers (const Element & w) {
                typename Field::Compute_t t;
                /* compute w^i and set first subarray */
                this->fld->assign (pow_w[0], this->fld->one);
                this->fld->precomp_b (t, pow_w[0]);
                pow_wp[0] = static_cast<Element> (t);
                for (size_t i = 1; i < n/2; i++) {
                    this->fld->mul (pow_w[i], pow_w[i-1], w);
                    this->fld->precomp_b (t, pow_w[i]);
                    pow_wp[i] = static_cast<Element> (t);
                }

                /* Other elements can be set from previously computed values */
                size_t idx = n/2; /* index for next value to be written in */
                for (size_t k=2; k <= n/2; k<<=1)
                    for(size_t i = 0; i < n/2; i+=k, idx++) {
                        pow_w[idx] = pow_w[i];
                        pow_wp[idx] = pow_wp[i];
                    }

                /* init powers in bitreverse order */
                size_t l = l2n, len = n >> 1, base_idx = 0;
                do
                {
                    l--;
                    for (size_t i = 0; i < len; i++)
                    {
                        size_t i_br = FFT_utils::bitreverse (i, l);
                        pow_w_br[base_idx + i] = pow_w[base_idx + i_br];
                        pow_wp_br[base_idx + i] = pow_wp[base_idx + i_br];
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
                typename std::enable_if<Field::is_elt_integral_v
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
            const Residu_t p; /* p = field characteristic */
            const Residu_t p2; /* p2 = 2*p */

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
            typename Simd::aligned_vector pow_wp;
            typename Simd::aligned_vector pow_wp_br;

            /******************************************************************/
            /* constructor ****************************************************/
            /******************************************************************/
            FFT_multi_base (const Field& F, size_t k, Element w)
                                            : fld(&F), l2n(k), n(1UL << l2n),
                                              p(F.characteristic()), p2(p << 1),
                                              pow_w(n-1), pow_w_br(n-1),
                                              pow_wp(n-1), pow_wp_br(n-1) {
                // TODO check size of p for integral
                init_powers (w);
            }

        public:
            void
            DIF (Element *coeffs, size_t stride) const {
                FFTSimdHelper<Simd::vect_size> h;
                /* w = n/2, f = 1 */
                DIF_core (coeffs, n >> 1, 1, stride, pow_w.data(),
                                                     pow_wp.data(), h);
                reduce_coeffs_2p (coeffs, stride, h);
            }

            void
            DIT_reversed (Element *coeffs, size_t stride) const {
                FFTSimdHelper<Simd::vect_size> h;
                /* w = n/2, f = 1 */
                DIT_reversed_core (coeffs, n >> 1, 1, stride,
                                                    pow_w_br.data()+ (n-2),
                                                    pow_wp_br.data()+ (n-2), h);
                reduce_coeffs_4p (coeffs, stride, h);
            }

            void
            DIT (Element *coeffs, size_t stride) const {
                FFTSimdHelper<Simd::vect_size> h;
                /* w = 1, f = n / 2 */
                DIT_core (coeffs, 1, n >> 1, stride, pow_w.data() + (n-2),
                                                     pow_wp.data() + (n-2), h);
                reduce_coeffs_4p (coeffs, stride, h);
            }

            void
            DIF_reversed (Element *coeffs, size_t stride) const {
                FFTSimdHelper<Simd::vect_size> h;
                /* w = 1, f = n / 2 */
                DIF_reversed_core (coeffs, 1, n >> 1, stride, pow_w_br.data(),
                                                        pow_wp_br.data(), h);
                reduce_coeffs_2p (coeffs, stride, h);
            }

       protected:
            /******************************************************************/
            /* reduce *********************************************************/
            /******************************************************************/
            /* NoSimd */
            void
            reduce (Element &v, const Residu_t m) const {
                v -= (v >= m) ? m : 0;
            }

            void
            reduce_coeffs_2p (Element *coeffs, size_t stride,
                                                    FFTSimdHelper<1>) const {
                for (size_t i = 0; i < n*stride; i++)
                    reduce (coeffs[i], p);
            }

            void
            reduce_coeffs_4p (Element *coeffs, size_t stride,
                                                    FFTSimdHelper<1>) const {
                for (size_t i = 0; i < n*stride; i++) {
                    reduce (coeffs[i], p2);
                    reduce (coeffs[i], p);
                }
            }

            /* Simd */
            template <size_t VecSize>
            void
            reduce_coeffs_2p (Element *coeffs, size_t stride,
                                                FFTSimdHelper<VecSize>) const {
                simd_vect_t P = Simd::set1 (p);
                size_t i = 0;
                for ( ; i+Simd::vect_size <= n*stride ; i += Simd::vect_size) {
                    simd_vect_t T = Simd::loadu (coeffs+i);
                    T = SimdExtra::reduce (T, P);
                    Simd::storeu (coeffs+i, T);
                }
                for ( ; i < n*stride; i++)
                    reduce (coeffs[i], p);
            }

            template <size_t VecSize>
            void
            reduce_coeffs_4p (Element *coeffs, size_t stride,
                                                FFTSimdHelper<VecSize>) const {
                simd_vect_t P = Simd::set1 (p);
                simd_vect_t P2 = Simd::set1 (p2);
                size_t i = 0;
                for ( ; i+Simd::vect_size <= n*stride ; i += Simd::vect_size) {
                    simd_vect_t T = Simd::loadu (coeffs+i);
                    T = SimdExtra::reduce (T, P2);
                    T = SimdExtra::reduce (T, P);
                    Simd::storeu (coeffs+i, T);
                }
                for ( ; i < n*stride; i++) {
                    reduce (coeffs[i], p2);
                    reduce (coeffs[i], p);
                }
            }

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
                                        const Element *pow, const Element *powp,
                                        FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                for ( ; w > 0; pow += w, powp += w, f <<= 1, w >>= 1) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws)
                        for (size_t j = 0; j < w; j += 1) {
                            simd_vect_t alpha = Simd::set1 (pow[j]);
                            simd_vect_t alphap = Simd::set1 (powp[j]);
                            size_t l = 0;
                            for ( ; l + Simd::vect_size <= stride;
                                                        l += Simd::vect_size,
                                                        Aptr += Simd::vect_size,
                                                        Bptr += Simd::vect_size)
                                Butterfly_DIF (Aptr, Bptr, alpha, alphap, P,P2);
                            for ( ; l < stride; l++, Aptr++, Bptr++)
                                Butterfly_DIF (*Aptr, *Bptr, pow[j],powp[j],p2);
                        }
                }
            }

            /* NoSimd */
            void
            DIF_core (Element *coeffs, size_t w, size_t f, size_t stride,
                                        const Element *pow, const Element *powp,
                                        FFTSimdHelper<1>) const {
                for ( ; w > 0; pow += w, powp += w, f <<= 1, w >>= 1) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws)
                        for (size_t j = 0; j < w; j += 1)
                            for (size_t l = 0; l < stride; l++, Aptr++, Bptr++)
                                Butterfly_DIF (*Aptr, *Bptr, pow[j],powp[j],p2);
                }
            }

            /* DIT reversed ***************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                        size_t stride,
                                        const Element *pow, const Element *powp,
                                        FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f, powp -= f) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        simd_vect_t alphap = Simd::set1 (powp[i]);
                        size_t j = 0;
                        for ( ; j + Simd::vect_size <= ws; j += Simd::vect_size,
                                                        Aptr += Simd::vect_size,
                                                        Bptr += Simd::vect_size)
                            Butterfly_DIT (Aptr, Bptr, alpha, alphap, P, P2);
                        for ( ; j < ws; j++, Aptr++, Bptr++)
                            Butterfly_DIT (*Aptr, *Bptr, pow[i], powp[i], p2);
                    }
                }
            }

            /* NoSimd */
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                        size_t stride,
                                        const Element *pow, const Element *powp,
                                        FFTSimdHelper<1>) const {
                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f, powp -= f) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws) {
                        Element alpha = pow[i];
                        Element alphap = powp[i];
                        for (size_t j = 0; j < ws; j += 1, Aptr++, Bptr++)
                            Butterfly_DIT (*Aptr, *Bptr, alpha, alphap, p2);
                    }
                }
            }

            /* DIT ************************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIT_core (Element *coeffs, size_t w, size_t f, size_t stride,
                                        const Element *pow, const Element *powp,
                                        FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                for ( ; w < n; w <<= 1, f >>= 1, pow -= w, powp -= w) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws)
                        for (size_t j = 0; j < w; j += 1) {
                            simd_vect_t alpha = Simd::set1 (pow[j]);
                            simd_vect_t alphap = Simd::set1 (powp[j]);
                            size_t l = 0;
                            for ( ; l + Simd::vect_size <= stride;
                                                        l += Simd::vect_size,
                                                        Aptr += Simd::vect_size,
                                                        Bptr += Simd::vect_size)
                                Butterfly_DIT (Aptr, Bptr, alpha, alphap, P,P2);
                            for ( ; l < stride; l++, Aptr++, Bptr++)
                                Butterfly_DIT (*Aptr, *Bptr, pow[j],powp[j],p2);
                        }
                }
            }

            /* NoSimd */
            void
            DIT_core (Element *coeffs, size_t w, size_t f, size_t stride,
                                        const Element *pow, const Element *powp,
                                        FFTSimdHelper<1>) const {
                for ( ; w < n; w <<= 1, f >>= 1, pow -= w, powp -= w) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws)
                        for (size_t j = 0; j < w; j += 1)
                            for (size_t l = 0; l < stride; l++, Aptr++, Bptr++)
                                Butterfly_DIT (*Aptr, *Bptr, pow[j],powp[j],p2);
                }
            }

            /* DIF reversed ***************************************************/
            /* Simd */
            template<size_t VecSize>
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                        size_t stride,
                                        const Element *pow, const Element *powp,
                                        FFTSimdHelper<VecSize> h) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                for ( ; w < n; pow += f, powp += f, w <<= 1, f >>= 1) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        simd_vect_t alphap = Simd::set1 (powp[i]);
                        size_t j = 0;
                        for ( ; j + Simd::vect_size <= ws; j += Simd::vect_size,
                                                        Aptr += Simd::vect_size,
                                                        Bptr += Simd::vect_size)
                            Butterfly_DIF (Aptr, Bptr, alpha, alphap, P, P2);
                        for ( ; j < ws; j++, Aptr++, Bptr++)
                            Butterfly_DIF (*Aptr, *Bptr, pow[i], powp[i], p2);
                    }
                }
            }

            /* NoSimd */
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                        size_t stride,
                                        const Element *pow, const Element *powp,
                                        FFTSimdHelper<1>) const {
                for ( ; w < n; pow += f, powp += f, w <<= 1, f >>= 1) {
                    size_t ws = w*stride;
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + ws;
                    for (size_t i = 0; i < f; i++, Aptr += ws, Bptr += ws) {
                        Element alpha = pow[i];
                        Element alphap = powp[i];
                        for (size_t j = 0; j < ws; j += 1, Aptr++, Bptr++)
                            Butterfly_DIF (*Aptr, *Bptr, alpha, alphap, p2);
                    }
                }
            }

            /******************************************************************/
            /* Butterflies ****************************************************/
            /******************************************************************/
            /* Compute A[i]+B[i], (A[i]-B[i])*alpha[i] using Harvey's algorithm,
             * for 0 <= i < Simd::vect_size
             * Input must satisfy:
             *  - 0 <= A[i],B[i] < 2*p
             *  - 0 <= alpha[i] < p
             *  - p < 2^#nbits(Element) / 4
             *  - alphap[i] = Floor(alpha[i] * 2^#nbits(Element) / p)
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < 2*p
             *
             * Note: maybe 2^#nbits(Element) should be maxCardinality ? (in p<)
             */
            void
            Butterfly_DIF (Element& A, Element& B, const Element& alpha,
                                                const Element& alphap,
                                                const Element& p2) const {
                Element tmp = A;
                A += B;
                reduce (A, p2); /* A -= 2p if A >= 2p */
                B = tmp + (p2 - B);
                this->fld->mul_precomp_b_without_reduction (B, B, alpha, alphap);
            }

            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
                simd_vect_t A, B, T1, T2, T3, T4;
                A = Simd::loadu (Aptr);
                B = Simd::loadu (Bptr);

                /* A+B mod 2p and store it in Aptr */
                T1 = SimdExtra::add_mod (A, B, P2);
                Simd::storeu (Aptr, T1);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha and store it in Bptr */
                T4 = SimdExtra::mul_mod (T3, alpha, P, alphap);
                Simd::storeu (Bptr, T4);
            }

            /* Compute A[i]+B[i]*alpha[i], A[i]-B[i]*alpha[i] using Harvey's
             * algorithm, for 0 <= i < simd::vect_size.
             * Input must satisfy:
             *  - 0 <= A[i],B[i] < 4*p
             *  - 0 <= alpha[i] < p
             *  - p < 2^#nbits(Element) / 4
             *  - alphap[i] = Floor(alpha[i] * 2^#nbits(Element) / p)
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < 4*p
             *
             * Note: maybe 2^#nbits(Element) should be maxCardinality ? (in p<)
             */
            void
            Butterfly_DIT (Element& A, Element& B, const Element& alpha,
                                                const Element& alphap,
                                                const Element& p2) const {
                reduce (A, p2); /* A -= 2p if A >= 2p */
                Element tmp;
                this->fld->mul_precomp_b_without_reduction (tmp, B, alpha, alphap);
                B = A + (p2 - tmp);
                A += tmp;
            }

            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::loadu (Aptr);
                B = Simd::loadu (Bptr);

                T1 = SimdExtra::reduce (A, P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdExtra::mul_mod (B, alpha, P, alphap);
                /* A+B*alpha */
                A = Simd::add (T1, T2);
                Simd::storeu (Aptr, A);
                /* A-B*alpha (computed as A+(2p-B*alpha)) */
                T3 = Simd::sub (P2, T2);
                B = Simd::add (T1, T3);
                Simd::storeu (Bptr, B);
            }

            /******************************************************************/
            /* Utils **********************************************************/
            /******************************************************************/
            void
            init_powers (const Element & w) {
                typename Field::Compute_t t;
                /* compute w^i and set first subarray */
                this->fld->assign (pow_w[0], this->fld->one);
                this->fld->precomp_b (t, pow_w[0]);
                pow_wp[0] = static_cast<Element> (t);
                for (size_t i = 1; i < n/2; i++) {
                    this->fld->mul (pow_w[i], pow_w[i-1], w);
                    this->fld->precomp_b (t, pow_w[i]);
                    pow_wp[i] = static_cast<Element> (t);
                }

                /* Other elements can be set from previously computed values */
                size_t idx = n/2; /* index for next value to be written in */
                for (size_t k=2; k <= n/2; k<<=1)
                    for(size_t i = 0; i < n/2; i+=k, idx++) {
                        pow_w[idx] = pow_w[i];
                        pow_wp[idx] = pow_wp[i];
                    }

                /* init powers in bitreverse order */
                size_t l = l2n, len = n >> 1, base_idx = 0;
                do
                {
                    l--;
                    for (size_t i = 0; i < len; i++)
                    {
                        size_t i_br = FFT_utils::bitreverse (i, l);
                        pow_w_br[base_idx + i] = pow_w[base_idx + i_br];
                        pow_wp_br[base_idx + i] = pow_wp[base_idx + i_br];
                    }
                    base_idx += len;
                    len >>= 1;
                }
                while (l) ;
            }


    };
}
#endif /* __LINBOX_fft_integral_INL */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
