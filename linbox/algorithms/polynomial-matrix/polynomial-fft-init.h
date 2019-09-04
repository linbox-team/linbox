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


#ifndef __LINBOX_polynomial_fft_init_H
#define __LINBOX_polynomial_fft_init_H


#include <iostream>
#include <type_traits>
#include "linbox/linbox-config.h"
#include "linbox/util/error.h"
#include <givaro/givpower.h>
#include <fflas-ffpack/fflas/fflas_simd.h>
#include <fflas-ffpack/utils/align-allocator.h>

#include "linbox/algorithms/polynomial-matrix/simd-additional-functions.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-utils.h"

//#define _FFT_NO_INLINED_BUTTERFLY 1

namespace LinBox {

    // Should these function be static method of the FFT class ???
    namespace FFT_internals {
        /* return i with its first n bits reversed */
        static inline size_t
        bitreverse (size_t i, size_t n)
        {
            size_t r = 0;
            for (size_t j = 0; j < n; j++, i>>=1)
            {
                r <<= 1;
                r |= i & 1;
            }
            return r;
        }

        /* return e-p if e >= p, return p otherwise */
        template<typename Element, typename Residu>
        static inline void
        sub_if_greater (Element &e, const Residu p)
        {
            e -= (e >= p) ? p : 0;
        }

        /* Return true if w is a primitive 2^k-root of the multiplicative group
         * of the field, return false otherwise.
         */
        template <typename Field>
        bool
        is_primitive_root (const Field& fld, typename Field::Element w,
                                                                    size_t k) {
            for (size_t i = 0; i < k; i++)
            {
                if (fld.isOne (w))
                    return false;
                else
                    fld.mul (w, w, w);
            }
            return true;
        }

        /* Return a primitive 2^k-root of the multiplicative group of the field.
         */
        template <typename Field>
        typename Field::Element
        compute_primitive_root (const Field& fld, size_t k)
        {
            /* write p-1 as 2^val2p * m where m is odd */
            size_t val2p = 0;
            typename Field::Residu_t m = fld.characteristic() - 1;
            for ( ; (m & 1) == 0 ; m >>= 1, val2p++) ;

            if (k > val2p)
                throw LinBoxError ("FFT: k is larger that log2(p-1)");

            typename Field::Element t, w;
            for (t = 2; ; t++) {
                Givaro::dom_power (w, t, m, fld); /* w <- t^m */
                for (size_t i = 0; i < val2p-k; i++)
                    fld.mul (w, w, w);
                /* now w = t^(2^(val2p-k)*m) */
                if (is_primitive_root (fld, w, k))
                    return w;
            }
        }

    }

    // TODO check size of p for integral

    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
    template<typename Field, typename Simd = Simd<typename Field::Element>,
                             typename Enable = void>
    class FFT;

    template <typename Field, typename Simd>
    class FFT<Field, Simd,
                FFT_utils::enable_if_floating_compatible_simd_t<Field, Simd>> {
        /* Types **************************************************************/
        private:
            using Element = typename Field::Element;
            using simd_vect_t = typename std::conditional<std::is_same<Simd, NoSimd<typename Field::Element>>::value, Element, typename Simd::vect_t>::type;
            using elt_vect_t = typename Simd::aligned_vector;
            using SimdMod = SimdModular<Field, Simd>;
            using SimdExtra = FFT_utils::SimdExtra<Field, Simd>;

        public:
            /* constructor ****************************************************/
            FFT (const Field& F, size_t k, Element w = 0)
                                            : fld(&F), l2n(k), n(1UL << l2n),
                                              pow_w(n-1) , pow_w_br(n-1) {
                if (!k)
                    throw LinBoxError ("FFT: k must be positive");

                if (w == 0)
                    w = compute_primitive_root ();
                else if (!is_primitive_root (w))
                    throw LinBoxError ("FFT: w is not a primitive k-root of unity");
                init_powers (w);
            }

            /* getters */
            const Field &
            field () const {
                return *fld;
            }

            const Element &
            getRoot() const
            {
                if (l2n == 1) /* if n=2^ln equals 2, return -1 */
                    return fld->mOne;
                else
                    return pow_w[1];
            }


            /******************************************************************/
            /* main functions: FFT direct and inverse *************************/
            /******************************************************************/

            /* Perform a FFT in place on the array 'coeffs' of size n.
             * Input:
             *  - must be < p
             *  - is read in natural order
             * Output:
             *  - are < p
             *  - is written in bitreversed order.
             */
            void
            FFT_direct (Element *coeffs) const {
                DIF (coeffs); /* or DIT_reversed */
            }

            /* Perform an inverse FFT in place on the array 'coeffs' of size n.
             * Input:
             *  - must be < p
             *  - is read in bitreversed order.
             * Output:
             *  - are < p
             *  - is written in natural order
             */
            void
            FFT_inverse (Element *coeffs) const {
                DIT (coeffs); /* or DIF_reversed */
            }

        protected:
            /******************************************************************/
            /* Attributes *****************************************************/
            /******************************************************************/
            const Field *fld;
            size_t l2n; /* log2 of size */
        //uint64_t n; /* 2^l2n */
        size_t n;
            elt_vect_t pow_w; /* Table of roots of unity.
                               * If w = primitive n-th root, then the table is:
                               * 1, w, w^2, ..., w^{n/2-1},
                               * 1, w^2, w^4, ..., w^{n/2-2},
                               * 1, w^4, w^8, ..., w^{n/2-4}
                               * ...
                               * 1, w^{n/8}, w^{n/4}, w^{3n/8},
                               * 1, w^{n/4},
                               * 1.
                               * Its size is n-1.
                               */
            elt_vect_t pow_w_br; /* Same as above with each subarray in
                                  * bitreverse order.
                                  */

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
            void
            DIF (Element *coeffs) const {
                /* w = n/2, f = 1 */
                DIF_core (coeffs, n >> 1, 1, pow_w.data());
            }

            void
            DIF_core_no_simd (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                for ( ; w > 0; pow += w, f <<= 1, w >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++)
                        for (size_t j = 0; j < w; j += 1)
                            Butterfly_DIF (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                             pow[j]);
                }
            }

            /* NoSimd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_core (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                DIF_core_no_simd (coeffs, w, f, pow);
            }

            /* Simd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIF_core (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                for ( ; w >= Simd::vect_size; pow += w, f <<= 1, w >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++)
                        for (size_t j = 0; j < w; j += Simd::vect_size)
                            Butterfly_DIF (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                                            pow+j, P, U);
                }

                DIF_core_laststeps (coeffs, w, f, pow, P, U);
            }

            /* DIT reversed ***************************************************/
            void
            DIT_reversed (Element *coeffs) const {
                /* w = n/2, f = 1 */
                DIT_reversed_core (coeffs, n >> 1, 1, pow_w_br.data()+ (n-2));
            }

            void
            DIT_reversed_core_no_simd (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++) {
                        Element alpha = pow[i];
                        for (size_t j = 0; j < w; j += 1)
                            Butterfly_DIT (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                             alpha);
                    }
                }
            }

            /* NoSimd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                DIT_reversed_core_no_simd (coeffs, w, f, pow);
            }

            /* Simd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                for ( ; w >= Simd::vect_size; f <<= 1, w >>= 1, pow -= f) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        for (size_t j = 0; j < w; j += Simd::vect_size)
                            Butterfly_DIT (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                                            alpha, P, U);
                    }
                }

                DIT_reversed_core_laststeps (coeffs, w, f, pow, P, U);
            }

            /* DIT ************************************************************/
            void
            DIT (Element *coeffs) const {
                /* w = 1, f = n / 2 */
                DIT_core (coeffs, 1, n >> 1, pow_w.data() + (n-2));
            }

            void
            DIT_core_no_simd (Element *coeffs, size_t w, size_t f,
                              const Element *pow, size_t bound = 0) const {
                bound = bound ? std::min (bound, n) : n;
                for ( ; w < bound; w <<= 1, f >>= 1, pow -= w) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++)
                        for (size_t j = 0; j < w; j += 1)
                            Butterfly_DIT (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                             pow[j]);
                }
            }

            /* NoSimd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_core (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                DIT_core_no_simd (coeffs, w, f, pow);
            }

            /* Simd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIT_core (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                DIT_core_firststeps (coeffs, w, f, pow, P, U);

                for ( ; w < n; w <<= 1, f >>= 1, pow -= w) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++)
                        for (size_t j = 0; j < w; j += Simd::vect_size)
                            Butterfly_DIT (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                            pow+j, P, U);
                }
            }

            /* DIF reversed ***************************************************/
            void
            DIF_reversed (Element *coeffs) const {
                /* w = 1, f = n / 2 */
                DIF_reversed_core (coeffs, 1, n >> 1, pow_w_br.data());
            }

            void
            DIF_reversed_core_no_simd (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow,
                                                    size_t bound = 0) const {
                bound = bound ? std::min (bound, n) : n;
                for ( ; w < bound; pow += f, w <<= 1, f >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++) {
                        Element alpha = pow[i];
                        for (size_t j = 0; j < w; j += 1)
                            Butterfly_DIF (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                             alpha);
                    }
                }
            }

            /* NoSimd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                DIF_reversed_core_no_simd (coeffs, w, f, pow);
            }

            /* Simd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t U = Simd::set1 (1.0/fld->characteristic());

                DIF_reversed_core_firststeps (coeffs, w, f, pow, P, U);

                for ( ; w < n; pow += f, w <<= 1, f >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        for (size_t j = 0; j < w; j += Simd::vect_size)
                            Butterfly_DIF (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                                            alpha, P, U);
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

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIF (simd_vect_t& A, simd_vect_t& B,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
                simd_vect_t T1, T2;

                /* A+B mod p */
                T1 = SimdMod::add_mod (A, B, P);
                /* A-B mod p */
                T2 = SimdMod::sub_mod (A, B, P);
                /* multiply A-B by alpha and store it in Bptr */
                B = SimdMod::mul_mod (T2, alpha, P, U);

                A = T1;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
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
                T1 = SimdMod::add_mod (A, B, P);
                Simd::store (Aptr, T1);
                /* A-B mod p */
                T2 = SimdMod::sub_mod (A, B, P);
                /* multiply A-B by alpha and store it in Bptr */
                T3 = SimdMod::mul_mod (T2, alpha, P, U);
                Simd::store (Bptr, T3);
#else
                simd_vect_t alpha = Simd::load (alpha_ptr);
                Butterfly_DIF (Aptr, Bptr, alpha, P, U);
#endif
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* A+B mod p */
                T1 = SimdMod::add_mod (A, B, P);
                Simd::store (Aptr, T1);
                /* A-B mod p */
                T2 = SimdMod::sub_mod (A, B, P);
                /* multiply A-B by alpha and store it in Bptr */
                T3 = SimdMod::mul_mod (T2, alpha, P, U);
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

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (simd_vect_t& A, simd_vect_t& B,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
                simd_vect_t T1, T2;

                /* B*alpha mod P */
                T1 = SimdMod::mul_mod (B, alpha, P, U);
                /* A+B*alpha */
                T2 = SimdMod::add_mod (A, T1, P);
                /* A-B*alpha */
                B = SimdMod::sub_mod (A, T1, P);

                A = T2;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
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
                T1 = SimdMod::mul_mod (B, alpha, P, U);
                /* A+B*alpha */
                T2 = SimdMod::add_mod (A, T1, P);
                Simd::store (Aptr, T2);
                /* A-B*alpha */
                T3 = SimdMod::sub_mod (A, T1, P);
                Simd::store (Bptr, T3);
#else
                simd_vect_t alpha = Simd::load (alpha_ptr);
                Butterfly_DIT (Aptr, Bptr, alpha, P, U);
#endif
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& P,
                            const simd_vect_t& U) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* B*alpha mod P */
                T1 = SimdMod::mul_mod (B, alpha, P, U);
                /* A+B*alpha */
                T2 = SimdMod::add_mod (A, T1, P);
                Simd::store (Aptr, T2);
                /* A-B*alpha */
                T3 = SimdMod::sub_mod (A, T1, P);
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

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 2>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                if (n < (Simd::vect_size << 1)) {
                    DIF_core_no_simd (coeffs, w, f, pow);
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C], V2 = [B D] */
                        V1 = Simd::set (coeffs[0], coeffs[2]);
                        V2 = Simd::set (coeffs[1], coeffs[3]);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, P);
                        V2 = SimdMod::sub_mod (V1, V2, P);

                        /* Result in T = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        SimdExtra::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 4>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                if (n < (Simd::vect_size << 1)) {
                    DIF_core_no_simd (coeffs, w, f, pow);
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;
                        /* V1 = [A E B F], V2 = [C G D H] */
                        V1 = Simd::set (coeffs[0], coeffs[4], coeffs[1],
                                                                    coeffs[5]);
                        V2 = Simd::set (coeffs[2], coeffs[6], coeffs[3],
                                                                    coeffs[7]);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W, P, U);
                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, P);
                        V2 = SimdMod::sub_mod (V1, V2, P);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        SimdExtra::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 8>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                if (n < (Simd::vect_size << 1)) {
                    DIF_core_no_simd (coeffs, w, f, pow);
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1],
                                               pow[2], pow[2], pow[3], pow[3]);
                    simd_vect_t W2 = Simd::set (pow[4], pow[4], pow[4], pow[4],
                                                pow[5], pow[5], pow[5], pow[5]);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A I B J C K D L], V2 = [E M F N G O H P] */
                        V1 = Simd::set (coeffs[0], coeffs[8], coeffs[1],
                                        coeffs[9], coeffs[2], coeffs[10],
                                        coeffs[3], coeffs[11]);
                        V2 = Simd::set (coeffs[4], coeffs[12], coeffs[5],
                                        coeffs[13], coeffs[6], coeffs[14],
                                        coeffs[7], coeffs[15]);

                        /*** step *********************************************/
                        Butterfly_DIF (V1, V2, W, P, U);
                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W2, P, U);
                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, P);
                        V2 = SimdMod::sub_mod (V1, V2, P);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        SimdExtra::unpacklohi (V1, V2, T, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 16>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                /* P and U are unused with no_simd fallback */
                DIF_core_no_simd (coeffs, w, f, pow);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 2>* = nullptr>
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                /* P and U are unused with no_simd fallback */
                DIT_reversed_core_no_simd (coeffs, w, f, pow);
            }

            /******************************************************************/
            /* Firststeps for DIT and DIF reversed ****************************/
            /******************************************************************/

            /* Firststeps perform the first log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = 1, f = n/2, and the correct pow pointer.
             */

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 2>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                if (n < (Simd::vect_size << 1)) {
                    DIT_core_no_simd (coeffs, w, f, pow);
                    for ( ; w < n; w <<= 1, f >>= 1, pow -= w) ;
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C], V2 = [B D] */
                        V1 = Simd::set (coeffs[0], coeffs[2]);
                        V2 = Simd::set (coeffs[1], coeffs[3]);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, P);
                        V2 = SimdMod::sub_mod (V1, V2, P);

                        /* Result in T = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        SimdExtra::pack (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    w <<= 1;
                    f >>= 1;
                    pow -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 4>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                if (n < (Simd::vect_size << 1)) {
                    DIT_core_no_simd (coeffs, w, f, pow);
                    for ( ; w < n; w <<= 1, f >>= 1, pow -= w) ;
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    f >>= 2;
                    w <<= 2;
                    pow -= 2;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C E G], V2 = [B D F H] */
                        V1 = Simd::set (coeffs[0], coeffs[2], coeffs[4],
                                                                    coeffs[6]);
                        V2 = Simd::set (coeffs[1], coeffs[3], coeffs[5],
                                                                    coeffs[7]);


                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, P);
                        V2 = SimdMod::sub_mod (V1, V2, P);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        SimdExtra::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, P, U);

                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        SimdExtra::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 8>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                if (n < (Simd::vect_size << 1)) {
                    DIT_core_no_simd (coeffs, w, f, pow);
                    for ( ; w < n; w <<= 1, f >>= 1, pow -= w) ;
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
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

                        /* V1 = [A C E G I K M O], V2 = [B D F H J L N P] */
                        V1 = Simd::set (coeffs[0], coeffs[2], coeffs[4],
                                        coeffs[6], coeffs[8], coeffs[10],
                                                   coeffs[12], coeffs[14]);
                        V2 = Simd::set (coeffs[1], coeffs[3], coeffs[5],
                                        coeffs[7], coeffs[9], coeffs[11],
                                                   coeffs[13], coeffs[15]);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, P);
                        V2 = SimdMod::sub_mod (V1, V2, P);

                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        SimdExtra::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, P, U);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        SimdExtra::pack (V1, V2, V1, V2);

                        /*** third step ***************************************/
                        Butterfly_DIT (V1, V2, W2, P, U);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        SimdExtra::pack (V1, V2, V1, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 16>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const simd_vect_t& P,
                                    const simd_vect_t& U) const {
                /* P and U are unused with no_simd fallback */
                DIT_core_no_simd (coeffs, w, f, pow, Simd::vect_size);
                for ( ; w < Simd::vect_size; w <<= 1, f >>= 1, pow -= w) ;
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 1>* = nullptr>
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                          const Element *&pow,
                                          const simd_vect_t& P,
                                          const simd_vect_t& U) const {
                /* P and U are unused with no_simd fallback */
                DIF_reversed_core_no_simd (coeffs, w, f, pow, Simd::vect_size);
                for ( ; w < Simd::vect_size; pow+=f, w <<= 1, f >>= 1) ;
            }

            /******************************************************************/
            /* Utils **********************************************************/
            /******************************************************************/
            bool
            is_primitive_root (Element w) {
                return FFT_internals::is_primitive_root (*this->fld, w, l2n);
            }

            Element
            compute_primitive_root () {
                return FFT_internals::compute_primitive_root (*this->fld, l2n);
            }

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
                        size_t i_br = FFT_internals::bitreverse (i, l);
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
    class FFT<Field, Simd,
                FFT_utils::enable_if_integral_compatible_simd_t<Field, Simd>> {
        /* Types **************************************************************/
        private:
            using Element = typename Field::Element;
            using Residu_t = typename Field::Residu_t;
            using simd_vect_t = typename std::conditional<std::is_same<Simd, NoSimd<typename Field::Element>>::value, Element, typename Simd::vect_t>::type;
            using elt_vect_t = typename Simd::aligned_vector;
            using SimdMod = SimdModular<Field, Simd>;
            using SimdExtra = FFT_utils::SimdExtra<Field, Simd>;

        public:
            /* constructor ****************************************************/
            FFT (const Field& F, size_t k, Element w = 0)
                                            : fld(&F), l2n(k), n(1UL << l2n),
                                              p(F.characteristic()), p2(p << 1),
                                              pow_w(n-1), pow_w_br(n-1),
                                              pow_wp(n-1), pow_wp_br(n-1) {
                if (!k)
                    throw LinBoxError ("FFT: k must be positive");

                if (w == 0)
                    w = compute_primitive_root ();
                else if (!is_primitive_root (w))
                    throw LinBoxError ("FFT: w is not a primitive k-root of unity");
                init_powers (w);
            }

            /* getters */
            const Field &
            field () const {
                return *fld;
            }

            const Element &
            getRoot() const
            {
                if (l2n == 1) /* if n=2^ln equals 2, return -1 */
                    return fld->mOne;
                else
                    return pow_w[1];
            }


            /******************************************************************/
            /* main functions: FFT direct and inverse *************************/
            /******************************************************************/

            /* Perform a FFT in place on the array 'coeffs' of size n.
             * Input:
             *  - must be < p
             *  - is read in natural order
             * Output:
             *  - are < p
             *  - is written in bitreversed order.
             */
            void
            FFT_direct (Element *coeffs) const {
                DIF (coeffs); /* or DIT_reversed */
            }

            /* Perform an inverse FFT in place on the array 'coeffs' of size n.
             * Input:
             *  - must be < p
             *  - is read in bitreversed order.
             * Output:
             *  - are < p
             *  - is written in natural order
             */
            void
            FFT_inverse (Element *coeffs) const {
                DIT (coeffs); /* or DIF_reversed */
            }

        protected:
            /******************************************************************/
            /* Attributes *****************************************************/
            /******************************************************************/
            const Field *fld;
            size_t l2n; /* log2 of size */
        //uint64_t n; /* 2^l2n */
        size_t n;
            const Residu_t p; /* p = field characteristic */
            const Residu_t p2; /* p2 = 2*p */

            elt_vect_t pow_w; /* Table of roots of unity.
                               * If w = primitive n-th root, then the table is:
                               * 1, w, w^2, ..., w^{n/2-1},
                               * 1, w^2, w^4, ..., w^{n/2-2},
                               * 1, w^4, w^8, ..., w^{n/2-4}
                               * ...
                               * 1, w^{n/8}, w^{n/4}, w^{3n/8},
                               * 1, w^{n/4},
                               * 1.
                               * Its size is n-1.
                               */
            elt_vect_t pow_w_br; /* Same as above with each subarray in
                                  * bitreverse order.
                                  */
            elt_vect_t pow_wp; /* pow_wp[i] := precomp_b (pow_w[i]) */
            elt_vect_t pow_wp_br; /* Same as above with each subarray in
                                   * bitreverse order.
                                   */


            /* NoSimd/Simd operations *****************************************/
            /* reduce */
            void
            reduce_no_simd (Element &v, const Residu_t m) const {
               FFT_internals::sub_if_greater (v, m);
            }
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            reduce (Element *v, const Element m) const {
                reduce_no_simd (*v, m);
            }
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            reduce (Element *v, const simd_vect_t m) const {
                SimdMod::reduce (v, m);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            reduce_coeffs_2p (Element *coeffs) const {
                for (uint64_t i = 0; i < n; i++)
                    reduce_no_simd (coeffs[i], p);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            reduce_coeffs_2p (Element *coeffs) const {
                simd_vect_t P = Simd::set1 (p);
                uint64_t i = 0;
                if (n >= Simd::vect_size)
                    for (; i < n; i += Simd::vect_size)
                        reduce (&coeffs[i], P);
                else
                    for (uint64_t i = 0; i < n; i++)
                        reduce_no_simd (coeffs[i], p);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            reduce_coeffs_4p (Element *coeffs) const {
                for (uint64_t i = 0; i < n; i++) {
                    reduce_no_simd (coeffs[i], p2);
                    reduce_no_simd (coeffs[i], p);
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            reduce_coeffs_4p (Element *coeffs) const {
                simd_vect_t P = Simd::set1 (p);
                simd_vect_t P2 = Simd::set1 (p2);
                uint64_t i = 0;
                if (n >= Simd::vect_size)
                    for (; i < n; i += Simd::vect_size) {
                        reduce (&coeffs[i], P2);
                        reduce (&coeffs[i], P);
                    }
                else
                    for (uint64_t i = 0; i < n; i++) {
                        reduce_no_simd (coeffs[i], p2);
                        reduce_no_simd (coeffs[i], p);
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
            void
            DIF (Element *coeffs) const {
                /* w = n/2, f = 1 */
                DIF_core (coeffs, n >> 1, 1, pow_w.data(), pow_wp.data());
                reduce_coeffs_2p (coeffs);
            }

            void
            DIF_core_no_simd (Element *coeffs, size_t w, size_t f,
                              const Element *pow, const Element *powp) const {
                for ( ; w > 0; pow += w, powp += w, f <<= 1, w >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++)
                        for (size_t j = 0; j < w; j += 1)
                            Butterfly_DIF (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                        pow[j], powp[j], p2);
                }
            }

            /* NoSimd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_core (Element *coeffs, size_t w, size_t f, const Element *pow,
                                                    const Element *powp) const {
                DIF_core_no_simd (coeffs, w, f, pow, powp);
            }

            /* Simd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIF_core (Element *coeffs, size_t w, size_t f, const Element *pow,
                                                    const Element *powp) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                for ( ; w >= Simd::vect_size; pow += w, powp += w, f <<= 1,
                                                                   w >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++)
                        for (size_t j = 0; j < w; j += Simd::vect_size)
                            Butterfly_DIF (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                                        pow+j, powp+j, P, P2);
                }

                DIF_core_laststeps (coeffs, w, f, pow, powp, P, P2);
            }

            /* DIT reversed ***************************************************/
            void
            DIT_reversed (Element *coeffs) const {
                /* w = n/2, f = 1 */
                DIT_reversed_core (coeffs, n >> 1, 1, pow_w_br.data()+ (n-2),
                                                      pow_wp_br.data()+ (n-2));
                reduce_coeffs_4p (coeffs);
            }

            void
            DIT_reversed_core_no_simd (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow,
                                                    const Element *powp) const {
                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f, powp -= f) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++) {
                        Element alpha = pow[i];
                        Element alphap = powp[i];
                        for (size_t j = 0; j < w; j += 1)
                            Butterfly_DIT (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                             alpha, alphap, p2);
                    }
                }
            }

            /* NoSimd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                               const Element *pow, const Element *powp) const {
                DIT_reversed_core_no_simd (coeffs, w, f, pow, powp);
            }

            /* Simd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIT_reversed_core (Element *coeffs, size_t w, size_t f,
                               const Element *pow, const Element *powp) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                for ( ; w >= Simd::vect_size; f <<= 1, w >>= 1, pow -= f,
                                                                powp -= f) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        simd_vect_t alphap = Simd::set1 (powp[i]);
                        for (size_t j = 0; j < w; j += Simd::vect_size)
                            Butterfly_DIT (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                                        alpha, alphap, P, P2);
                    }
                }

                DIT_reversed_core_laststeps (coeffs, w, f, pow, powp, P, P2);
            }

            /* DIT ************************************************************/
            void
            DIT (Element *coeffs) const {
                /* w = 1, f = n / 2 */
                DIT_core (coeffs, 1, n >> 1, pow_w.data() + (n-2),
                                             pow_wp.data() + (n-2));
                reduce_coeffs_4p (coeffs);
            }

            void
            DIT_core_no_simd (Element *coeffs, size_t w, size_t f,
                              const Element *pow, const Element *powp,
                              size_t bound = 0) const {
                bound = bound ? std::min (bound, n) : n;
                for ( ; w < bound; w <<= 1, f >>= 1, pow -= w, powp -= w) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++)
                        for (size_t j = 0; j < w; j += 1)
                            Butterfly_DIT (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                        pow[j], powp[j], p2);
                }
            }

            /* NoSimd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_core (Element *coeffs, size_t w, size_t f, const Element *pow,
                                                    const Element *powp) const {
                DIT_core_no_simd (coeffs, w, f, pow, powp);
            }

            /* Simd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIT_core (Element *coeffs, size_t w, size_t f, const Element *pow,
                                                    const Element *powp) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                DIT_core_firststeps (coeffs, w, f, pow, powp, P, P2);

                for ( ; w < n; w <<= 1, f >>= 1, pow -= w, powp -= w) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++)
                        for (size_t j = 0; j < w; j += Simd::vect_size)
                            Butterfly_DIT (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                            pow+j, powp+j, P, P2);
                }
            }

            /* DIF reversed ***************************************************/
            void
            DIF_reversed (Element *coeffs) const {
                /* w = 1, f = n / 2 */
                DIF_reversed_core (coeffs, 1, n >> 1, pow_w_br.data(),
                                                      pow_wp_br.data());
                reduce_coeffs_2p (coeffs);
            }

            void
            DIF_reversed_core_no_simd (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow,
                                                    const Element *powp,
                                                    size_t bound = 0) const {
                bound = bound ? std::min (bound, n) : n;
                for ( ; w < bound; pow += f, powp += f, w <<= 1, f >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++) {
                        Element alpha = pow[i];
                        Element alphap = powp[i];
                        for (size_t j = 0; j < w; j += 1)
                            Butterfly_DIF (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                             alpha, alphap, p2);
                    }
                }
            }

            /* NoSimd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                               const Element *pow, const Element *powp) const {
                DIF_reversed_core_no_simd (coeffs, w, f, pow, powp);
            }

            /* Simd */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIF_reversed_core (Element *coeffs, size_t w, size_t f,
                               const Element *pow, const Element *powp) const {
                simd_vect_t P = Simd::set1 (fld->characteristic());
                simd_vect_t P2 = Simd::set1 (fld->characteristic() << 1);

                DIF_reversed_core_firststeps (coeffs, w, f, pow, powp, P, P2);

                for ( ; w < n; pow += f, powp += f, w <<= 1, f >>= 1) {
                    Element *Aptr = coeffs;
                    Element *Bptr = coeffs + w;
                    for (size_t i = 0; i < f; i++) {
                        simd_vect_t alpha = Simd::set1 (pow[i]);
                        simd_vect_t alphap = Simd::set1 (powp[i]);
                        for (size_t j = 0; j < w; j += Simd::vect_size)
                            Butterfly_DIF (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                                        alpha, alphap, P, P2);
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
                FFT_internals::sub_if_greater (A, p2); /* A -= 2*p if A >= 2p */
                B = tmp + (p2 - B);
                this->fld->mul_precomp_b_without_reduction (B, B, alpha, alphap);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIF (simd_vect_t& A, simd_vect_t& B,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
                simd_vect_t T1, T2, T3;

                /* A+B mod 2p */
                T1 = SimdMod::add_mod (A, B, P2);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha */
                B = SimdMod::mul_mod (T3, alpha, P, alphap);

                A = T1;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
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
                T1 = SimdMod::add_mod (A, B, P2);
                Simd::store (Aptr, T1);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha and store it in Bptr */
                T4 = SimdMod::mul_mod (T3, alpha, P, alphap);
                Simd::store (Bptr, T4);
#else
                simd_vect_t alpha = Simd::load (alpha_ptr);
                simd_vect_t alphap = Simd::load (alphap_ptr);
                Butterfly_DIF (Aptr, Bptr, alpha, alphap, P, P2);
#endif
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3, T4;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* A+B mod 2p and store it in Aptr */
                T1 = SimdMod::add_mod (A, B, P2);
                Simd::store (Aptr, T1);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha and store it in Bptr */
                T4 = SimdMod::mul_mod (T3, alpha, P, alphap);
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
                FFT_internals::sub_if_greater (A, p2); /* A -= 2*p if A >= 2p */
                Element tmp;
                this->fld->mul_precomp_b_without_reduction (tmp, B, alpha, alphap);
                B = A + (p2 - tmp);
                A += tmp;
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (simd_vect_t& A, simd_vect_t& B,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
                simd_vect_t T1, T2, T3;
                T1 = SimdMod::reduce (A, P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdMod::mul_mod (B, alpha, P, alphap);
                /* A+B*alpha */
                A = Simd::add (T1, T2);
                /* A-B*alpha (computed as A+(2p-B*alpha)) */
                T3 = Simd::sub (P2, T2);
                B = Simd::add (T1, T3);
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
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

                T1 = SimdMod::reduce (A, P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdMod::mul_mod (B, alpha, P, alphap);
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

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                            const simd_vect_t& alpha, const simd_vect_t& alphap,
                            const simd_vect_t& P, const simd_vect_t& P2) const {
#ifndef _FFT_NO_INLINED_BUTTERFLY
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                T1 = SimdMod::reduce (A, P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdMod::mul_mod (B, alpha, P, alphap);
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

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 2>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                if (n < (Simd::vect_size << 1)) {
                    DIF_core_no_simd (coeffs, w, f, pow, powp);
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C], V2 = [B D] */
                        V1 = Simd::set (coeffs[0], coeffs[2]);
                        V2 = Simd::set (coeffs[1], coeffs[3]);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, P2);
                        V2 = SimdMod::sub_mod (V1, V2, P2);

                        /* Result in T = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        SimdExtra::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 4>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                if (n < (Simd::vect_size << 1)) {
                    DIF_core_no_simd (coeffs, w, f, pow, powp);
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    simd_vect_t Wp = Simd::set(powp[0],powp[0],powp[1],powp[1]);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;
                        /* V1 = [A E B F], V2 = [C G D H] */
                        V1 = Simd::set (coeffs[0], coeffs[4], coeffs[1],
                                                                    coeffs[5]);
                        V2 = Simd::set (coeffs[2], coeffs[6], coeffs[3],
                                                                    coeffs[7]);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);
                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, P2);
                        V2 = SimdMod::sub_mod (V1, V2, P2);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        SimdExtra::unpacklohi (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 8>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                if (n < (Simd::vect_size << 1)) {
                    DIF_core_no_simd (coeffs, w, f, pow, powp);
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
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

                        /* V1 = [A I B J C K D L], V2 = [E M F N G O H P] */
                        V1 = Simd::set (coeffs[0], coeffs[8], coeffs[1],
                                        coeffs[9], coeffs[2], coeffs[10],
                                        coeffs[3], coeffs[11]);
                        V2 = Simd::set (coeffs[4], coeffs[12], coeffs[5],
                                        coeffs[13], coeffs[6], coeffs[14],
                                        coeffs[7], coeffs[15]);

                        /*** step *********************************************/
                        Butterfly_DIF (V1, V2, W, Wp, P, P2);
                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W2, W2p, P, P2);
                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, P2);
                        V2 = SimdMod::sub_mod (V1, V2, P2);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        SimdExtra::unpacklohi (V1, V2, T, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 16>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                /* P and P2 are unused with no_simd fallback */
                DIF_core_no_simd (coeffs, w, f, pow, powp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 2>* = nullptr>
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                                const Element *&pow,
                                                const Element *&powp,
                                                const simd_vect_t& P,
                                                const simd_vect_t& P2) const {
                /* P and P2 are unused with no_simd fallback */
                DIT_reversed_core_no_simd (coeffs, w, f, pow, powp);
            }

            /******************************************************************/
            /* Firststeps for DIT and DIF reversed ****************************/
            /******************************************************************/

            /* Firststeps perform the first log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = 1, f = n/2, and the correct pow pointer.
             */

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 2>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                if (n < (Simd::vect_size << 1)) {
                    DIT_core_no_simd (coeffs, w, f, pow, powp);
                    for ( ; w < n; w <<= 1, f >>= 1, pow -= w, powp -= w) ;
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C], V2 = [B D] */
                        V1 = Simd::set (coeffs[0], coeffs[2]);
                        V2 = Simd::set (coeffs[1], coeffs[3]);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, P2);
                        V2 = SimdMod::sub_mod (V1, V2, P2);

                        /* Result in T = [A C]  and V2 = [B D]
                         * Transform to V1 = [A B], V2 = [C D] and store
                         */
                        SimdExtra::pack (V1, V2, T, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    w <<= 1;
                    f >>= 1;
                    pow -= w;
                    powp -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 4>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                if (n < (Simd::vect_size << 1)) {
                    DIT_core_no_simd (coeffs, w, f, pow, powp);
                    for ( ; w < n; w <<= 1, f >>= 1, pow -= w, powp -= w) ;
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
                    f >>= 2;
                    w <<= 2;
                    pow -= 2;
                    powp -= 2;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    simd_vect_t Wp = Simd::set (powp[0],powp[0],powp[1],powp[1]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C E G], V2 = [B D F H] */
                        V1 = Simd::set (coeffs[0], coeffs[2], coeffs[4],
                                                                    coeffs[6]);
                        V2 = Simd::set (coeffs[1], coeffs[3], coeffs[5],
                                                                    coeffs[7]);


                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, P2);
                        V2 = SimdMod::sub_mod (V1, V2, P2);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        SimdExtra::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);

                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        SimdExtra::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                    powp -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 8>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                if (n < (Simd::vect_size << 1)) {
                    DIT_core_no_simd (coeffs, w, f, pow, powp);
                    for ( ; w < n; w <<= 1, f >>= 1, pow -= w, powp -= w) ;
                } else {
                    const constexpr size_t incr = Simd::vect_size << 1;
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

                        /* V1 = [A C E G I K M O], V2 = [B D F H J L N P] */
                        V1 = Simd::set (coeffs[0], coeffs[2], coeffs[4],
                                        coeffs[6], coeffs[8], coeffs[10],
                                                   coeffs[12], coeffs[14]);
                        V2 = Simd::set (coeffs[1], coeffs[3], coeffs[5],
                                        coeffs[7], coeffs[9], coeffs[11],
                                                   coeffs[13], coeffs[15]);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, P2);
                        V2 = SimdMod::sub_mod (V1, V2, P2);

                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        SimdExtra::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, Wp, P, P2);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        SimdExtra::pack (V1, V2, V1, V2);

                        /*** third step ***************************************/
                        Butterfly_DIT (V1, V2, W2, W2p, P, P2);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        SimdExtra::pack (V1, V2, V1, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                    powp -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 16>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                /* P and P2 are unused with no_simd fallback */
                DIT_core_no_simd (coeffs, w, f, pow, powp, Simd::vect_size);
                for ( ; w < Simd::vect_size; w <<= 1, f >>= 1, pow-=w, powp-=w);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 2>* = nullptr>
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                    const Element *&pow, const Element *&powp,
                                    const simd_vect_t& P,
                                    const simd_vect_t& P2) const {
                /* P and P2 are unused with no_simd fallback */
                DIF_reversed_core_no_simd (coeffs, w, f, pow, powp,
                                                            Simd::vect_size);
                for ( ; w < Simd::vect_size; pow+=f, powp+=f, w <<= 1, f >>= 1);
            }

            /******************************************************************/
            /* Utils **********************************************************/
            /******************************************************************/
            bool
            is_primitive_root (Element w) {
                return FFT_internals::is_primitive_root (*this->fld, w, l2n);
            }

            Element
            compute_primitive_root () {
                return FFT_internals::compute_primitive_root (*this->fld, l2n);
            }

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
                        size_t i_br = FFT_internals::bitreverse (i, l);
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
#endif // __LINBOX_polynomial_fft_init_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
