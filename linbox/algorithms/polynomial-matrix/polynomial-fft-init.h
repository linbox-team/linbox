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

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-utils.h"

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
    class FFT__;

    template <typename Field, typename Simd>
    class FFT__<Field, Simd,
                FFT_utils::enable_if_floating_compatible_simd_t<Field, Simd>> {
        /* Types **************************************************************/
        private:
            using Element = typename Field::Element;
            using Residu_t = typename Field::Residu_t;
            using simd_vect_t = typename std::conditional<std::is_same<Simd, NoSimd<typename Field::Element>>::value, Element, typename Simd::vect_t>::type;
            using elt_vect_t = typename Simd::aligned_vector;
            using SimdMod = SimdModular<Field, Simd>;
            using SimdExtra = FFT::SimdExtra<Field, Simd>;

        public:
            /* constructor ****************************************************/
            FFT__ (const Field& F, size_t k, Element w = 0)
                                            : fld(&F), l2n(k), n(1UL << l2n),
                                              p(F.characteristic()), P(set1(p)),
                                              U(set1(1.0/p)), pow_w(n-1) ,
                                              pow_w_br(n-1) {
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


            /* main functions: FFT direct and inverse *************************/
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
            FFT_inverse (Element *coeffs, bool final_division = true) const {
                DIT (coeffs); /* or DIF_reversed */
                if (final_division) {
                    /* divide by n */
                    Element inv_n;
                    this->fld->init (inv_n, n);
                    this->fld->invin (inv_n);
                    scal_coeffs (coeffs, inv_n);
                }
            }

        protected:
            /* Attributes *****************************************************/
            const Field *fld;
            size_t l2n; /* log2 of size */
            uint64_t n; /* 2^l2n */
            const Residu_t p; /* p = field characteristic */
            simd_vect_t P; /* same as above, for vectorized computation */
            simd_vect_t U; /* U = 1/P, for vectorized computation */

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

            /* NoSimd/Simd operations *****************************************/
            /* set1 */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            simd_vect_t
            set1 (const Element v) const {
                return Simd::set1 (v);
            }
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            Element
            set1 (const Element v) const {
                return v;
            }

            /* scal_coeffs */
            void
            scal_coeffs_no_simd (Element *coeffs, const Element& m) const {
                for (size_t i = 0; i < n; i++)
                    this->fld->mulin (coeffs[i], m);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            scal_coeffs (Element *coeffs, const Element& m) const {
                scal_coeffs_no_simd (coeffs, m);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            scal_coeffs (Element *coeffs, const Element& m) const {
                if (n < Simd::vect_size) {
                    scal_coeffs_no_simd (coeffs, m);
                } else {
                    simd_vect_t M = Simd::set1 (m);
                    Element *cptr = coeffs;
                    for (size_t i = 0; i < n; i += Simd::vect_size,
                                                    cptr += Simd::vect_size) {
                        simd_vect_t C = Simd::load (cptr);
                        C = SimdMod::mul_mod (C, M, this->P, this->U);
                        Simd::store (cptr, C);
                    }
                }
            }

            /* DIF ************************************************************/

            /* Compute A[i]+B[i], (A[i]-B[i])*alpha[i],
             * for 0 <= i < simd::vect_size.
             * Input must satisfy:
             *  - 0 <= A[i],B[i],alpha[i] < p
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < p
             */
            /* no_simd version is defined for every type of Simd */
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
                                            const simd_vect_t& alpha) const {
                simd_vect_t T1, T2;

                /* A+B mod p */
                T1 = SimdMod::add_mod (A, B, this->P);
                /* A-B mod p */
                T2 = SimdMod::sub_mod (A, B, this->P);
                /* multiply A-B by alpha and store it in Bptr */
                B = SimdMod::mul_mod (T2, alpha, this->P, this->U);

                A = T1;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                                            const simd_vect_t& alpha) const {
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* A+B mod p */
                T1 = SimdMod::add_mod (A, B, this->P);
                Simd::store (Aptr, T1);
                /* A-B mod p */
                T2 = SimdMod::sub_mod (A, B, this->P);
                /* multiply A-B by alpha and store it in Bptr */
                T3 = SimdMod::mul_mod (T2, alpha, this->P, this->U);
                Simd::store (Bptr, T3);
            }
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                                            const Element* alpha_ptr) const {
                simd_vect_t A, B, T1, T2, T3, alpha;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);
                alpha = Simd::load (alpha_ptr);

                /* A+B mod p */
                T1 = SimdMod::add_mod (A, B, this->P);
                Simd::store (Aptr, T1);
                /* A-B mod p */
                T2 = SimdMod::sub_mod (A, B, this->P);
                /* multiply A-B by alpha and store it in Bptr */
                T3 = SimdMod::mul_mod (T2, alpha, this->P, this->U);
                Simd::store (Bptr, T3);
            }
            /* See DIF_core for the meaning of the parameter */
            /* no_simd version is defined for every type of Simd */
            void
            DIF_core_one_step_no_simd (Element *coeffs, const size_t w,
                                       const size_t f, const Element *pow)
                                                                        const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    for (size_t j = 0; j < w; j += 1) {
                        Butterfly_DIF (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                                        pow[j]);
                    }
                }
            }

            void
            DIF_reversed_core_one_step_no_simd (Element *coeffs, const size_t w,
                                             const size_t f, const Element *pow)
                                                                        const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    const Element alpha = pow[i];
                    for (size_t j = 0; j < w; j += 1) {
                        Butterfly_DIF (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                                        alpha);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_core_one_step (Element *coeffs, const size_t w, const size_t f,
                                                    const Element *pow) const {
                DIF_core_one_step_no_simd (coeffs, w, f, pow);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIF_core_one_step (Element *coeffs, const size_t w, const size_t f,
                                                    const Element *pow) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    for (size_t j = 0; j < w; j += Simd::vect_size) {
                        Butterfly_DIF (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j, pow+j);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_reversed_core_one_step (Element *coeffs, const size_t w,
                                        const size_t f, const Element *pow)
                                                                        const {
                DIF_reversed_core_one_step_no_simd (coeffs, w, f, pow);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIF_reversed_core_one_step (Element *coeffs, const size_t w,
                                        const size_t f, const Element *pow)
                                                                        const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    simd_vect_t alpha = Simd::set1 (pow[i]);
                    for (size_t j = 0; j < w; j += Simd::vect_size) {
                        Butterfly_DIF (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j, alpha);
                    }
                }
            }

            /* Laststeps perform the last log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = Simd::vect_size/2, f = 4*n/Simd::vect_size,
             * pow = pow_w.data() - Simd::vect_size.
             */
            void
            DIF_core_laststeps_fallback (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow) const {
                for ( ; w > 0; pow+=w, f <<= 1, w >>= 1)
                    DIF_core_one_step_no_simd (coeffs, w, f, pow);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                                 const Element *pow) const {
                /* nothing to do for laststeps of S::vect_size == 1 */
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 2>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                                 const Element *pow) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIF_core_laststeps_fallback (coeffs, w, f, pow);
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C], V2 = [B D] */
                        V1 = Simd::set (coeffs[0], coeffs[2]);
                        V2 = Simd::set (coeffs[1], coeffs[3]);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, this->P);
                        V2 = SimdMod::sub_mod (V1, V2, this->P);

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
                                                 const Element *pow) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIF_core_laststeps_fallback (coeffs, w, f, pow);
                } else {
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;
                        /* V1 = [A E B F], V2 = [C G D H] */
                        V1 = Simd::set (coeffs[0], coeffs[4], coeffs[1],
                                                                    coeffs[5]);
                        V2 = Simd::set (coeffs[2], coeffs[6], coeffs[3],
                                                                    coeffs[7]);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W);
                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, this->P);
                        V2 = SimdMod::sub_mod (V1, V2, this->P);

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
                                                 const Element *pow) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIF_core_laststeps_fallback (coeffs, w, f, pow);
                } else {
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
                        Butterfly_DIF (V1, V2, W);
                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W2);
                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, this->P);
                        V2 = SimdMod::sub_mod (V1, V2, this->P);

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
                                                 const Element *pow) const {
                DIF_core_laststeps_fallback (coeffs, w, f, pow);
            }

            /* Firststeps perform the first log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = 1, f = n/2, pow = pow_w.data().
             * Note that for now, there is no special code for first steps with
             * SIMD, we just call the no_simd code.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 1>* = nullptr>
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                 const Element *&pow) const {
                size_t bound = (n < Simd::vect_size) ? n : Simd::vect_size;
                for ( ; w < bound; pow+=f, w <<= 1, f >>= 1)
                    DIF_reversed_core_one_step_no_simd (coeffs, w, f, pow);
            }

            /* core functions for DIF */
            void
            DIF_core (Element *coeffs) const {
                /* n : length of the array 'coeffs' (always of power of 2)
                 * f : number of families of butterflies
                 * w : width of butterflies
                 * (outmost) loop invariant : 2*f*w == n
                 */
                size_t w = n >> 1, f = 1;
                const Element *pow_ptr = pow_w.data();
                for ( ; w >= Simd::vect_size; pow_ptr += w, f <<= 1, w >>= 1)
                    DIF_core_one_step (coeffs, w, f, pow_ptr);

                DIF_core_laststeps (coeffs, w, f, pow_ptr);
            }

            void
            DIF_reversed_core (Element *coeffs) const {
                /* n : length of the array 'coeffs' (always of power of 2)
                 * f : number of families of butterflies
                 * w : width of butterflies
                 * (outmost) loop invariant : 2*f*w == n
                 */
                size_t w = 1, f = n >> 1;
                const Element *pow_ptr = pow_w_br.data();

                DIF_reversed_core_firststeps (coeffs, w, f, pow_ptr);

                for ( ; w < n; pow_ptr += f, w <<= 1, f >>= 1)
                    DIF_reversed_core_one_step (coeffs, w, f, pow_ptr);
            }

            void
            DIF (Element *coeffs) const {
                DIF_core (coeffs);
            }

            void
            DIF_reversed (Element *coeffs) const {
                DIF_reversed_core (coeffs);
            }

            /* DIT ************************************************************/
            /* Compute A[i]+B[i]*alpha[i], A[i]-B[i]*alpha[i],
             * for 0 <= i < simd::vect_size.
             * Input must satisfy:
             *  - 0 <= A[i],B[i],alpha[i] < p
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < p
             */
            /* no_simd version is defined for every type of Simd */
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
                                            const simd_vect_t& alpha) const {
                simd_vect_t T1, T2;

                /* B*alpha mod P */
                T1 = SimdMod::mul_mod (B, alpha, this->P, this->U);
                /* A+B*alpha */
                T2 = SimdMod::add_mod (A, T1, this->P);
                /* A-B*alpha */
                B = SimdMod::sub_mod (A, T1, this->P);

                A = T2;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                                            const simd_vect_t& alpha) const {
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* B*alpha mod P */
                T1 = SimdMod::mul_mod (B, alpha, this->P, this->U);
                /* A+B*alpha */
                T2 = SimdMod::add_mod (A, T1, this->P);
                Simd::store (Aptr, T2);
                /* A-B*alpha */
                T3 = SimdMod::sub_mod (A, T1, this->P);
                Simd::store (Bptr, T3);
            }
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                                            const Element *alpha_ptr) const {
                simd_vect_t A, B, T1, T2, T3, alpha;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);
                alpha = Simd::load (alpha_ptr);

                /* B*alpha mod P */
                T1 = SimdMod::mul_mod (B, alpha, this->P, this->U);
                /* A+B*alpha */
                T2 = SimdMod::add_mod (A, T1, this->P);
                Simd::store (Aptr, T2);
                /* A-B*alpha */
                T3 = SimdMod::sub_mod (A, T1, this->P);
                Simd::store (Bptr, T3);
            }
            /* See DIT_core for the meaning of the parameter */
            /* No simd version is defined for every type of Simd */
            void
            DIT_core_one_step_no_simd (Element *coeffs, const size_t w,
                                       const size_t f, const Element *pow)
                                                                        const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    for (size_t j = 0; j < w; j += 1) {
                        Butterfly_DIT (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                                        pow[j]);
                    }
                }
            }

            void
            DIT_reversed_core_one_step_no_simd (Element *coeffs, const size_t w,
                                             const size_t f, const Element *pow)
                                                                        const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    const Element alpha = pow[i];
                    for (size_t j = 0; j < w; j += 1) {
                        Butterfly_DIT (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                                        alpha);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_core_one_step (Element *coeffs, const size_t w, const size_t f,
                                                    const Element *pow) const {
                DIT_core_one_step_no_simd (coeffs, w, f, pow);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIT_core_one_step (Element *coeffs, const size_t w, const size_t f,
                                                    const Element *pow) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    for (size_t j = 0; j < w; j += Simd::vect_size) {
                        Butterfly_DIT (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j, pow+j);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_reversed_core_one_step (Element *coeffs, const size_t w,
                                        const size_t f, const Element *pow)
                                                                        const {
                DIT_reversed_core_one_step_no_simd (coeffs, w, f, pow);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIT_reversed_core_one_step (Element *coeffs, const size_t w,
                                        const size_t f, const Element *pow)
                                                                        const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    simd_vect_t alpha = Simd::set1 (pow[i]);
                    for (size_t j = 0; j < w; j += Simd::vect_size) {
                        Butterfly_DIT (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j, alpha);
                    }
                }
            }

            /* Firststeps perform the first log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = 1, f = n/2, pow = pow_w.data() + (n-2).
             */
            void
            DIT_core_firststeps_fallback (Element *coeffs, size_t &w, size_t &f,
                                                 const Element *&pow) const {
                size_t bound = (n < Simd::vect_size) ? n : Simd::vect_size;
                for ( ; w < bound; w <<= 1, f >>= 1, pow -= w)
                    DIT_core_one_step_no_simd (coeffs, w, f, pow);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                 const Element *&pow) const {
                /* nothing to do for firststeps of S::vect_size == 1 */
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 2>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                 const Element *&pow) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIT_core_firststeps_fallback (coeffs, w, f, pow);
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C], V2 = [B D] */
                        V1 = Simd::set (coeffs[0], coeffs[2]);
                        V2 = Simd::set (coeffs[1], coeffs[3]);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, this->P);
                        V2 = SimdMod::sub_mod (V1, V2, this->P);

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
                                                 const Element *&pow) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIT_core_firststeps_fallback (coeffs, w, f, pow);
                } else {
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
                        T = SimdMod::add_mod (V1, V2, this->P);
                        V2 = SimdMod::sub_mod (V1, V2, this->P);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        SimdExtra::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W);

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
                                                    const Element *&pow) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIT_core_firststeps_fallback (coeffs, w, f, pow);
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

                        /* V1 = [A C E G I K M O], V2 = [B D F H J L N P] */
                        V1 = Simd::set (coeffs[0], coeffs[2], coeffs[4],
                                        coeffs[6], coeffs[8], coeffs[10],
                                                   coeffs[12], coeffs[14]);
                        V2 = Simd::set (coeffs[1], coeffs[3], coeffs[5],
                                        coeffs[7], coeffs[9], coeffs[11],
                                                   coeffs[13], coeffs[15]);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, this->P);
                        V2 = SimdMod::sub_mod (V1, V2, this->P);

                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        SimdExtra::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        SimdExtra::pack (V1, V2, V1, V2);

                        /*** third step ***************************************/
                        Butterfly_DIT (V1, V2, W2);

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
                                                 const Element *&pow) const {
                DIT_core_firststeps_fallback (coeffs, w, f, pow);
            }

            /* Laststeps perform the last log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = Simd::vect_size/2, f = 4*n/Simd::vect_size,
             * pow = pow_w.data() - Simd::vect_size.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 1>* = nullptr>
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                                 const Element *pow) const {
                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f)
                    DIT_reversed_core_one_step_no_simd (coeffs, w, f, pow);
            }

            /* core functions for DIT */
            void
            DIT_core (Element *coeffs) const {
                /* n : length of the array 'coeffs' (always of power of 2)
                 * f : number of families of butterflies
                 * w : width of butterflies
                 * (outmost) loop invariant : 2*f*w == n
                 */
                size_t w = 1, f = n >> 1;
                const Element *pow_ptr = pow_w.data() + (n-2);

                DIT_core_firststeps (coeffs, w, f, pow_ptr);

                for ( ; w < n; w <<= 1, f >>= 1, pow_ptr -= w)
                    DIT_core_one_step (coeffs, w, f, pow_ptr);
            }

            void
            DIT_reversed_core (Element *coeffs) const {
                /* n : length of the array 'coeffs' (always of power of 2)
                 * f : number of families of butterflies
                 * w : width of butterflies
                 * (outmost) loop invariant : 2*f*w == n
                 */
                size_t w = n >> 1, f = 1;
                const Element *pow_ptr = pow_w_br.data() + (n-2);
                for ( ; w >= Simd::vect_size; f <<= 1, w >>= 1, pow_ptr -= f)
                    DIT_reversed_core_one_step (coeffs, w, f, pow_ptr);

                DIT_reversed_core_laststeps (coeffs, w, f, pow_ptr);
            }

            void
            DIT (Element *coeffs) const {
                DIT_core (coeffs);
            }

            void
            DIT_reversed (Element *coeffs) const {
                DIT_reversed_core (coeffs);
            }

            /* Utils **********************************************************/
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
    class FFT__<Field, Simd,
                FFT_utils::enable_if_integral_compatible_simd_t<Field, Simd>> {
        /* Types **************************************************************/
        private:
            using Element = typename Field::Element;
            using Residu_t = typename Field::Residu_t;
            using simd_vect_t = typename std::conditional<std::is_same<Simd, NoSimd<typename Field::Element>>::value, Element, typename Simd::vect_t>::type;
            using elt_vect_t = typename Simd::aligned_vector;
            using SimdMod = SimdModular<Field, Simd>;
            using SimdExtra = FFT::SimdExtra<Field, Simd>;

        public:
            /* constructor ****************************************************/
            FFT__ (const Field& F, size_t k, Element w = 0)
                                            : fld(&F), l2n(k), n(1UL << l2n),
                                              p(F.characteristic()), P(set1(p)),
                                              p2(p << 1), P2(set1(p2)),
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


            /* main functions: FFT direct and inverse *************************/
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
            FFT_inverse (Element *coeffs, bool final_division = true) const {
                DIT (coeffs); /* or DIF_reversed */
                if (final_division) {
                    /* divide by n */
                    Element inv_n;
                    typename Field::Compute_t precomp;
                    this->fld->init (inv_n, n);
                    this->fld->invin (inv_n);
                    this->fld->precomp_b (precomp, inv_n);
                    scal_coeffs (coeffs, inv_n, static_cast<Element>(precomp));
                }
            }

        protected:
            /* Attributes *****************************************************/
            const Field *fld;
            size_t l2n; /* log2 of size */
            uint64_t n; /* 2^l2n */
            const Residu_t p; /* p = field characteristic */
            simd_vect_t P; /* same as above, for vectorized computation */
            const Residu_t p2; /* p2 = 2*p */
            simd_vect_t P2; /* same as above, for vectorized computation */

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
            /* set1 */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            simd_vect_t
            set1 (const Element v) const {
                return Simd::set1 (v);
            }
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            Element
            set1 (const Element v) const {
                return v;
            }

            /* scal_coeffs */
            /* TODO use mul_precomp */
            void
            scal_coeffs_no_simd (Element *coeffs, const Element& m,
                                                  const Element& mp) const {
                for (size_t i = 0; i < n; i++)
                    this->fld->mul_precomp_b_without_reduction (coeffs[i],
                                                                coeffs[i], m,
                                                                            mp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            scal_coeffs (Element *coeffs, const Element& m,
                                                  const Element& mp) const {
                scal_coeffs_no_simd (coeffs, m, mp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            scal_coeffs (Element *coeffs, const Element& m,
                                                  const Element& mp) const {
                if (n < Simd::vect_size) {
                    scal_coeffs_no_simd (coeffs, m, mp);
                } else {
                    simd_vect_t M = Simd::set1 (m);
                    simd_vect_t Mp = Simd::set1 (mp);
                    Element *cptr = coeffs;
                    for (size_t i = 0; i < n; i += Simd::vect_size,
                                                    cptr += Simd::vect_size) {
                        simd_vect_t C = Simd::load (cptr);
                        C = SimdMod::mul_mod (C, M, this->P, Mp);
                        Simd::store (cptr, C);
                    }
                }
            }

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

            void
            reduce_coeffs_2p (Element *coeffs) const {
                uint64_t i = 0;
                if (n >= Simd::vect_size)
                    for (; i < n; i += Simd::vect_size)
                        reduce (&coeffs[i], P);
                else
                    for (uint64_t i = 0; i < n; i++)
                        reduce_no_simd (coeffs[i], p);
            }

            void
            reduce_coeffs_4p (Element *coeffs) const {
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

            /* DIF ************************************************************/

            /* Compute A[i]+B[i], (A[i]-B[i])*alpha[i] using Harvey's algorithm.
             * Input must satisfy:
             *  - 0 <= A[i],B[i] < 4*p
             *  - 0 <= alpha[i] < p
             *  - p < 2^#nbits(Element) / 4
             *  - alphap[i] = Floor(alpha[i] * 2^#nbits(Element) / p)
             * Ensure that output satisfy:
             *  - 0 <= A[i],B[i] < 2*p
             *
             * Note: maybe 2^#nbits(Element) should be maxCardinality ? (in p<)
             */
            /* no_simd version is defined for every type of Simd */
            void
            Butterfly_DIF (Element& A, Element& B, const Element& alpha,
                                                const Element& alphap) const {
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
                                            const simd_vect_t& alpha,
                                            const simd_vect_t& alphap) const {
                simd_vect_t T1, T2, T3;

                /* A+B mod 2p */
                T1 = SimdMod::add_mod (A, B, this->P2);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (this->P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha and store it in Bptr */
                B = SimdMod::mul_mod (T3, alpha, this->P, alphap);

                A = T1;
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                                            const simd_vect_t& alpha,
                                            const simd_vect_t& alphap) const {
                simd_vect_t A, B, T1, T2, T3, T4;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                /* A+B mod 2p and store it in Aptr */
                T1 = SimdMod::add_mod (A, B, this->P2);
                Simd::store (Aptr, T1);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (this->P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha and store it in Bptr */
                T4 = SimdMod::mul_mod (T3, alpha, this->P, alphap);
                Simd::store (Bptr, T4);
            }
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIF (Element *Aptr, Element *Bptr,
                                            const Element* alpha_ptr,
                                            const Element* alphap_ptr) const {
                simd_vect_t A, B, T1, T2, T3, T4, alpha, alphap;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);
                alpha = Simd::load (alpha_ptr);
                alphap = Simd::load (alphap_ptr);

                /* A+B mod 2p and store it in Aptr */
                T1 = SimdMod::add_mod (A, B, this->P2);
                Simd::store (Aptr, T1);
                /* A-B mod 2p (computed as A+(2p-B)) */
                T2 = Simd::sub (this->P2, B);
                T3 = Simd::add (A, T2);
                /* multiply A-B by alpha and store it in Bptr */
                T4 = SimdMod::mul_mod (T3, alpha, this->P, alphap);
                Simd::store (Bptr, T4);
            }

            /* See DIF_core for the meaning of the parameter */
            /* no_simd version is defined for every type of Simd */
            void
            DIF_core_one_step_no_simd (Element *coeffs, const size_t w,
                                       const size_t f, const Element *pow,
                                                    const Element *wp) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    for (size_t j = 0; j < w; j += 1) {
                        Butterfly_DIF (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                                pow[j], wp[j]);
                    }
                }
            }

            void
            DIF_reversed_core_one_step_no_simd (Element *coeffs, const size_t w,
                                             const size_t f, const Element *pow,
                                                    const Element *wp) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    const Element alpha = pow[i];
                    const Element alphap = wp[i];
                    for (size_t j = 0; j < w; j += 1) {
                        Butterfly_DIF (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                                alpha, alphap);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_core_one_step (Element *coeffs, const size_t w, const size_t f,
                                                    const Element *pow,
                                                    const Element *wp) const {
                DIF_core_one_step_no_simd (coeffs, w, f, pow, wp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIF_core_one_step (Element *coeffs, const size_t w, const size_t f,
                                                    const Element *pow,
                                                    const Element *wp) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    for (size_t j = 0; j < w; j += Simd::vect_size) {
                        Butterfly_DIF (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j, pow+j,
                                                                         wp+j);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_reversed_core_one_step (Element *coeffs, const size_t w,
                                        const size_t f, const Element *pow,
                                                    const Element *wp) const {
                DIF_reversed_core_one_step_no_simd (coeffs, w, f, pow, wp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIF_reversed_core_one_step (Element *coeffs, const size_t w,
                                        const size_t f, const Element *pow,
                                                    const Element *wp) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    simd_vect_t alpha = Simd::set1 (pow[i]);
                    simd_vect_t alphap = Simd::set1 (wp[i]);
                    for (size_t j = 0; j < w; j += Simd::vect_size) {
                        Butterfly_DIF (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j, alpha,
                                                                        alphap);
                    }
                }
            }

            /* Laststeps perform the last log2(Simd::vect_size) step(s) (if n is
             * large enough). w, f, and pow are passed to avoid recomputation,
             * it is the responsability of the caller to pass the correct value,
             * i.e., if n is large enough, w = Simd::vect_size/2,
             * f = 4*n/Simd::vect_size, pow = pow_w.data() - Simd::vect_size.
             */
            void
            DIF_core_laststeps_fallback (Element *coeffs, size_t w, size_t f,
                                                            const Element *pow,
                                                    const Element *wp) const {
                for ( ; w > 0; pow+=w, wp+=w, f <<= 1, w >>= 1)
                    DIF_core_one_step_no_simd (coeffs, w, f, pow, wp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow,
                                                    const Element *wp) const {
                /* nothing to do for laststeps of S::vect_size == 1 */
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 2>* = nullptr>
            void
            DIF_core_laststeps (Element *coeffs, size_t w, size_t f,
                                                    const Element *pow,
                                                    const Element *wp) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIF_core_laststeps_fallback (coeffs, w, f, pow, wp);
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C], V2 = [B D] */
                        V1 = Simd::set (coeffs[0], coeffs[2]);
                        V2 = Simd::set (coeffs[1], coeffs[3]);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, this->P2);
                        V2 = SimdMod::sub_mod (V1, V2, this->P2);

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
                                                    const Element *pow,
                                                    const Element *wp) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIF_core_laststeps_fallback (coeffs, w, f, pow, wp);
                } else {
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    simd_vect_t Wp = Simd::set (wp[0], wp[0], wp[1], wp[1]);
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;
                        /* V1 = [A E B F], V2 = [C G D H] */
                        V1 = Simd::set (coeffs[0], coeffs[4], coeffs[1],
                                                                    coeffs[5]);
                        V2 = Simd::set (coeffs[2], coeffs[6], coeffs[3],
                                                                    coeffs[7]);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W, Wp);
                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A C E G], V2 = [B D F H]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, this->P2);
                        V2 = SimdMod::sub_mod (V1, V2, this->P2);

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
                                                    const Element *pow,
                                                    const Element *wp) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIF_core_laststeps_fallback (coeffs, w, f, pow, wp);
                } else {
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1],
                                               pow[2], pow[2], pow[3], pow[3]);
                    simd_vect_t Wp = Simd::set (wp[0], wp[0], wp[1], wp[1],
                                                wp[2], wp[2], wp[3], wp[3]);
                    simd_vect_t W2 = Simd::set (pow[4], pow[4], pow[4], pow[4],
                                                pow[5], pow[5], pow[5], pow[5]);
                    simd_vect_t W2p = Simd::set (wp[4], wp[4], wp[4], wp[4],
                                                 wp[5], wp[5], wp[5], wp[5]);
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
                        Butterfly_DIF (V1, V2, W, Wp);
                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last but one step ********************************/
                        Butterfly_DIF (V1, V2, W2, W2p);
                        /* transform into
                         *      V1 = [A C E G I K M O], V2 = [B D F H J L N P]
                         */
                        SimdExtra::unpacklohi (V1, V2, V1, V2);

                        /*** last step (special butterfly with mul by 1) ******/
                        T = SimdMod::add_mod (V1, V2, this->P2);
                        V2 = SimdMod::sub_mod (V1, V2, this->P2);

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
                                                    const Element *pow,
                                                    const Element *wp) const {
                DIF_core_laststeps_fallback (coeffs, w, f, pow, wp);
            }


            /* Firststeps perform the first log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = 1, f = n/2, pow = pow_w.data() + (n-2).
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 1>* = nullptr>
            void
            DIF_reversed_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                    const Element *&pow,
                                                    const Element *&wp) const {
                size_t bound = (n < Simd::vect_size) ? n : Simd::vect_size;
                for ( ; w < bound; pow+=f, wp+=f, w <<= 1, f >>= 1)
                    DIF_reversed_core_one_step_no_simd (coeffs, w, f, pow, wp);
            }

            /* core functions for DIF */
            void
            DIF_core (Element *coeffs) const {
                /* n : length of the array 'coeffs' (always of power of 2)
                 * f : number of families of butterflies
                 * w : width of butterflies
                 * (outmost) loop invariant : 2*f*w == n
                 */
                size_t w = n >> 1, f = 1;
                const Element *pow_ptr = pow_w.data();
                const Element *wp_ptr = pow_wp.data();
                for ( ; w >= Simd::vect_size; pow_ptr += w, wp_ptr += w,
                                                            f <<= 1, w >>= 1)
                    DIF_core_one_step (coeffs, w, f, pow_ptr, wp_ptr);

                DIF_core_laststeps (coeffs, w, f, pow_ptr, wp_ptr);
            }

            void
            DIF_reversed_core (Element *coeffs) const {
                /* n : length of the array 'coeffs' (always of power of 2)
                 * f : number of families of butterflies
                 * w : width of butterflies
                 * (outmost) loop invariant : 2*f*w == n
                 */
                size_t w = 1, f = n >> 1;
                const Element *pow_ptr = pow_w_br.data();
                const Element *wp_ptr = pow_wp_br.data();

                DIF_reversed_core_firststeps (coeffs, w, f, pow_ptr, wp_ptr);

                for ( ; w < n; pow_ptr += f, wp_ptr += f, w <<= 1, f >>= 1)
                    DIF_reversed_core_one_step (coeffs, w, f, pow_ptr, wp_ptr);
            }

            void
            DIF (Element *coeffs) const {
                DIF_core (coeffs);
                reduce_coeffs_2p (coeffs);
            }

            void
            DIF_reversed (Element *coeffs) const {
                DIF_reversed_core (coeffs);
                reduce_coeffs_2p (coeffs);
            }

            /* DIT ************************************************************/
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
            /* no_simd version is defined for every type of Simd */
            void
            Butterfly_DIT (Element& A, Element& B, const Element& alpha,
                                                const Element& alphap) const {
                FFT_internals::sub_if_greater (A, p2); /* A -= 2*p if A >= 2p */
                Element tmp;
                this->fld->mul_precomp_b_without_reduction (tmp, B, alpha, alphap);
                B = A + (p2 - tmp);
                A += tmp;
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (simd_vect_t& A, simd_vect_t &B,
                                            const simd_vect_t& alpha,
                                            const simd_vect_t& alphap) const {
                simd_vect_t T1, T2, T3;
                T1 = SimdMod::reduce (A, this->P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdMod::mul_mod (B, alpha, this->P, alphap);
                /* A+B*alpha */
                A = Simd::add (T1, T2);
                /* A-B*alpha (computed as A+(2p-B*alpha)) */
                T3 = Simd::sub (this->P2, T2);
                B = Simd::add (T1, T3);
            }

            /* Same as above but with input of type Element*, so we can
             * interleave the store and some of the computation.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                                            const simd_vect_t& alpha,
                                            const simd_vect_t& alphap) const {
                simd_vect_t A, B, T1, T2, T3;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);

                T1 = SimdMod::reduce (A, this->P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdMod::mul_mod (B, alpha, this->P, alphap);
                /* A+B*alpha */
                A = Simd::add (T1, T2);
                Simd::store (Aptr, A);
                /* A-B*alpha (computed as A+(2p-B*alpha)) */
                T3 = Simd::sub (this->P2, T2);
                B = Simd::add (T1, T3);
                Simd::store (Bptr, B);
            }
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            Butterfly_DIT (Element *Aptr, Element *Bptr,
                                            const Element* alpha_ptr,
                                            const Element* alphap_ptr) const {
                simd_vect_t A, B, T1, T2, T3, alpha, alphap;
                A = Simd::load (Aptr);
                B = Simd::load (Bptr);
                alpha = Simd::load (alpha_ptr);
                alphap = Simd::load (alphap_ptr);

                T1 = SimdMod::reduce (A, this->P2); /* A - 2*p if A >= 2p */
                /* B*alpha */
                T2 = SimdMod::mul_mod (B, alpha, this->P, alphap);
                /* A+B*alpha */
                A = Simd::add (T1, T2);
                Simd::store (Aptr, A);
                /* A-B*alpha (computed as A+(2p-B*alpha)) */
                T3 = Simd::sub (this->P2, T2);
                B = Simd::add (T1, T3);
                Simd::store (Bptr, B);
            }

            /* See DIT_core for the meaning of the parameter */
            /* No simd version is defined for every type of Simd */
            void
            DIT_core_one_step_no_simd (Element *coeffs, const size_t w,
                                       const size_t f, const Element *pow,
                                                    const Element *wp) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    for (size_t j = 0; j < w; j += 1) {
                        Butterfly_DIT (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                                pow[j], wp[j]);
                    }
                }
            }

            void
            DIT_reversed_core_one_step_no_simd (Element *coeffs, const size_t w,
                                             const size_t f, const Element *pow,
                                                    const Element *wp) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    const Element alpha = pow[i];
                    const Element alphap = wp[i];
                    for (size_t j = 0; j < w; j += 1) {
                        Butterfly_DIT (Aptr[(i<<1)*w+j], Bptr[(i<<1)*w+j],
                                                                alpha, alphap);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_core_one_step (Element *coeffs, const size_t w, const size_t f,
                                                    const Element *pow,
                                                    const Element *wp) const {
                DIT_core_one_step_no_simd (coeffs, w, f, pow, wp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIT_core_one_step (Element *coeffs, const size_t w, const size_t f,
                                                    const Element *pow,
                                                    const Element *wp) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    for (size_t j = 0; j < w; j += Simd::vect_size) {
                        Butterfly_DIT (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j, pow+j,
                                                                         wp+j);
                    }
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_reversed_core_one_step (Element *coeffs, const size_t w,
                                        const size_t f, const Element *pow,
                                                    const Element *wp) const {
                DIT_reversed_core_one_step_no_simd (coeffs, w, f, pow, wp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size != 1>* = nullptr>
            void
            DIT_reversed_core_one_step (Element *coeffs, const size_t w,
                                        const size_t f, const Element *pow,
                                                    const Element *wp) const {
                Element *Aptr = coeffs;
                Element *Bptr = coeffs + w;
                for (size_t i = 0; i < f; i++) {
                    simd_vect_t alpha = Simd::set1 (pow[i]);
                    simd_vect_t alphap = Simd::set1 (wp[i]);
                    for (size_t j = 0; j < w; j += Simd::vect_size) {
                        Butterfly_DIT (Aptr+(i<<1)*w+j, Bptr+(i<<1)*w+j,
                                                                alpha, alphap);
                    }
                }
            }

            /* Firststeps perform the first log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = 1, f = n/2, pow = pow_w.data() + (n-2).
             * Note that for now, there is no special code for first steps with
             * SIMD, we just call the no_simd code.
             */
            void
            DIT_core_firststeps_fallback (Element *coeffs, size_t &w, size_t &f,
                                                    const Element *&pow,
                                                    const Element *&wp) const {
                size_t bound = (n < Simd::vect_size) ? n : Simd::vect_size;
                for ( ; w < bound; w <<= 1, f >>= 1, pow -= w, wp -= w)
                    DIT_core_one_step_no_simd (coeffs, w, f, pow, wp);
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 1>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                    const Element *&pow,
                                                    const Element *&wp) const {
                /* nothing to do for firststeps of S::vect_size == 1 */
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 2>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                    const Element *&pow,
                                                    const Element *&wp) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIT_core_firststeps_fallback (coeffs, w, f, pow, wp);
                } else {
                    for (size_t i = 0; i < f; i += 2, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C], V2 = [B D] */
                        V1 = Simd::set (coeffs[0], coeffs[2]);
                        V2 = Simd::set (coeffs[1], coeffs[3]);

                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, this->P2);
                        V2 = SimdMod::sub_mod (V1, V2, this->P2);

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
                    wp -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 4>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                    const Element *&pow,
                                                    const Element *&wp) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIT_core_firststeps_fallback (coeffs, w, f, pow, wp);
                } else {
                    f >>= 2;
                    w <<= 2;
                    pow -= 2;
                    wp -= 2;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[1], pow[1]);
                    simd_vect_t Wp = Simd::set (wp[0], wp[0], wp[1], wp[1]);
                    for (size_t i = 0; i < f; i++, coeffs += incr) {
                        simd_vect_t V1, V2, T;

                        /* V1 = [A C E G], V2 = [B D F H] */
                        V1 = Simd::set (coeffs[0], coeffs[2], coeffs[4],
                                                                    coeffs[6]);
                        V2 = Simd::set (coeffs[1], coeffs[3], coeffs[5],
                                                                    coeffs[7]);


                        /*** first step (special butterfly with mul by 1) *****/
                        T = SimdMod::add_mod (V1, V2, this->P2);
                        V2 = SimdMod::sub_mod (V1, V2, this->P2);

                        /* transform  T = [A C E G], V2 = [B D F H]
                         *      into V1 = [A E B F], V2 = [C G D H]
                         */
                        SimdExtra::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, Wp);

                        /* transform V1 = [A E B F], V2 = [C G D H]
                         *      into V1 = [A B C D], V2 = [E F G H] and store
                         */
                        SimdExtra::pack (V1, V2, V1, V2);
                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                    wp -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size == 8>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                    const Element *&pow,
                                                    const Element *&wp) const {
                const constexpr size_t incr = 2 * Simd::vect_size;
                if (n < incr) {
                    DIT_core_firststeps_fallback (coeffs, w, f, pow, wp);
                } else {
                    f >>= 3;
                    w <<= 3;
                    pow -= 2;
                    wp -= 2;
                    simd_vect_t W = Simd::set (pow[0], pow[0], pow[0], pow[0],
                                               pow[1], pow[1], pow[1], pow[1]);
                    simd_vect_t Wp = Simd::set (wp[0], wp[0], wp[0], wp[0],
                                                wp[1], wp[1], wp[1], wp[1]);
                    pow -= 4;
                    wp -= 4;
                    simd_vect_t W2 = Simd::set (pow[0], pow[0], pow[1], pow[1],
                                                pow[2], pow[2], pow[3], pow[3]);
                    simd_vect_t W2p = Simd::set (wp[0], wp[0], wp[1], wp[1],
                                                 wp[2], wp[2], wp[3], wp[3]);
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
                        T = SimdMod::add_mod (V1, V2, this->P2);
                        V2 = SimdMod::sub_mod (V1, V2, this->P2);

                        /* transform into
                         *      V1 = [A E I M B F J N], V2 = [C G K O D H L P]
                         */
                        SimdExtra::pack (V1, V2, T, V2);

                        /*** second step **************************************/
                        Butterfly_DIT (V1, V2, W, Wp);

                        /* transform into
                         *      V1 = [A I B J C K D L], V2 = [E M F N G O H P]
                         */
                        SimdExtra::pack (V1, V2, V1, V2);

                        /*** third step ***************************************/
                        Butterfly_DIT (V1, V2, W2, W2p);

                        /* transform into
                         *      V1 = [A B C D E F G H], V2 = [I J K L M N O P]
                         */
                        SimdExtra::pack (V1, V2, V1, V2);

                        Simd::store (coeffs, V1);
                        Simd::store (coeffs + Simd::vect_size, V2);
                    }
                    pow -= w;
                    wp -= w;
                }
            }

            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 16>* = nullptr>
            void
            DIT_core_firststeps (Element *coeffs, size_t &w, size_t &f,
                                                    const Element *&pow,
                                                    const Element *&wp) const {
                DIT_core_firststeps_fallback (coeffs, w, f, pow, wp);
            }

            /* Laststeps perform the last log2(Simd::vect_size) step(s).
             * w, f, and pow are passed to avoid recomputation, it is the
             * responsability of the caller to pass the correct value, i.e.,
             * w = Simd::vect_size/2, f = 4*n/Simd::vect_size,
             * pow = pow_w.data() - Simd::vect_size.
             */
            template <typename S=Simd,
                      FFT_utils::enable_if_t<S::vect_size >= 1>* = nullptr>
            void
            DIT_reversed_core_laststeps (Element *coeffs, size_t w, size_t f,
                                                 const Element *pow,
                                                 const Element *wp) const {
                for ( ; w > 0; f <<= 1, w >>= 1, pow -= f, wp -= f)
                    DIT_reversed_core_one_step_no_simd (coeffs, w, f, pow, wp);
            }

            /* core functions for DIT */
            void
            DIT_core (Element *coeffs) const {
                /* n : length of the array 'coeffs' (always of power of 2)
                 * f : number of families of butterflies
                 * w : width of butterflies
                 * (outmost) loop invariant : 2*f*w == n
                 */
                size_t w = 1, f = n >> 1;
                const Element *pow_ptr = pow_w.data() + (n-2);
                const Element *wp_ptr = pow_wp.data() + (n-2);

                DIT_core_firststeps (coeffs, w, f, pow_ptr, wp_ptr);

                for ( ; w < n; w <<= 1, f >>= 1, pow_ptr -= w, wp_ptr -= w)
                    DIT_core_one_step (coeffs, w, f, pow_ptr, wp_ptr);
            }

            void
            DIT_reversed_core (Element *coeffs) const {
                /* n : length of the array 'coeffs' (always of power of 2)
                 * f : number of families of butterflies
                 * w : width of butterflies
                 * (outmost) loop invariant : 2*f*w == n
                 */
                size_t w = n >> 1, f = 1;
                const Element *pow_ptr = pow_w_br.data() + (n-2);
                const Element *wp_ptr = pow_wp_br.data() + (n-2);

                for ( ; w >= Simd::vect_size; f <<= 1, w >>= 1, pow_ptr -= f,
                                                                wp_ptr -= f)
                    DIT_reversed_core_one_step (coeffs, w, f, pow_ptr, wp_ptr);

                DIT_reversed_core_laststeps (coeffs, w, f, pow_ptr, wp_ptr);
            }

            void
            DIT (Element *coeffs) const {
                DIT_core (coeffs);
                reduce_coeffs_4p (coeffs);
            }

            void
            DIT_reversed (Element *coeffs) const {
                DIT_reversed_core (coeffs);
                reduce_coeffs_4p (coeffs);
            }

            /* Utils **********************************************************/
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

    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
	// TODO : A rendre générique / Simd si on doit faire des précalculs dans des Simd::vect_t
	template <class Field>
	class FFT_init {
    private:
        using Element = typename Field::Element;
        using Compute_t = typename Field::Compute_t;
        using Residu_t = typename Field::Residu_t;

		/* TODO we should align on the max possible in order for the vect to be
         * used by all SIMD avalaible */
		typedef AlignedAllocator<Element, Alignment::DEFAULT> EltAlignedAllocator;
		typedef std::vector<Element, EltAlignedAllocator> aligned_vect;


	public:
		const Field                *fld;
		Residu_t              _pl, _dpl;
		size_t                       ln;
		uint64_t                      n;
        aligned_vect pow_w; /* Table of roots of unity.
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
        aligned_vect pow_w_br; /* Same as above with each subarray in bitreverse
                                * order.
                                */
        aligned_vect pow_wp; /* pow_wp[i] := precomp_b (pow_w[i])
                              * Used in integral modular multiplication.
                              */
        aligned_vect pow_wp_br; /* Same as above with each subarray in
                                 * bitreverse order.
                                 */

        /* ctor */
        FFT_init (const Field& F, size_t k, Element w = 0)
            : fld(&F), _pl(fld->characteristic()), _dpl(_pl << 1),
                ln(k), n(1UL << ln), pow_w(n-1), pow_w_br(n-1), pow_wp(n-1),
                pow_wp_br (n-1) {

            if (!k)
                throw LinBoxError ("FFT: k must be positive");

            if (!check_cardinality())
                throw LinBoxError ("FFT: 4p must be <= maxCardinality");

            if (w == 0)
                w = FFT_internals::compute_primitive_root (*fld, ln);
            else if (!FFT_internals::is_primitive_root (*fld, w, ln))
                throw LinBoxError ("FFT: w is not a primitive k-root of unity");

			init_powers (w);
		}

        /* getters */
        const Field &
        field () const
        {
            return *fld;
        }

        const Element &
        getRoot() const
        {
            if (ln == 1) /* if n=2^ln equals 2, return -1 */
                return fld->mOne;
            else
                return pow_w[1];
        }

        /* Return the inverse of w. We use the following facts:
         *    w^n = 1
         *    w^(n/2) = -1
         *    w^(-1) = w^(n-1) = w^(n/2+n/2-1) = -w^(n/2-1)
         *  - pow_w[i] = w^i for 0 <= i <= n/2-1
         */
        Element
        getInvRoot() const
        {
            Element invw;
            fld.init (invw);
            fld.neg (invw, pow_w[(n>>1)-1]);
            return invw;
        }

        /* Butterflies */
        void butterfly_DIF_mod2p (Element &A, Element &B, const Element &alpha,
                                  const Element &alphap);
        template <typename SIMD>
        typename enable_if<is_same<typename SIMD::scalar_t, Element>::value
                            && is_integral<Element>::value>::type
        butterfly_DIF_mod2p (Element *, Element *, const Element *,
                             const Element *, const typename SIMD::vect_t &,
                             const typename SIMD::vect_t &);

        
#if 0
        template <typename SIMD>
        typename enable_if<is_same<typename SIMD::scalar_t, Element>::value
                            && is_integral<Element>::value>::type
        butterfly_DIT_mod4p (Element *, Element *, const Element *,
                             const Element *, const typename SIMD::vect_t &,
                             const typename SIMD::vect_t &);
#endif

    private:
		template<typename T=Element>
		typename std::enable_if<std::is_integral<T>::value>::type
        init_powers (const Element & w) {

			size_t pos = 0;
			//uint64_t wi = 1;
			Element wi = 1;

			if (ln>0){
//				using simd=Simd<uint32_t>;
//				using vect_t =typename simd::vect_t;

				size_t tpts = 1 << (ln - 1);
				size_t i=0;
//				for( ;i<std::min(simd::vect_size+1, tpts);i++,pos++){
				// Precompute pow_wp[1] for faster mult by pow_w[1]
				for( ;i<std::min((size_t) 2, tpts);i++,pos++){
					pow_w[pos] = wi;

					// Fake conversion since precomp_b will be used as a Compute_t in mul_precomp_b
					Compute_t temp;
					fld->precomp_b(temp, wi); //(((Compute_t)wi*invp)>>(fld->_bitsizep));
					pow_wp[pos] = static_cast<Element>(temp);

					fld->mulin(wi, w);
				}
				/*
				vect_t wp_vect, Q_vect,BAR_vect,w_vect,pow_w_vect,pow_wp_vect, pl_vect;
				BAR_vect= simd::set1(BAR);
				wp_vect = simd::set1(pow_wp[simd::vect_size]);
				w_vect  = simd::set1(pow_w[simd::vect_size]);
				pl_vect = simd::set1(_pl);
				for (; i < ROUND_DOWN(tpts,simd::vect_size);
					 i+=simd::vect_size,pos+=simd::vect_size) {
					pow_w_vect  = simd::loadu((int32_t*)pow_w.data()+pos-simd::vect_size);
					Q_vect=simd::mulhi(pow_w_vect,wp_vect);
					pow_w_vect = simd::sub(simd::mullo(pow_w_vect,w_vect),simd::mullo(Q_vect,pl_vect));
					pow_w_vect=simd::sub(pow_w_vect, simd::vandnot(simd::greater(pow_w_vect,pl_vect),pl_vect));
					simd::storeu((int32_t*)pow_w.data()+pos,pow_w_vect);
					pow_wp_vect= simd::mulhi(simd::sll(pow_w_vect,32-_logp),BAR_vect);
					simd::storeu((int32_t*)pow_wp.data()+pos,pow_wp_vect);
				}
				*/
				// Use pow_wp[1] for speed-up mult by pow_w[1]
				for( ;i<tpts;i++,pos++){
					pow_w[pos] = wi;

					// Fake conversion since precomp_b will be used as a Compute_t in mul_precomp_b
					Compute_t temp;
					fld->precomp_b(temp, wi); //(((Compute_t)wi*invp)>>(fld->_bitsizep));
					pow_wp[pos] = static_cast<Element>(temp);

					fld->mul_precomp_b(wi, wi, w, static_cast<Compute_t>(pow_wp[1]));
				}

				// Other pow_w elements can be read from previously computed pow_w
				for(size_t k=2;k<=tpts;k<<=1)
					for(size_t i=0;i<tpts;i+=k,pos++){
						pow_w[pos]  = pow_w[i];
						pow_wp[pos] = pow_wp[i];
					}

//				std::cout << "Check precomputations : pow_w, pow_wp \n";
//				std::cout << "[";
//				for (size_t i=0; i < tpts; i++) std::cout << pow_w[i] << ", ";
//				std::cout << "]\n";
//				std::cout << "[";
//				for (size_t i=0; i < tpts; i++) std::cout << pow_wp[i] << ", ";
//				std::cout << "]\n\n";
			}

            /* init powers in bitreverse order */
            size_t l = ln, len = n >> 1, base_idx = 0;
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

		template<typename T=Element>
		typename std::enable_if<std::is_floating_point<T>::value>::type
        init_powers (const Element & w) {

			size_t pos = 0;
			//uint64_t wi = 1;
			Element wi = 1;

			if (ln>0){
				size_t tpts = 1 << (ln - 1);

				for(size_t i=0; i<tpts;i++,pos++){
					pow_w[pos] = wi;
					fld->mulin(wi,w);
				}

				// Other pow_w elements can be read from previously computed pow_w
				for(size_t k=2;k<=tpts;k<<=1)
					for(size_t i=0;i<tpts;i+=k,pos++){
						pow_w[pos]  = pow_w[i];
					}

			}

            /* init powers in bitreverse order */
            size_t l = ln, len = n >> 1, base_idx = 0;
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

        template<typename T=Element>
        typename std::enable_if<std::is_integral<T>::value, bool>::type
        check_cardinality ()
        {
            /* for integral types, we use Harvey's butterflies that do
             * computation with residues <= 4*p, so we need
             * 4*p <= maxCardinality
             */
            return _pl <= (field().maxCardinality() >> 2);
        }
        template<typename T=Element>
        typename std::enable_if<std::is_floating_point<T>::value, bool>::type
        check_cardinality ()
        {
            /* for floating point types, the above restriction is always
             * satisfied.
             */
            return true;
        }

	};

template <class Field>
template <typename SIMD>
typename enable_if<is_same<typename SIMD::scalar_t, typename Field::Element>::value
                   && is_integral<typename Field::Element>::value>::type
FFT_init<Field>::butterfly_DIF_mod2p (Element *a, Element *b,
                                      const Element *alpha,
                                      const Element *alphap,
                                      const typename SIMD::vect_t & p,
                                      const typename SIMD::vect_t & p2)
{
}

/*
template <class Field>
template <>
void
FFT_init<Field>::butterfly_DIT_mod4p<NoSimd<typename Field::Element>>
    (Element *a, Element *b, const Element *alpha, const Element *alphap,
     const typename NoSimd<typename Field::Element>::vect_t & p,
     const typename NoSimd<typename Field::Element>::vect_t & p2)
{
}
*/

}
#endif // __LINBOX_polynomial_fft_init_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
