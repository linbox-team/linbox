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


#ifndef __LINBOX_fft_H
#define __LINBOX_fft_H


#include "linbox/linbox-config.h"
#include "linbox/util/error.h"

#include "fflas-ffpack/fflas/fflas_simd.h"

/* This file contains the specialization of FFT_base for Modular based on
 * floating point types */
#include "fft-floating.inl"

/* This file contains the specialization of FFT_base for Modular based on
 * integral types */
#include "fft-integral.inl"

namespace LinBox {

    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
    template <typename Field, typename Simd= Simd<typename Field::Element> >
    class FFT : public FFT_base<Field, Simd>
    {
        private:
            using Element = typename Field::Element;

        public:
            FFT (const Field& F, size_t k, Element w = 0)
                        : FFT_base<Field, Simd>(F, k, check_args (F, k, w)) {
            }

            /******************************************************************/
            /* main functions: FFT_direct and FFT_inverse *********************/
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
                this->DIF (coeffs); /* or DIT_reversed */
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
                this->DIT (coeffs); /* or DIF_reversed */
            }


            /******************************************************************/
            /* getters ********************************************************/
            /******************************************************************/
            const Field &
            field () const {
                return *(this->fld);
            }

            size_t
            size () const {
                return this->n;
            }

            size_t
            log2_size () const {
                return this->l2n;
            }

            const Element &
            root() const
            {
                if (this->l2n == 1) /* if n=2^ln equals 2, return -1 */
                    return this->fld->mOne;
                else
                    return this->pow_w[1];
            }

            const Element
            invroot() const
            {
                if (this->l2n == 1) /* if n=2^ln equals 2, return -1 */
                    return this->fld->mOne;
                else { /* w^-1 = w^(n/2) * w^(n/2-1) = -1 * w^(n/2-1) */
                    Element r;
                    this->fld->assign (r, this->pow_w[(this->n>>1)-1]);
                    this->fld->negin (r);
                    return r;
                }
            }

        private:
            static Element
            check_args (const Field& fld, size_t k, Element w) {
                if (!k)
                    throw LinBoxError ("FFT: k must be positive");

                if (w == 0)
                    return FFT_utils::compute_primitive_root (fld, k);
                else if (FFT_utils::is_primitive_root (fld, w, k))
                    return w;
                else
                    throw LinBoxError ("w is not a primitive k-root of unity");
            }

    };

    /**************************************************************************/
    /**************************************************************************/
    /**************************************************************************/
    template <typename Field, typename Simd= Simd<typename Field::Element> >
    class FFT_multi : public FFT_multi_base<Field, Simd>
    {
        private:
            using Element = typename Field::Element;

        public:
            FFT_multi (const Field& F, size_t k, Element w = 0)
                : FFT_multi_base<Field, Simd>(F, k, check_args(F, k, w)) {
            }

            /******************************************************************/
            /* main functions: FFT_direct and FFT_inverse *********************/
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
            FFT_direct (Element *coeffs, size_t stride) const {
                this->DIF (coeffs, stride); /* or DIT_reversed */
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
            FFT_inverse (Element *coeffs, size_t stride) const {
                this->DIT (coeffs, stride); /* or DIF_reversed */
            }


            /******************************************************************/
            /* getters ********************************************************/
            /******************************************************************/
            const Field &
            field () const {
                return *(this->fld);
            }

            size_t
            size () const {
                return this->n;
            }

            size_t
            log2_size () const {
                return this->l2n;
            }

            const Element &
            root() const
            {
                if (this->l2n == 1) /* if n=2^ln equals 2, return -1 */
                    return this->fld->mOne;
                else
                    return this->pow_w[1];
            }

            const Element
            invroot() const
            {
                if (this->l2n == 1) /* if n=2^ln equals 2, return -1 */
                    return this->fld->mOne;
                else { /* w^-1 = w^(n/2) * w^(n/2-1) = -1 * w^(n/2-1) */
                    Element r;
                    this->fld->assign (r, this->pow_w[(this->n>>1)-1]);
                    this->fld->negin (r);
                    return r;
                }
            }

        private:
            static Element
            check_args (const Field& fld, size_t k, Element w) {
                if (!k)
                    throw LinBoxError ("FFT_multi: k must be positive");

                if (w == 0)
                    return FFT_utils::compute_primitive_root (fld, k);
                else if (FFT_utils::is_primitive_root (fld, w, k))
                    return w;
                else
                    throw LinBoxError ("w is not a primitive k-root of unity");
            }

    };
}
#endif // __LINBOX_fft_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
