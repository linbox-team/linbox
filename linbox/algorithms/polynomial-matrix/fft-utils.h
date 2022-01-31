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

#ifndef __LINBOX_fft_utils_H
#define __LINBOX_fft_utils_H

#include <type_traits>

namespace LinBox {

    /* Forward declaration of FFT_base class */
    template<typename Field, typename Simd = Simd<typename Field::Element>,
                             typename Enable = void>
    class FFT_base;

    /* Forward declaration of FFT_multi_base class */
    template<typename Field, typename Simd = Simd<typename Field::Element>,
                             typename Enable = void>
    class FFT_multi_base;

    /* Helper to define both Simd and NoSimd methods */
    template<size_t VectSize>
    struct FFTSimdHelper { size_t size = VectSize; };

    /* Some utils functions */
    namespace FFT_utils {

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
}
#endif /* __LINBOX_fft_utils_H */

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
