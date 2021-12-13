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
#include "givaro/modular.h"
#include "givaro/modular-extended.h"

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
        template <class T=vect_t, enable_if_floating_point_t<T>* = nullptr,
                  enable_if_t<std::is_same<Field, typename Givaro::ModularExtended<Element>>::value
                              || sizeof(Element) != sizeof(Compute_t), T>* =nullptr>
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
        /* mul mod for floating type */
        template <class T=vect_t, enable_if_floating_point_t<T>* = nullptr,
                  enable_if_compute_same_size_t<T>* = nullptr,
                  enable_if_t<!std::is_same<Field, typename Givaro::ModularExtended<Element>>::value, T>* = nullptr>
        static inline T
        mul_mod (const vect_t& x, const vect_t& y, const vect_t& p,
                 const vect_t& u) {
            // u = 1/p
            vect_t h = Simd::mul(x,y);
            //vect_t l = Simd::fmsub(h,x,y); // Beware of the order!
            vect_t b = Simd::mul(h,u);
            vect_t c = Simd::floor(b);
            vect_t d = Simd::fnmadd(h,c,p); // Beware of the order!
            //vect_t g = Simd::add(d,l);
            vect_t t = Simd::sub(d,p);
            d = Simd::blendv(t,d,t);
            t = Simd::add(d,p);
            return Simd::blendv(d,t,d);
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
