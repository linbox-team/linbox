/*
 * Copyright (C) 2013  Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@lirmm.fr>
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

#include "linbox/linbox-config.h"

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-algorithms.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/ring/modular.h"

#include <givaro/modular-extended.h>

#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace LinBox;
using Givaro::Modular;
using Givaro::ModularExtended;

/* For pretty printing type */
template<typename...> const char *TypeName();
template<template <typename...> class> const char *TypeName();

#define REGISTER_TYPE_NAME(type) \
    template<> const char *TypeName<type>(){return #type;}

REGISTER_TYPE_NAME(double);
REGISTER_TYPE_NAME(uint16_t);
REGISTER_TYPE_NAME(uint32_t);
REGISTER_TYPE_NAME(uint64_t);
REGISTER_TYPE_NAME(uint128_t);
REGISTER_TYPE_NAME(NoSimd);
REGISTER_TYPE_NAME(Simd128);
REGISTER_TYPE_NAME(Simd256);
REGISTER_TYPE_NAME(Modular);
REGISTER_TYPE_NAME(ModularExtended);

/******************************************************************************/
template<typename Field>
struct Checker {
    typedef typename Field::Element Elt;
    typedef AlignedAllocator<Elt, Alignment::DEFAULT> _Alignement;
    typedef typename std::vector<Elt, _Alignement> EltVector;

    const Field& _F;
    size_t _k;
    size_t _n;
    FFT_init<Field> _fft;

    /* ctor */
    Checker (const Field& F, size_t k) : _F(F), _k(k), _n(1<<k), _fft(F,k) {
    }

    /* compute r = P evaluated at x using Horner method */
    void horner_eval (Elt& r, EltVector& P, Elt& x) {
        _F.assign (r, _F.zero);
        for (auto ptr = P.rbegin(); ptr != P.rend(); ptr++ ) {
            _F.mulin (r, x);
            _F.addin (r, *ptr);
        }
    }

    /* return i with its first _k bits reversed */
    size_t bitreversed (size_t i) {
        size_t r = 0;
        for (size_t j = 0; j < _k; j++, i>>=1)
        {
            r <<= 1;
            r |= i & 1;
        }
        return r;
    }

    /* check DIF and DIT */
    template<template <typename Elt> class SimdType>
    bool actual_check (const EltVector& in, const EltVector& in_br,
                            const EltVector& out, const EltVector& out_br) {
        FFT_algorithms<Field, SimdType<Elt>> fft_algo (_fft);
        EltVector v(_n);
        string s;
        s.append ("<"); s.append (TypeName<SimdType>()); s.append ("> ");
        s.append (string (80-(s.size()+16), '.'));

        /* DIF : natural order => bitreversed order */
        v = in;
        fft_algo.DIF (v.data());
        bool bf = equal (v.begin(), v.end(), out_br.begin());
        cout << "     DIF" << s << (bf ? " success" : " failure") << endl;

        /* DIT : bitreversed order => natural order */
        v = in_br;
        fft_algo.DIT (v.data());
        bool bt = equal (v.begin(), v.end(), out.begin());
        cout << "     DIT" << s << (bt ? " success" : " failure") << endl;

        bool b = bt & bf;
        return b;
    }

    /* draw random vector and check DIF and DIT for all available SIMD implem */
    bool check (unsigned long seed) {
        bool passed = true;
        EltVector in(_n), in_br(_n), out(_n), out_br(_n);

        /* Generate random input */
        typename Field::RandIter Gen (_F, 0, seed+_k); /*0 has no meaning here*/
        for (auto v = in.begin(); v < in.end(); v++)
            Gen.random (*v);

        /* Compute the out vector using a naive polynomial evaluation and set
         * the bitreversed version of in and out */
        Elt x(1);
        for (size_t i = 0; i < _n; _F.mulin (x, _fft._w), i++) {
            size_t i_br = bitreversed (i);
            in_br[i_br] = in[i];
            horner_eval (out[i], in, x);
            out_br[i_br] = out[i];
        }

        /* NoSimd */
        passed &= actual_check<NoSimd> (in, in_br, out, out_br);

        /* Simd128 */
#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
        if (Simd128<Elt>::vect_size == 4 || Simd128<Elt>::vect_size == 8)
            passed &= actual_check<Simd128> (in, in_br, out, out_br);
#endif

        /* Simd256 */
#if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
        if (Simd256<Elt>::vect_size == 4 || Simd256<Elt>::vect_size == 8)
            passed &= actual_check<Simd256> (in, in_br, out, out_br);
#endif

        return passed;
    }
};

/* Test FFT on polynomial with coefficients in Modular<T1,T2> with a random
 * prime p with the required number of bits and such that 2^k divides p-1
 */
template<template<typename, typename...> class ModImplem, typename Elt, typename... C>
bool test_one_modular_implem (uint64_t bits, size_t k, unsigned long seed)
{

    RandomFFTPrime RandomGen (integer(1)<<bits, seed);
    integer p = RandomGen.randomPrime (k);
    ModImplem<Elt, C...> GFp ((Elt) p);

    cout << endl << string (80, '*') << endl;
    cout << "Test FFT with " << TypeName<ModImplem>() << "<" << TypeName<Elt>();
    if (sizeof...(C) > 0)
        cout << ", " << TypeName<C...>();
    cout << ">, p=" << GFp.cardinality() << " (" << bits << " bits, 2^" << k;
    cout << " divides p-1)" << endl;

    bool b = true;
    for (size_t kc = 1; kc <= k; kc++) {
        cout << "** with n=2^" << kc << endl;
        Checker<ModImplem<Elt, C...>> Test(GFp, kc);
        b &= Test.check (seed);
    }

    return b;
}

/******************************************************************************/
/************************************ main ************************************/
/******************************************************************************/
int main (int argc, char *argv[]) {
    bool pass = true;
    unsigned long seed = time (NULL);

    Argument args[] = {
        { 's', "-s seed", "set the seed.", TYPE_INT, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments (argc, argv, args);

    cout << "# To rerun this test: test-fft-new -s " << seed << endl;
    cout << "# seed = " << seed << endl;

    /* Test with Modular<double>, and a 22-bit prime and k=10 */
    pass &= test_one_modular_implem<Modular,double> (22, 10, seed);

    /* Test with ModularExtended<double>, and a 51-bit prime and k=11 */
    pass &= test_one_modular_implem<ModularExtended,double> (51, 11, seed);

    /* Test with Modular<uint16_t,uint32_t>, and 11-bit prime and k=7 */
    pass &= test_one_modular_implem<Modular,uint16_t,uint32_t> (11, 7, seed);

    /* Test with Modular<uint32_t, uint64_t>, and 27-bit prime and k=10 */
    pass &= test_one_modular_implem<Modular,uint32_t,uint64_t> (27, 10, seed);

#ifdef __FFLASFFPACK_HAVE_INT128
    /* Test with Modular<uint64_t,uint128_t>, and a 59-bit prime and k=11 */
    pass &= test_one_modular_implem<Modular,uint64_t,uint128_t> (59, 11, seed);
#endif

    cout << endl << "Test " << (pass ? "passed" : "failed") << endl;
    return pass ? 0 : 1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
