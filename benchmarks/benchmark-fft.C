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

#include "linbox/algorithms/polynomial-matrix/fft.h"
#include "linbox/randiter/random-fftprime.h"

#include <givaro/modular.h>
#include <givaro/modular-extended.h>
#include <fflas-ffpack/utils/args-parser.h>

#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace LinBox;
using Givaro::Modular;
using Givaro::ModularExtended;
using FFLAS::parseArguments;

/* For pretty printing type */
template<typename ...> const char *TypeName();
template<template <typename ...> class> const char *TypeName();

#define REGISTER_TYPE_NAME(type) \
    template<> const char *TypeName<type>(){return #type;}

REGISTER_TYPE_NAME(float);
REGISTER_TYPE_NAME(double);
REGISTER_TYPE_NAME(uint16_t);
REGISTER_TYPE_NAME(uint32_t);
REGISTER_TYPE_NAME(uint64_t);
#ifdef __FFLASFFPACK_HAVE_INT128
REGISTER_TYPE_NAME(uint128_t);
#endif
REGISTER_TYPE_NAME(Modular);
REGISTER_TYPE_NAME(ModularExtended);

/******************************************************************************/
template<typename Field, typename Simd>
struct BenchFFT {
    static const size_t min_run = 4; /* do at least this number of run of fft */

    FFT<Field, Simd> fft;

    BenchFFT (const Field &F, size_t k, unsigned long seed) : fft (F, k) {

        typename Simd::aligned_vector in (fft.size()), v(fft.size());
        Timer chrono;
        double time;
        size_t cnt;

        /* Generate random input */
        typename Field::RandIter Gen (F, 0, seed+k); /* 0 has no meaning here */
        for (auto elt = in.begin(); elt < in.end(); elt++)
            Gen.random (*elt);

        /* DIF */
        v = in;
        chrono.start();
        for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
            fft.DIF (v.data());
        time = chrono.userElapsedTime()/cnt; /* time per iteration */
        print_result_line ("DIF", time);

        /* DIT_reversed */
        v = in;
        chrono.start();
        for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
            fft.DIT_reversed (v.data());
        time = chrono.userElapsedTime()/cnt; /* time per iteration */
        print_result_line ("DIT_reversed", time);

        /* FFT direct */
        v = in;
        chrono.start();
        for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
            fft.FFT_direct (v.data());
        time = chrono.userElapsedTime()/cnt; /* time per iteration */
        print_result_line ("FFT_direct", time);

        cout << endl;

        /* DIT */
        v = in;
        chrono.start();
        for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
            fft.DIT (v.data());
        time = chrono.userElapsedTime()/cnt; /* time per iteration */
        print_result_line ("DIT", time);

        /* DIF_reversed */
        v = in;
        chrono.start();
        for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
            fft.DIF_reversed (v.data());
        time = chrono.userElapsedTime()/cnt; /* time per iteration */
        print_result_line ("DIF_reversed", time);

        /* FFT inverse (without div) */
        v = in;
        chrono.start();
        for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
            fft.FFT_inverse (v.data());
        time = chrono.userElapsedTime()/cnt; /* time per iteration */
        print_result_line ("FFT_inverse", time);
    }

    void
    print_result_line (const char *name, double time) {
        /* Miops = #operation (~3/2 n log n) / time / 1e6 */
        size_t l2n = fft.log2_size();
        double Miops = 17 * (l2n<<(l2n-1)) / (1e6 * time);
        size_t t = strlen(name) + Simd::type_string().size() + 15 + 35;
        cout << "  " << string (80-t, ' ') << name << "<";
        cout << Simd::type_string() << "> " << string (15, '.');
        cout.precision(2); cout.width(10); cout<< scientific << time << " s, ";
        cout.precision(2); cout.width(10); cout<<fixed<<Miops << " Miops";
        cout << endl;
    }
};

/* Bench FFT on polynomial with coefficients in ModImplem<Elt, C...> (i.e.,
 * Modular<Elt, C> or ModularExtended<Elt>) with a random prime p with the
 * required number of bits and such that 2^k divides p-1
 */
template<template<typename, typename...> class ModImplem, typename Elt, typename... C>
void bench_one_modular_implem (uint64_t bits, size_t k, unsigned long seed)
{
    /* First, we check if this Modular implem can handle this many bits */
    bool ok;
    typename ModImplem<Elt, C...>::Residu_t s = 1;
    for (uint64_t i = 0; i+1 < bits; i++, s<<=1); /* at the end s=2^(bits-1) */
    if (s == 0)
        ok = false;
    else {
        s <<= 1; s--; /* now s=2^bits-1 */
        ok = (s <= ModImplem<Elt, C...>::maxCardinality());
    }
    if (!ok) {
        cout << endl << "# Skipping bench with " << TypeName<ModImplem>();
        cout << "<" << TypeName<Elt>();
        if (sizeof...(C) > 0)
            cout << ", " << TypeName<C...>();
        cout << ">, bits = " << bits << " is too large" << endl;
        return;
    }

    integer p;
    RandomFFTPrime::seeding (seed);
    if (!RandomFFTPrime::randomPrime (p, integer(1)<<bits, k))
        throw LinboxError ("RandomFFTPrime::randomPrime failed");
    ModImplem<Elt, C...> GFp ((Elt) p);

    cout << endl << string (80, '*') << endl;
    cout << "Bench FFT with " << TypeName<ModImplem>() << "<" <<TypeName<Elt>();
    if (sizeof...(C) > 0)
        cout << ", " << TypeName<C...>();
    cout << ">, p=" << GFp.cardinality() << " (" << bits << " bits, ";
    cout << "n=2^" << k << " divides p-1)" << endl;

    BenchFFT<ModImplem<Elt, C...>, NoSimd<Elt>> BenchNoSimd (GFp, k, seed);

    /* Simd128 */
#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
    cout << string (15, ' ') << string (50, '-') << string (15, ' ') << endl;
    BenchFFT<ModImplem<Elt, C...>, Simd128<Elt>> BenchSimd128 (GFp, k, seed);
#endif

    /* Simd256 */
#if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
    cout << string (15, ' ') << string (50, '-') << string (15, ' ') << endl;
    BenchFFT<ModImplem<Elt, C...>, Simd256<Elt>> BenchSimd256 (GFp, k, seed);
#endif

    // FIXME Simd512: Which macro should be used to discriminate against
    // simd512 for both floating and integral type ?
}

/******************************************************************************/
/************************************ main ************************************/
/******************************************************************************/
int main (int argc, char *argv[]) {
    unsigned long bits = 27;
    unsigned long k = 16;
    unsigned long seed = time (NULL);

    Argument args[] = {
        { 'b', "-b nbits", "number of bits of prime.", TYPE_INT, &bits },
        { 'k', "-k k", "bench FFT for n=2^k.", TYPE_INT, &k },
        { 's', "-s seed", "set the seed.", TYPE_INT, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments (argc, argv, args);

    cout << "# command: ";
    FFLAS::writeCommandString (cout, args, "benchmark-fft-new") << endl;

    if (k >= bits) {
        cerr << "Error, k=" << k << " must be smaller than nbits=" << bits;
        cerr << endl;
        return 1;
    }

    /* Bench with Modular<double, double> */
    bench_one_modular_implem<Modular, float, double> (bits, k, seed);

    /* Bench with Modular<double, double> */
    bench_one_modular_implem<Modular, double> (bits, k, seed);

    /* Bench with ModularExtended<double> */
    bench_one_modular_implem<ModularExtended, double> (bits, k, seed);

    /* Bench with Modular<uint16_t,uint32_t> */
    bench_one_modular_implem<Modular, uint16_t, uint32_t> (bits, k, seed);

    /* Bench with Modular<uint32_t, uint64_t> */
    bench_one_modular_implem<Modular, uint32_t, uint64_t> (bits, k, seed);

#ifdef __FFLASFFPACK_HAVE_INT128
    /* Bench with Modular<uint64_t,uint128_t> */
    bench_one_modular_implem<Modular, uint64_t, uint128_t> (bits, k, seed);
#endif

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
