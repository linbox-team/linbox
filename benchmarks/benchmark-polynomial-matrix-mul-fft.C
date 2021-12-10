/*
 * Copyright (C) 2020 Cyril Bouvier, Pascal Giorgi
 *
 * Written by Cyril Bouvier <cyril.bouvier@lirmm.fr>
 *            Pascal Giorgi <pascal.giorgi@lirmm.fr>
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




#include <iostream>
#include <iomanip>
using namespace std;

//#define FFT_PROFILER
#include <linbox/ring/modular.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/randiter/random-fftprime.h>
#include <givaro/zring.h>
#include <recint/rint.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/matrix/polynomial-matrix.h>
#include <linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h>
#include <linbox/algorithms/polynomial-matrix/matpoly-mult-naive.h>
#include <linbox/algorithms/polynomial-matrix/matpoly-mult-kara.h>

#include <givaro/modular.h>
#include <givaro/modular-extended.h>
#include <fflas-ffpack/utils/args-parser.h>

using namespace std;
using namespace LinBox;
using Givaro::Modular;
using Givaro::ModularExtended;

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
REGISTER_TYPE_NAME(uint128_t);
REGISTER_TYPE_NAME(Modular);
REGISTER_TYPE_NAME(ModularExtended);



using namespace LinBox;



template<typename PolMatType, typename PolMatMulDomain>
double bench_one (const PolMatMulDomain &PMMD, unsigned int m, unsigned int n,
                  unsigned int k, unsigned int d, unsigned long seed) {
    Timer chrono;
    size_t cnt;
    double time;
    typename PolMatType::Field::RandIter G(PMMD.field(), seed);
    
    PolMatType M(PMMD.field(),m,n,d), N(PMMD.field(),n,k,d);
    PolMatType R(PMMD.field(),m,k,2*d-1);

    M.random(G);
    N.random(G);
    
    chrono.start();
    for (cnt = 0; chrono.realElapsedTime() < 1 ; cnt++)
        PMMD.mul (R, M, N);
    time = chrono.userElapsedTime()/cnt; /* time per iteration */

    return time;
}

/* Bench multiplication via FFT on polynomial matrices with coefficients in
 * ModImplem<Elt, C...> (i.e., Modular<Elt, C> or ModularExtended<Elt>) with a
 * random prime p with the required number of bits and such that 2^k divides p-1
 */
template<template<typename, typename...> class ModImplem, typename Elt, typename... C>
void bench_one_modular_implem_fft (uint64_t bits, unsigned int m,
                                   unsigned int n, unsigned int k,
                                   unsigned int d, unsigned long seed)
{
    typedef ModImplem<Elt, C...> Field;
    double time;

    /* First, we check if this Modular implem can handle this many bits */
    bool ok;
    typename Field::Residu_t s = 1;
    for (uint64_t i = 0; i+1 < bits; i++, s<<=1); /* at the end s=2^(bits-1) */
    if (s == 0)
        ok = false;
    else {
        s <<= 1; s--; /* now s=2^bits-1 */
        ok = (s <= Field::maxCardinality());
    }
    if (!ok) {
        cout << endl << "# Skipping bench with " << TypeName<ModImplem>();
        cout << "<" << TypeName<Elt>();
        if (sizeof...(C) > 0)
            cout << ", " << TypeName<C...>();
        cout << ">, bits = " << bits << " is too large" << endl;
        return;
    }

    size_t lpts = 0;
    size_t pts = 1; while (pts <= 2*d-2) { pts<<=1; lpts++; }
    integer p;
    RandomFFTPrime::seeding (seed);
    if (!RandomFFTPrime::randomPrime (p, integer(1)<<bits, lpts)) {
        cout << endl << "# Skipping bench with " << TypeName<ModImplem>();
        cout << "<" << TypeName<Elt>();
        if (sizeof...(C) > 0)
            cout << ", " << TypeName<C...>();
        cout << ">, RandomFFTPrime::randomPrime failed to find a prime" << endl;
        return;
    }
    Field GFp ((Elt) p);

    cout << endl << string (80, '*') << endl;
    cout << "Bench polynomial matrix multiplication with "
         << TypeName<ModImplem>() << "<" <<TypeName<Elt>();
    if (sizeof...(C) > 0)
        cout << ", " << TypeName<C...>();
    cout << ">, p=" << GFp.cardinality() << " (" << bits << " bits, ";
    cout << "n=2^" << k << " divides p-1)" << endl;

    typedef PolynomialMatrixFFTPrimeMulDomain<Field> PMMDType;
    PMMDType PMMD_fft (GFp);

    typedef PolynomialMatrix<Field, PMType::polfirst> MatrixP;
    typedef PolynomialMatrix<Field, PMType::matfirst> PMatrix;
    typedef PolynomialMatrix<Field, PMType::matrowfirst> PMatrixP;

    /* polfirst */
    time = bench_one<MatrixP> (PMMD_fft, m, n, k, d, seed);
    cout << "  polfirst " << string (80-25, ' ');
    cout.precision(2); cout.width(10); cout<< scientific << time << " s";
    cout << endl;

    /* matfirst */
    time = bench_one<PMatrix> (PMMD_fft, m, n, k, d, seed);
    cout << "  matfirst " << string (80-25, ' ');
    cout.precision(2); cout.width(10); cout<< scientific << time << " s";
    cout << endl;

    /* matrowfirst */
    /* TODO: does not compile */
    time = bench_one<PMatrixP> (PMMD_fft, m, n, k, d, seed);
    cout << "  matrowfirst " << string (80-25-3, ' ');
    cout.precision(2); cout.width(10); cout<< scientific << time << " s";
    cout << endl;
}

/******************************************************************************/
/************************************ main ************************************/
/******************************************************************************/
int main (int argc, char* argv[]) {
    unsigned long bits = 23;
    unsigned long seed = time (NULL);
    unsigned long m = 10;
    unsigned long n = 30;
    unsigned long k = 20;
    unsigned long d = 2000;

    Argument args[] = {
        { 'b', "-b nbits", "number of bits of prime.", TYPE_INT, &bits },
        { 'm', "-m m", "number of rows of first matrix.", TYPE_INT, &m },
        { 'n', "-n n", "number of columns of first matrix and rows of second matrix.", TYPE_INT, &n },
        { 'k', "-k k", "number of columns of second matrix.", TYPE_INT, &k },
        { 'd', "-d d", "strict bound on degree of matrices.", TYPE_INT, &d },
        { 's', "-s seed", "set the seed.", TYPE_INT, &seed },
        END_OF_ARGUMENTS
    };

    parseArguments (argc, argv, args);

    cout << "# command: ";
    FFLAS::writeCommandString (cout, args, "benchmark-polynomial-matrix-mul-fft") << endl;

    if (d >> (bits-1)) {
        cerr << "Error, d=" << d << " must be smaller than 2^nbits-1";
        cerr << endl;
        return 1;
    }

    /* Bench with Modular<float, double> */
    bench_one_modular_implem_fft<Modular, float, double> (bits, m, n, k, d, seed);

    /* Bench with Modular<double, double> */
    bench_one_modular_implem_fft<Modular, double> (bits, m, n, k, d, seed);

    /* Bench with ModularExtended<double> */
    /* TODO: does not compile */
    //bench_one_modular_implem_fft<ModularExtended, double> (bits, m, n, k, d, seed);

    /* Bench with Modular<uint16_t,uint32_t> */
    bench_one_modular_implem_fft<Modular, uint16_t, uint32_t> (bits, m, n, k, d, seed);

    /* Bench with Modular<uint32_t> */
    bench_one_modular_implem_fft<Modular, uint32_t> (bits, m, n, k, d, seed);

    /* Bench with Modular<uint32_t, uint64_t> */
    bench_one_modular_implem_fft<Modular, uint32_t, uint64_t> (bits, m, n, k, d, seed);

    /* Bench with Modular<uint64_t> */
    bench_one_modular_implem_fft<Modular, uint64_t> (bits, m, n, k, d, seed);

#ifdef __FFLASFFPACK_HAVE_INT128
    /* Bench with Modular<uint64_t,uint128_t> */
    bench_one_modular_implem_fft<Modular, uint64_t, uint128_t> (bits, m, n, k, d, seed);
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
