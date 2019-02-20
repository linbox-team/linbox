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

#include <functional>
#include <iostream>
#include <string>
#include <vector>

using namespace std;
using namespace LinBox;

/* For pretty printing type */
template<typename T> const char *TypeName();
template<template <typename S> typename T> const char *TypeName();

#define REGISTER_TYPE_NAME(type) \
	template<> const char *TypeName<type>(){return #type;}

REGISTER_TYPE_NAME(double);
REGISTER_TYPE_NAME(uint16_t);
REGISTER_TYPE_NAME(uint32_t);
REGISTER_TYPE_NAME(uint64_t);
REGISTER_TYPE_NAME(uint128_t);

/******************************************************************************/
template<typename Field>
struct Benchs {
	typedef typename Field::Element Elt;
	typedef AlignedAllocator<Elt, Alignment::DEFAULT> _Alignement;
	typedef typename std::vector<Elt, _Alignement> EltVector;

	static const size_t min_run = 4; /* do at least this number of run of fft */
	const Field& _F;
	size_t _k;
	size_t _n;
	FFT_init<Field> _fft;

	/* ctor */
	Benchs (const Field& F, size_t k) : _F(F), _k(k), _n(1<<k), _fft(F,k) {
	}

	/* bench DIF and DIT */
	template<template <typename Elt> typename SimdType>
	void actual_bench (const EltVector& in) {
		Timer chrono;
		double time, Miops;
		size_t cnt;
		FFT_algorithms<Field, SimdType<Elt>> fft_algo (_fft);
		EltVector v(_n);
		string s;
		s.append ("<"); s.append (SimdType<Elt>::type_name); s.append ("> ");
		s.append (string (80-(s.size()+36), '.'));

		/* DIF */
		v = in;
		chrono.start();
		for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
			fft_algo.DIF (v.data());
		time = chrono.userElapsedTime()/cnt;
		Miops = 17 * (_k<<(_k-1)) / (1e6 * time); /* 3/2 n log n */
		cout << "  DIF" << s;
		cout.precision(2); cout.width(10); cout<< scientific << time << " s, ";
		cout.precision(2); cout.width(10); cout<<fixed<<Miops << " Miops";
		cout << endl;

		/* DIT */
		v = in;
		chrono.start();
		for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
			fft_algo.DIT (v.data());
		time = chrono.userElapsedTime()/cnt;
		Miops = 17 * (_k<<(_k-1)) / (1e6 * time); /* 3/2 n log n */
		cout << "  DIT" << s;
		cout.precision(2); cout.width(10); cout<< scientific << time << " s, ";
		cout.precision(2); cout.width(10); cout<<fixed<<Miops << " Miops";
		cout << endl;
	}

	/* draw random vector and bench DIF and DIT for all available SIMD implem */
	void bench (unsigned long seed) {
		EltVector in(_n), v(_n);

		/* Generate random input */
		typename Field::RandIter Gen (_F, 0, seed+_k); /*0 has no meaning here*/
		for (auto elt = in.begin(); elt < in.end(); elt++)
			Gen.random (*elt);

		/* NoSimd */
		actual_bench<NoSimd> (in);

		/* Simd128 */
#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
		if (Simd128<Elt>::vect_size == 4 || Simd128<Elt>::vect_size == 8)
			actual_bench<Simd128> (in);
#endif

		/* Simd256 */
#if defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
		if (Simd256<Elt>::vect_size == 4 || Simd256<Elt>::vect_size == 8)
			actual_bench<Simd256> (in);
#endif
	}
};

/* Bench FFT on polynomial with coefficients in Modular<T1,T2> with a random
 * prime p with the required number of bits and such that 2^k divides p-1
 */
template<typename T1, typename T2>
void bench_one_modular_implem (uint64_t bits, size_t k, unsigned long seed)
{
	typedef typename Givaro::Modular<T1, T2> ModImplem;
	typedef typename ModImplem::Residu_t res_t;

	/* First, we check if this Modular implem can handle this many bits */
	bool ok;
	res_t s = 1;
	for (uint64_t i = 0; i+1 < bits; i++, s<<=1); /* at the end s=2^(bits-1) */
	if (s == 0)
		ok = false;
	else {
		s <<= 1; s--; /* now s=2^bits-1 */
		ok = (s <= ModImplem::maxCardinality());
	}
	if (!ok) {
		cout << endl << "# Skipping bench with Modular<" << TypeName<T1>();
		cout << ", " << TypeName<T2>() << ">, bits = " << bits;
		cout << " is too large" << endl;
		return;
	}

	RandomFFTPrime RandomGen (1<<bits, seed);
	T1 p = (T1) RandomGen.randomPrime (k);
	ModImplem GFp(p);

	cout << endl << string (80, '*') << endl;
	cout << "Bench FFT with Modular<" << TypeName<T1>() << ", ";
	cout << TypeName<T2>() << ">, p=" << GFp._p << " (" << bits << " bits, ";
	cout << "n=2^" << k << " divides p-1)" << endl;
	Benchs<ModImplem> C(GFp, k);
	C.bench (seed);
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
	bench_one_modular_implem<double, double> (bits, k, seed);

	/* Bench with Modular<uint16_t,uint32_t> */
	bench_one_modular_implem<uint16_t, uint32_t> (bits, k, seed);

	/* Bench with Modular<uint32_t, uint64_t> */
	bench_one_modular_implem<uint32_t, uint64_t> (bits, k, seed);

#ifdef __FFLASFFPACK_HAVE_INT128
	/* Bench with Modular<uint64_t,uint128_t> */
	bench_one_modular_implem<uint64_t, uint128_t> (bits, k, seed);
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
