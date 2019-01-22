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

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.h"
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
REGISTER_TYPE_NAME(NoSimd);
REGISTER_TYPE_NAME(Simd128);
REGISTER_TYPE_NAME(Simd256);

/******************************************************************************/
template<typename Field>
struct Benchs {
	typedef typename Field::Element Elt;
	typedef AlignedAllocator<Elt, Alignment::DEFAULT> _Alignement;
	typedef typename std::vector<Elt, _Alignement> EltVector;

	const Field& _F;
	size_t _k;
	size_t _n;
	FFT_transform<Field> _fft;

	/* ctor */
	Benchs (const Field& F, size_t k) : _F(F), _k(k), _n(1<<k), _fft(F,k) {
	}

	/* draw random vector and bench DIF and DIT for all available SIMD implem */
	void bench (unsigned long seed) {
		EltVector in(_n), v(_n);
		size_t min_run = 4; /* do at least this number of run of fft */

		/* Generate random input */
		typename Field::RandIter Gen (_F, 0, seed+_k); /*0 has no meaning here*/
		for (auto elt = in.begin(); elt < in.end(); elt++)
			Gen.random (*elt);

		/* DIF */
#define DIF_BENCH(method) do {												\
	string s; Timer chrono; double time, Miops; size_t cnt;					\
	v = in;																	\
	chrono.start();															\
	for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)	\
		_fft.method (v.data());												\
	time = chrono.userElapsedTime()/cnt;									\
	Miops = 17 * (_k<<(_k-1)) / (1e6 * time); /* 3/2 n log n */				\
	s = "  "; s.append (#method); s.append (" ");							\
	s.append (string (80-(s.size()+32), '.'));								\
	cout << s;																\
	cout.precision(2); cout.width(10); cout<< scientific << time << " s, ";	\
	cout.precision(2); cout.width(10); cout<<fixed<<Miops << " Miops";		\
	cout << endl;															\
	} while (0)

		/*
	for (auto& e: v)
		if (e > _F._p) e-=_F._p;
			*/

		DIF_BENCH(FFT_DIF_Harvey_mod2p_iterative);
		DIF_BENCH(FFT_DIF_Harvey_mod2p_iterative2x2);
		DIF_BENCH(FFT_DIF_Harvey_mod2p_iterative3x3);
		DIF_BENCH(FFT_DIF_Harvey_mod2p_iterative4x1_SSE);
		DIF_BENCH(FFT_DIF_Harvey_mod2p_iterative4x2_SSE);
#ifdef __LINBOX_HAVE_AVX2_INSTRUCTIONS
		DIF_BENCH(FFT_DIF_Harvey_mod2p_iterative8x1_AVX);
#endif
		DIF_BENCH(FFT_DIF_Harvey);

		/* DIT */
#define DIT_BENCH(method) DIF_BENCH(method)

		DIT_BENCH(FFT_DIT_Harvey_mod4p_iterative);
		DIT_BENCH(FFT_DIT_Harvey_mod4p_iterative2x2);
		DIT_BENCH(FFT_DIT_Harvey_mod4p_iterative3x3);
		DIT_BENCH(FFT_DIT_Harvey_mod4p_iterative4x1_SSE);
#ifdef __LINBOX_HAVE_AVX2_INSTRUCTIONS
		DIT_BENCH(FFT_DIT_Harvey_mod4p_iterative8x1_AVX);
#endif
		DIF_BENCH(FFT_DIT_Harvey);
	}
};

/* Bench FFT on polynomial with coefficients in Modular<T1,T2> with a random
 * prime p with the required number of bits and such that 2^k divides p-1
 */
template<typename T1, typename T2>
void bench_one_modular_implem (uint64_t bits, size_t k, unsigned long seed)
{
	typedef typename Givaro::Modular<T1, T2> ModImplem;

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
	FFLAS::writeCommandString (cout, args, "benchmark-fft-old") << endl;

	if (k >= bits) {
		cerr << "Error, k=" << k << " must be smaller than nbits=" << bits;
		cerr << endl;
		return 1;
	}

	if (bits > 27) {
		cerr << "Error, nbits=" << bits << "must be smaller or equal to 27";
		cerr << endl;
		return 1;
	}

	/* Bench with Modular<uint32_t, uint64_t>, */
	bench_one_modular_implem<uint32_t, uint64_t> (bits, k, seed);

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
