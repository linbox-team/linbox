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

/******************************************************************************/
template<typename Field>
struct Checker {
	typedef typename Field::Element Elt;
	typedef AlignedAllocator<Elt, Alignment::DEFAULT> _Alignement;
	typedef typename std::vector<Elt, _Alignement> EltVector;

	const Field& _F;
	size_t _k;
	size_t _n;
	FFT_transform<Field> _fft;

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

	/* draw random vector and check DIF and DIT for all available SIMD implem */
	bool check (unsigned long seed) {
		EltVector in(_n), in_br(_n), out(_n), out_br(_n), v(_n);
		string s;
		bool passed = true;

		/* Generate random input */
		typename Field::RandIter Gen (_F, 0, seed+_k); /*0 has no meaning here*/
		for (auto elt = in.begin(); elt < in.end(); elt++)
			Gen.random (*elt);

		/* Compute the out vector using a naive polynomial evaluation and set
		 * the bitreversed version of in and out */
		Elt x(1);
		for (size_t i = 0; i < _n; _F.mulin (x, _fft._w), i++) {
			size_t i_br = bitreversed (i);
			in_br[i_br] = in[i];
			horner_eval (out[i], in, x);
			out_br[i_br] = out[i];
		}


		/* DIF : natural order => bitreversed order */
#define DIF_CHECK(method) do {								\
	v = in;													\
	_fft.method (v.data());									\
	for (auto& e: v) /* reduce from < 2p to < p */			\
		if (e > _F._p) e-=_F._p;							\
	bool b = equal (v.begin(), v.end(), out_br.begin());	\
	s = "     "; s.append (#method); s.append (" ");		\
	s.append (string (80-(s.size()+16), '.'));				\
	cout << s << (b ? " success" : " failure") << endl;		\
	passed &= b;											\
	} while (0)

		DIF_CHECK(FFT_DIF_Harvey_mod2p_iterative);
		DIF_CHECK(FFT_DIF_Harvey_mod2p_iterative2x2);
		DIF_CHECK(FFT_DIF_Harvey_mod2p_iterative3x3);
		DIF_CHECK(FFT_DIF_Harvey_mod2p_iterative4x1_SSE);
		DIF_CHECK(FFT_DIF_Harvey_mod2p_iterative4x2_SSE);
#ifdef __LINBOX_HAVE_AVX2_INSTRUCTIONS
		DIF_CHECK(FFT_DIF_Harvey_mod2p_iterative8x1_AVX);
#endif

		/* DIT : bitreversed order => natural order */
#define DIT_CHECK(method) do {								\
	v = in_br;												\
	_fft.method (v.data());									\
	for (auto& e: v) { /* reduce from < 3p to < p */		\
		if (e > _F._p) e-=_F._p;							\
		if (e > _F._p) e-=_F._p;							\
		if (e > _F._p) e-=_F._p;							\
	}														\
	bool b = equal (v.begin(), v.end(), out.begin());		\
	s = "     "; s.append (#method); s.append (" ");		\
	s.append (string (80-(s.size()+16), '.'));				\
	cout << s << (b ? " success" : " failure") << endl;		\
	passed &= b;											\
	} while (0)

		DIT_CHECK(FFT_DIT_Harvey_mod4p_iterative);
		DIT_CHECK(FFT_DIT_Harvey_mod4p_iterative2x2);
		DIT_CHECK(FFT_DIT_Harvey_mod4p_iterative3x3);
		DIT_CHECK(FFT_DIT_Harvey_mod4p_iterative4x1_SSE);
#ifdef __LINBOX_HAVE_AVX2_INSTRUCTIONS
		DIT_CHECK(FFT_DIT_Harvey_mod4p_iterative8x1_AVX);
#endif

		return passed;
	}
};

/* Test FFT on polynomial with coefficients in Modular<T1,T2> with a random
 * prime p with the required number of bits and such that 2^k divides p-1
 */
template<typename T1, typename T2>
bool test_one_modular_implem (uint64_t bits, size_t k, unsigned long seed)
{
	typedef typename Givaro::Modular<T1, T2> ModImplem;

	integer p;
	if (!RandomFFTPrime::randomPrime (p, integer(1)<<bits, k))
		throw LinboxError ("RandomFFTPrime::randomPrime failed");
	ModImplem GFp ((T1) p);
	
	cout << endl << string (80, '*') << endl;
	cout << "Test FFT with Modular<" << TypeName<T1>() << ", ";
	cout << TypeName<T2>() << ">, p=" << GFp._p << " (" << bits << " bits, 2^";
	cout << k << " divides p-1)" << endl;

	bool b = true;
	for (size_t kc = 1; kc <= k; kc++) {
		cout << "** with n=2^" << kc << endl;
		Checker<ModImplem> C(GFp, kc);
		b &= C.check (seed);
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

	cout << "# To rerun this test: test-fft-old -s " << seed << endl;
	cout << "# seed = " << seed << endl;

    RandomFFTPrime::seeding (seed);

	/* Test with Modular<uint32_t, uint64_t>, and 27-bit prime and k=10 */
	pass &= test_one_modular_implem<uint32_t, uint64_t> (27, 10, seed);

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
