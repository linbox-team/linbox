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

void to_ntt(uint64_t* data, const size_t convolution_size,
            const uint64_t* root_pow, const uint64_t* root_pow_shoup,
            const uint64_t p, const uint64_t p2){
    uint64_t k = convolution_size;
    for(uint64_t m = 1 ; m < convolution_size ; m <<= 1){
        k >>= 1;
        for(uint64_t i = 0 ; i < m ; ++i){
            uint64_t j1 = i * k * 2;
            uint64_t j2 = j1 + k;
            uint64_t w = root_pow[m+i];
            uint64_t w_shoup = root_pow_shoup[m+i];
            for(uint64_t j = j1 ; j < j2 ; ++j){
                uint64_t u = data[j];
                u -= (u >= p2) ? p2 : 0;
                uint64_t v = data[j+k];
                uint64_t temp = static_cast<uint64_t>((static_cast<__uint128_t>(w_shoup)*(v)) >> (8*sizeof(uint64_t)));
                temp = w*v - temp*p;
                data[j] = u+temp;
                data[j+k] = u + (p2 - temp);
            }
        }
    }
    for(uint64_t i = 0 ; i < convolution_size ; ++i){
        data[i] -= (data[i] >= p2) ? p2 : 0;
        data[i] -= (data[i] >= p) ? p : 0;
    }
}


void to_ntt_unroll(uint64_t* data, const size_t convolution_size,
            const uint64_t* root_pow, const uint64_t* root_pow_shoup,
            const uint64_t p, const uint64_t p2){
    uint64_t k = convolution_size;
    for(uint64_t m = 1 ; m < convolution_size/4 ; m <<= 1){
        k >>= 1;
        for(uint64_t i = 0 ; i < m ; ++i){
            uint64_t j1 = i * k * 2;
            uint64_t j2 = j1 + k;
            uint64_t w = root_pow[m+i];
            uint64_t w_shoup = root_pow_shoup[m+i];
            for(uint64_t j = j1 ; j < j2 ; j+=4){
                uint64_t u1 = data[j];
                uint64_t u2 = data[j+1];
                uint64_t u3 = data[j+2];
                uint64_t u4 = data[j+3];
                u1 -= (u1 >= p2) ? p2 : 0;
                u2 -= (u2 >= p2) ? p2 : 0;
                u3 -= (u3 >= p2) ? p2 : 0;
                u4 -= (u4 >= p2) ? p2 : 0;
                uint64_t v1 = data[j+k];
                uint64_t v2 = data[j+k+1];
                uint64_t v3 = data[j+k+2];
                uint64_t v4 = data[j+k+3];
                uint64_t temp1 = static_cast<uint64_t>((static_cast<__uint128_t>(w_shoup)*(v1)) >> (8*sizeof(uint64_t)));
                temp1 = w*v1 - temp1*p;
                data[j] = u1+temp1;
                data[j+k] = u1 + (p2 - temp1);

                uint64_t temp2 = static_cast<uint64_t>((static_cast<__uint128_t>(w_shoup)*(v2)) >> (8*sizeof(uint64_t)));
                temp2 = w*v2 - temp2*p;
                data[j+1] = u2+temp2;
                data[j+k+1] = u2 + (p2 - temp2);

                uint64_t temp3 = static_cast<uint64_t>((static_cast<__uint128_t>(w_shoup)*(v3)) >> (8*sizeof(uint64_t)));
                temp3 = w*v3 - temp3*p;
                data[j+2] = u3+temp3;
                data[j+k+2] = u3 + (p2 - temp3);

                uint64_t temp4 = static_cast<uint64_t>((static_cast<__uint128_t>(w_shoup)*(v4)) >> (8*sizeof(uint64_t)));
                temp4 = w*v4 - temp4*p;
                data[j+3] = u4+temp4;
                data[j+k+3] = u4 + (p - temp4);
            }
        }
    }
    // k = 2, m = conv_size / 4
    for(uint64_t i = 0 ; i < convolution_size/4 ; ++i){
        auto w = root_pow[convolution_size/4+i];
        auto w_shoup = root_pow_shoup[convolution_size/4+i];
        uint64_t u1 = data[i*4];
        uint64_t u2 = data[i*4+1];
        u1 -= (u1 >= p2) ? p2 : 0;
        u2 -= (u2 >= p2) ? p2 : 0;
        uint64_t v1 = data[i*4+2];
        uint64_t v2 = data[i*4+3];
        uint64_t temp1 = static_cast<uint64_t>((static_cast<__uint128_t>(w_shoup)*(v1)) >> (8*sizeof(uint64_t)));
        temp1 = w*v1 - temp1*p;
        data[i*4] = u1+temp1;
        data[i*4+k] = u1 + (p2 - temp1);

        uint64_t temp2 = static_cast<uint64_t>((static_cast<__uint128_t>(w_shoup)*(v2)) >> (8*sizeof(uint64_t)));
        temp2 = w*v2 - temp2*p;
        data[i*4+2] = u2+temp2;
        data[i*4+3] = u2 + (p2 - temp2);
    }
    // k = 1, m = conv_size / 2
    for(uint64_t i = 0 ; i < convolution_size / 2 ; ++i){
        auto w = root_pow[convolution_size/2+i];
        auto w_shoup = root_pow_shoup[convolution_size/2+i];
        uint64_t u1 = data[i*2];
        u1 -= (u1 >= p2) ? p2 : 0;
        uint64_t v1 = data[i*2+k];
        uint64_t temp1 = static_cast<uint64_t>((static_cast<__uint128_t>(w_shoup)*(v1)) >> (8*sizeof(uint64_t)));
        temp1 = w*v1 - temp1*p;
        data[i*2] = u1+temp1;
        data[i*2] -= (data[i*2] >= p2) ? p2 : 0;
        data[i*2] -= (data[i*2] >= p) ? p : 0;
        data[i*2+1] = u1 + (p2 - temp1);
        data[i*2+1] -= (data[i*2+1] >= p2) ? p2 : 0;
        data[i*2+1] -= (data[i*2+1] >= p) ? p : 0;
    }
}
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
	void actual_bench (const EltVector& in) {
		Timer chrono;
		double time, Miops;
		size_t cnt;
		EltVector v(_n);
		string s (30, '.');

		/* to_ntt */
		v = in;
		chrono.start();
		for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
			to_ntt (v.data(), _n, _fft.pow_w_br.data(), _fft.pow_wp_br.data(),
                                                           _fft._pl, _fft._dpl);
		time = chrono.userElapsedTime()/cnt;
		Miops = 17 * (_k<<(_k-1)) / (1e6 * time); /* 3/2 n log n */
		cout << "  to_ntt ......." << s;
		cout.precision(2); cout.width(10); cout<< scientific << time << " s, ";
		cout.precision(2); cout.width(10); cout<<fixed<<Miops << " Miops";
		cout << endl;

		/* to_ntt_unroll */
		v = in;
		chrono.start();
		for (cnt = 0; cnt < min_run || chrono.realElapsedTime() < 1 ; cnt++)
			to_ntt_unroll (v.data(), _n, _fft.pow_w_br.data(),
                                    _fft.pow_wp_br.data(), _fft._pl, _fft._dpl);
		time = chrono.userElapsedTime()/cnt;
		Miops = 17 * (_k<<(_k-1)) / (1e6 * time); /* 3/2 n log n */
		cout << "  to_ntt_unroll " << s;
		cout.precision(2); cout.width(10); cout<< scientific << time << " s, ";
		cout.precision(2); cout.width(10); cout<<fixed<<Miops << " Miops";
		cout << endl;
	}

	/* draw random vector and bench DIF and DIT for all available SIMD implem */
	void bench (unsigned long seed) {
		EltVector in(_n);

		/* Generate random input */
		typename Field::RandIter Gen (_F, 0, seed+_k); /*0 has no meaning here*/
		for (auto elt = in.begin(); elt < in.end(); elt++)
			Gen.random (*elt);

		/* NoSimd */
		actual_bench (in);
	}
};

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

	/* First, we check if this Modular implem can handle this many bits */
	bool ok;
    typedef Modular<uint64_t, uint128_t> Mod64;
	typename Mod64::Residu_t s = 1;
	for (uint64_t i = 0; i+1 < bits; i++, s<<=1); /* at the end s=2^(bits-1) */
	if (s == 0)
		ok = false;
	else {
		s <<= 1; s--; /* now s=2^bits-1 */
		ok = (s <= Mod64::maxCardinality());
	}
    if (!ok) {
        cout << endl << "# Skipping bench with Modular<uint64_t, uint128_t>, "
                     << "bits = " << bits << " is too large" << endl;
        return 0;
    }

    integer p;
    RandomFFTPrime::seeding (seed);
	if (!RandomFFTPrime::randomPrime (p, integer(1)<<bits, k))
		throw LinboxError ("RandomFFTPrime::randomPrime failed");
    Mod64 GFp ((uint64_t) p);

    cout << endl << string (80, '*') << endl;
    cout << "Bench FFT with Modular<uint64_t, uint128_t>, ";
    cout << "p=" << GFp.cardinality() << " (" << bits << " bits, ";
    cout << "n=2^" << k << " divides p-1)" << endl;
    Benchs<Mod64> Bench(GFp, k);
    Bench.bench (seed);

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
