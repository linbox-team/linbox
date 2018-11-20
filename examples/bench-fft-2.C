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

#include <linbox/linbox-config.h>

#include <functional>
#include <iostream>
#include <vector>

#include <givaro/modular.h>
#include <givaro/givranditer.h>

using namespace std; 

#include "fflas-ffpack/fflas-ffpack.h"

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-butterflies.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-algorithms.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.h"

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/ring/modular.h"
#include "fflas-ffpack/utils/align-allocator.h"

using namespace LinBox;

uint64_t bit_reverse(uint64_t x){
	uint64_t r = x;
	for(uint64_t i = 0 ; i < 64 ; ++i){
		r <<= 1;
		r |= x & 1;
		x >>= 1;
	}
	return r;
}

template<class Vect>
void sort_vector_bit_reverse(Vect & v){
	std::vector<size_t> v_idx(v.size());
	for(size_t i = 0 ; i < v.size() ; ++i){
		v_idx[i] = i;
	}
	std::sort(v_idx.begin(), v_idx.end(), [](size_t a, size_t b){return bit_reverse(a) < bit_reverse(b);});
	Vect temp = v;
	for(size_t i = 0 ; i < v.size() ; ++i){
		v[i] = temp[v_idx[i]];
	}
}

template <typename Rand, typename Vect>
void randomVect (Rand& r, Vect& v) {
	size_t s = v.size();
	for (size_t i = 0; i < s; ++i)
		r.random(v[i]);
}

template<class Vect>
void print_vectors(const Vect &x, const Vect & y){
	std::ostream_iterator<typename Vect::value_type> out_it(std::cout, ", ");
	std::cout << "x: ";
	std::copy(x.begin(), x.end(), out_it);
	std::cout << std::endl;
	std::cout << "y: ";
	std::copy(y.begin(), y.end(), out_it);
	std::cout << std::endl;
}

template<class Field>
bool check_DFT(const Field &fld, size_t kmax, long seed){
	using Element = typename Field::Element;
	bool passed = true;
	std::cout << "*******************" << std::endl;
	std::cout << "*** Testing DFT ***" << std::endl;
	std::cout << "*******************" << std::endl;
	for(size_t log_convolution_size = 1 ; log_convolution_size < kmax ; ++log_convolution_size){
		size_t convolution_size = 1 << log_convolution_size;

		Element inv_convolution_size;
		fld.inv(inv_convolution_size, convolution_size);

		std::vector<Element, AlignedAllocator<Element, Alignment::DEFAULT>> x(convolution_size), y(convolution_size);

		typename Field::RandIter Gen(fld);
		randomVect(Gen, y);
		x = y;

		FFT_init<Field> fft_init(fld, log_convolution_size);

		FFT_algorithms<Field, NoSimd<Element>> fft_algo_nosimd(fft_init);
		fft_algo_nosimd.DIF_sort(y.data());
		passed &= !std::equal(x.begin(), x.end(), y.begin());
		fft_algo_nosimd.DIT_sort(y.data());
		FFLAS::fscalin(fld, convolution_size, inv_convolution_size, y.data(), 1);

		passed &= std::equal(x.begin(), x.end(), y.begin());

		std::cout << "size: " << log_convolution_size << " DFT sort: " << (passed ? "passed" : "error") << std::endl;
		if(!passed){
			print_vectors(x, y);
			break;
		} 
	}
	return passed;
}

/**************************************
 ****** DFT PERFORMANCE FUNCTION ******
 **************************************/
template<typename Funct, typename FFT, typename Vect>
void DFT_performance(FFT& FFTDom, Funct f, size_t lpts, const Vect& x, string msg){
	Vect z(x);
	auto Functor = bind(f, &FFTDom, &z[0]);
	Timer chrono;
	double time;
	double Miops;
	size_t ct,minct=4;
	ct = 0;
	chrono.start();
	while (chrono.realElapsedTime() < 1){
		for (size_t k=0;k<minct;k++)
			Functor();
		ct+=minct;
	}
	time = chrono.userElapsedTime()/ct;
	Miops = 17 * (lpts<<(lpts-1)) /(1e6 * time); // 3/2 n log n
	msg+="  ";
	msg.resize(45,'.');
	cout << "Timings ... " << msg <<" : ";
	cout.precision(2);
	cout.width(10);
	cout<<scientific<<time << " s, ";
	cout.precision(2);
	cout.width(10);
	cout<<fixed<<Miops << " Miops\n";
}

template<typename Field>
void bench_DIF(const Field& fld, size_t kmax, long seed) { 
	typedef typename Field::Element Element;
	for (size_t lpts = 5; lpts < kmax ; lpts++){
		uint64_t pts = 1UL << lpts;
		cout<<"*********************************************************"<<endl;
		cout<<"*** Benching polynomials of size 2^" << lpts <<endl;
		cout<<"*********************************************************"<<endl;
		vector<Element> x(pts);

		// Generate random inputs
		typename Field::RandIter Gen(fld,seed);
		randomVect(Gen,x); 

		FFT_init<Field> fft_init (fld, lpts);

		FFT_algorithms<Field,NoSimd<typename Field::Element> > fft_algo_nosimd (fft_init);
		using FFT_a = FFT_algorithms<Field,NoSimd<typename Field::Element> >;
		DFT_performance(fft_algo_nosimd,&FFT_a::DIF, lpts, x, "FFT_algorithms<Field,NoSimd>::DIF");
		DFT_performance(fft_algo_nosimd, &FFT_a::DIF_sort, lpts, x, "FFT_algorithms<Field, NoSimd>::DIF_sort");

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
		if (Simd128<typename Field::Element>::vect_size == 4 || Simd128<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd128<typename Field::Element> > fft_algo_simd128 (fft_init);
			using FFT_a128 = FFT_algorithms<Field,Simd128<typename Field::Element> >;
			DFT_performance(fft_algo_simd128,&FFT_a128::DIF, lpts, x, "FFT_algorithms<Field,Simd128>::DIF");   
		}
#endif

#if 0 //defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
		if (Simd256<typename Field::Element>::vect_size == 4 || Simd256<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd256<typename Field::Element> > fft_algo_simd256 (fft_init);
			using FFT_a256 = FFT_algorithms<Field,Simd256<typename Field::Element> >;
			DFT_performance(fft_algo_simd256,&FFT_a256::DIF, lpts, x, "FFT_algorithms<Field,Simd256>::DIF");
		}
#endif
		cout<<"---------------------------------------------------------------"<<endl;

		DFT_performance(fft_algo_nosimd,&FFT_a::DIT, lpts, x, "FFT_algorithms<Field,NoSimd>::DIT");
		DFT_performance(fft_algo_nosimd, &FFT_a::DIT_sort, lpts, x, "FFT_algorithms<Field, NoSimd>::DIT_sort");

#if defined(__FFLASFFPACK_HAVE_SSE4_1_INSTRUCTIONS)
		if (Simd128<typename Field::Element>::vect_size == 4 || Simd128<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd128<typename Field::Element> > fft_algo_simd128 (fft_init);
			using FFT_a128 = FFT_algorithms<Field,Simd128<typename Field::Element> >;
			DFT_performance(fft_algo_simd128,&FFT_a128::DIT, lpts, x, "FFT_algorithms<Field,Simd128>::DIT");
		}
#endif

#if 0 // defined(__FFLASFFPACK_HAVE_AVX2_INSTRUCTIONS)
		if (Simd256<typename Field::Element>::vect_size == 4 || Simd256<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd256<typename Field::Element> > fft_algo_simd256 (fft_init);
			using FFT_a256 = FFT_algorithms<Field,Simd256<typename Field::Element> >;
			DFT_performance(fft_algo_simd256,&FFT_a256::DIT, lpts, x, "FFT_algorithms<Field,Simd256>::DIT");
		}
#endif

		cout<<endl;
	}
}


int main(int argc, char** argv){
	//	if (argc < 2 || argc >3){
	//		cerr<<"usage : prime_bitsize , (seed)"<<endl;
	//		exit(0);
	//	}
	uint64_t bits = 0; //atoi(argv[1]);
	long seed=((argc>2)?atoi(argv[2]):time(NULL));
	size_t l2n = 24;
	size_t k = l2n;
	RandomFFTPrime Rd;
	uint32_t p;


	//Modular<double,double>
	bits = 59;
	Rd = RandomFFTPrime (uint64_t(1)<<bits,seed);
	p = (uint64_t)Rd.randomPrime(l2n);

	cout<<"prime : "<<p<<endl;
	cout<<endl;

	Givaro::Modular<uint64_t,uint128_t> Fi64(p);
	std::cout << "Test Modular<int64_t,uint128_t> : " << std::endl;
	std::cout << ((check_DFT(Fi64, k, seed)) ? "OK" : "KO!!!!") << endl;
	bench_DIF(Fi64,k,seed);

	//Modular<uint32_t,uint64_t>
	bits = 27;
	Rd = RandomFFTPrime (1<<bits,seed);
	p = (uint32_t)Rd.randomPrime(l2n);

	Givaro::Modular<uint32_t,uint64_t> Fi32(p);
	std::cout << "Test Modular<int32_t,uint32_t> : " << std::endl;
	std::cout << ((check_DFT(Fi32, k, seed)) ? "OK" : "KO!!!!") << endl;
	bench_DIF(Fi32, k, seed);

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
