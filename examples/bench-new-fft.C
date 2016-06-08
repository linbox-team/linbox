/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-butterflies.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-fft-algorithms.h"

#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/ring/modular.h"
#include "fflas-ffpack/utils/align-allocator.h"

using namespace LinBox;



template <typename Rand, typename Vect>
void randomVect (Rand& r, Vect& v) {
	size_t s = v.size();
	for (size_t i = 0; i < s; ++i)
		r.random(v[i]);
}


/**********************************
 ****** DFT CHECKING FUNCTION *****
 *********************************/
template<typename T>
struct congruent{
	T p;
	congruent(T _p): p(_p){}
	bool operator()(T a, T b) const { return ((uint64_t)a%(uint64_t)p) == ((uint64_t)b%(uint64_t)p);}
};
template<typename Funct, typename FFT, typename Vect>
bool DFT_sanity_check(FFT& FFTDom, Funct f, const Vect& x, const Vect& y, string msg){
	typedef typename FFT::Element Element ;
	Vect z(x);
	auto Functor = bind(f, &FFTDom, &z[0]);
	Functor();
	msg+="  ";
	msg.resize(45,'.');
	cout<<"  Checking ... "<<msg
	   << (equal(y.begin(),y.end(),z.begin(),congruent<Element>(FFTDom._p))?" done":" error")<<endl;

	if (!(equal(y.begin(),y.end(),z.begin(),congruent<Element>(FFTDom._p)))){
//		std::ostream_iterator<Element> out_it (std::cout,", ");
//		std::copy ( z.begin(), z.end(), out_it );
//		std::cout<<std::endl;
//		std::copy ( y.begin(), y.end(), out_it );
//		std::cout<<std::endl;
		return false;
	}
	return true;
}

template<typename Field>
bool check_DIF(const Field& fld, size_t kmax, long seed) {
	typedef typename Field::Element Element;
	bool passed = true;
	for (size_t lpts = 1; lpts < kmax ; lpts++){
		size_t pts = 1 << lpts;
		cout<<"********************************************************"<<endl;
		cout<<"*** Testing polynomials of size 2^" << lpts <<endl;
		cout<<"********************************************************"<<endl;
		//vector<Element> x(pts),y(pts);
		std::vector<Element,AlignedAllocator<Element, Alignment::DEFAULT>> x(pts),y(pts);

		// Generate random inputs
		typename Field::RandIter Gen(fld);//,fld.characteristic(),seed);
		randomVect(Gen,y);
		x=y;
		
//		FFT_transform<Field> MulDom(fld,lpts);
//		typedef FFT_transform<Field> FFT_t;

		FFT_init<Field> fft_init (fld, lpts);

		FFT_algorithms<Field,NoSimd<typename Field::Element> > fft_algo_nosimd (fft_init);
//		using FFT_a = FFT_algorithms<Field,NoSimd<typename Field::Element> >;


		/* CHECK DIF */
		// compute the correct result
		fft_algo_nosimd.DIF(y.data());

#if defined(__FFLASFFPACK_USE_SIMD)
		// check FFT_algorithms::DIF
		if (Simd128<typename Field::Element>::vect_size == 4 || Simd128<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd128<typename Field::Element> > fft_algo_simd128 (fft_init);
			using FFT_a128 = FFT_algorithms<Field,Simd128<typename Field::Element> >;
			passed &= DFT_sanity_check(fft_algo_simd128,&FFT_a128::DIF,x,y, "FFT_algorithms<Field,Simd128>::DIF");
		}
#endif

#if defined(__FFLASFFPACK_USE_AVX2)
		// check FFT_algorithms::DIF
		if (Simd256<typename Field::Element>::vect_size == 4 || Simd256<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd256<typename Field::Element> > fft_algo_simd256 (fft_init);
			using FFT_a256 = FFT_algorithms<Field,Simd256<typename Field::Element> >;
			passed &= DFT_sanity_check(fft_algo_simd256,&FFT_a256::DIF,x,y, "FFT_algorithms<Field,Simd256>::DIF");
		}
#endif
		cout<<"---------------------------------------------------------------"<<endl;

		/* CHECK DIT */
		// compute the correct result
		y=x;
		fft_algo_nosimd.DIT(y.data());

#if defined(__FFLASFFPACK_USE_SIMD)
		// check FFT_algorithms::DIT
		if (Simd128<typename Field::Element>::vect_size == 4 || Simd128<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd128<typename Field::Element> > fft_algo_simd128 (fft_init);
			using FFT_a128 = FFT_algorithms<Field,Simd128<typename Field::Element> >;
			passed &= DFT_sanity_check(fft_algo_simd128,&FFT_a128::DIT,x,y, "FFT_algorithms<Field,Simd128>::DIT");
		}
#endif

#if defined(__FFLASFFPACK_USE_AVX2)
		// check FFT_algorithms::DIT
		if (Simd256<typename Field::Element>::vect_size == 4 || Simd256<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd256<typename Field::Element> > fft_algo_simd256 (fft_init);
			using FFT_a256 = FFT_algorithms<Field,Simd256<typename Field::Element> >;
			passed &= DFT_sanity_check(fft_algo_simd256,&FFT_a256::DIT,x,y, "FFT_algorithms<Field,Simd256>::DIT");
		}
#endif

		cout<<endl;
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

#if defined(__FFLASFFPACK_USE_SIMD)
		if (Simd128<typename Field::Element>::vect_size == 4 || Simd128<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd128<typename Field::Element> > fft_algo_simd128 (fft_init);
			using FFT_a128 = FFT_algorithms<Field,Simd128<typename Field::Element> >;
			DFT_performance(fft_algo_simd128,&FFT_a128::DIF, lpts, x, "FFT_algorithms<Field,Simd128>::DIF");
		}
#endif

#if defined(__FFLASFFPACK_USE_AVX2)
		if (Simd256<typename Field::Element>::vect_size == 4 || Simd256<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd256<typename Field::Element> > fft_algo_simd256 (fft_init);
			using FFT_a256 = FFT_algorithms<Field,Simd256<typename Field::Element> >;
			DFT_performance(fft_algo_simd256,&FFT_a256::DIF, lpts, x, "FFT_algorithms<Field,Simd256>::DIF");
		}
#endif
		cout<<"---------------------------------------------------------------"<<endl;

		DFT_performance(fft_algo_nosimd,&FFT_a::DIT, lpts, x, "FFT_algorithms<Field,NoSimd>::DIT");

#if defined(__FFLASFFPACK_USE_SIMD)
		if (Simd128<typename Field::Element>::vect_size == 4 || Simd128<typename Field::Element>::vect_size == 8){
			FFT_algorithms<Field,Simd128<typename Field::Element> > fft_algo_simd128 (fft_init);
			using FFT_a128 = FFT_algorithms<Field,Simd128<typename Field::Element> >;
			DFT_performance(fft_algo_simd128,&FFT_a128::DIT, lpts, x, "FFT_algorithms<Field,Simd128>::DIT");
		}
#endif

#if defined(__FFLASFFPACK_USE_AVX2)
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
	size_t l2n = 12;
	size_t k = l2n;
	RandomFFTPrime Rd;
	uint32_t p;

	//Modular<double,double>
	bits = 22;
	Rd = RandomFFTPrime (1<<bits,seed);
	p = (double)Rd.randomPrime(l2n);

	cout<<"prime : "<<p<<endl;
	cout<<endl;

	Givaro::Modular<double,double> Fd(p);
//	cout << "Test Modular<double,double>: " << ((check_DIF(Fd,k,seed))?"OK":"KO!!!!") << endl;

	//Modular<int64_t,uint128_t>
	bits = 59;
	Rd = RandomFFTPrime (1ul<<bits,seed);
	p = (uint64_t)Rd.randomPrime(l2n);

	cout<<"prime : "<<p<<endl;
	cout<<endl;

//	Givaro::Modular<uint64_t,uint128_t> Fi64(p);
//	cout << "Test Modular<int64_t,uint128_t> : " << ((check_DIF(Fi64,k,seed))?"OK":"KO!!!!") << endl;

	//Modular<uint32_t,uint64_t>
	bits = 28;
	Rd = RandomFFTPrime (1<<bits,seed);
	p = (uint32_t)Rd.randomPrime(l2n);

	cout<<"prime : "<<p<<endl;
	cout<<endl;

	Givaro::Modular<uint32_t,uint64_t> Fi32(p);
	cout << "Test Modular<uint32_t,uint64_t>: " << ((check_DIF(Fi32,k,seed))?"OK":"KO!!!!") << endl;

	bench_DIF(Fi32,k,seed);


	//Modular<uint16_t,uint32_t>
	bits = 12;
	k = l2n = 8;
	Rd = RandomFFTPrime (1<<bits,seed);
	p = (uint16_t)Rd.randomPrime(l2n);

	cout<<"prime : "<<p<<endl;
	cout<<endl;

	Givaro::Modular<uint16_t,uint32_t> Fi16(p);
//	cout << "Test Modular<uint16_t,uint32_t> : " << ((check_DIF(Fi16,k,seed))?"OK":"KO!!!!") << endl;


	// Bench FFT

	//	cout << "Test : " << ((check_DIF(Fi16,k,seed))?"OK":"KO!!!!") << endl;
	//	cout << "Test : " << ((check_DIF(Fd,k,seed))?"OK":"KO!!!!") << endl;
	//	bench_DIF(Fi,k,seed);
	//	bench_DIF(Fd,k,seed);


	return 0;
}


