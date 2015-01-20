/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
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

#include <functional>
#include <iostream>
#include <vector>

#include <givaro/modular.h>
#include <givaro/givranditer.h>

using namespace std; 


#include "linbox/algorithms/polynomial-matrix/polynomial-fft-transform.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/field/modular.h"

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
	bool operator()(T a, T b) const { return (a%p) == (b%p);}
};
template<typename Funct, typename FFT, typename Vect>
void DFT_sanity_check(FFT& FFTDom, Funct f, const Vect& x, const Vect& y, string msg){
	typedef typename FFT::Element Element ;
	Vect z(x);
	auto Functor = bind(f, &FFTDom, ref(z));
	Functor();
	msg+="  ";
	msg.resize(45,'.');
	cout<<"  Checking ... "<<msg
	    << (equal(y.begin(),y.end(),z.begin(),congruent<Element>(FFTDom._pl))?" done":" error")<<endl;

	// if (!(equal(y.begin(),y.end(),z.begin(),congruent<Element>(FFTDom._pl)))){
	// 	std::ostream_iterator<Element> out_it (std::cout,", ");
	// 	std::copy ( z.begin(), z.end(), out_it );
	// 	std::cout<<std::endl;
	// 	std::copy ( y.begin(), y.end(), out_it );
	// 	std::cout<<std::endl; 
	// }

}

template<typename Field>
void check_DIF(const Field& fld, size_t kmax, long seed) {  
	typedef typename Field::Element Element;
	for (size_t lpts = 4; lpts < kmax ; lpts++){
		size_t pts = 1 << lpts;
		cout<<"********************************************************"<<endl;
		cout<<"*** Testing polynomials of size 2^" << lpts <<endl;
		cout<<"********************************************************"<<endl;
		vector<Element> x(pts),y(pts);

		// Generate random inputs
		typename Field::RandIter Gen(fld);//,fld.characteristic(),seed);
		randomVect(Gen,y);
		x=y;
		
		FFT_transform<Field> MulDom(fld,lpts);
		typedef FFT_transform<Field> FFT_t;
		typedef vector<Element>  Vect;

		/* CHECK DIF */
		// compute the correct result
		MulDom.FFT_DIF_Harvey_mod2p_iterative(y); 
		// check 2x2 		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative2x2<Vect>,x,y, "DIF_Harvey_mod2p_iterative2x2");
		// check 3x3 		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative3x3<Vect>,x,y, "DIF_Harvey_mod2p_iterative3x3");
		// check 4x1 SSE		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative4x1_SSE<Vect>,x,y, "DIF_Harvey_mod2p_iterative4x1_SSE");
		// check 4x2 SSE		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative4x2_SSE<Vect>,x,y, "DIF_Harvey_mod2p_iterative4x2_SSE");
#ifdef __AVX2__
		// check 8x1 AVX		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative8x1_AVX<Vect>,x,y, "DIF_Harvey_mod2p_iterative8x1_AVX");
#endif
		// check Harvey SSE		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIF_Harvey_SSE<Vect>,x,y, "DIF_Harvey_SSE");
		cout<<"---------------------------------------------------------------"<<endl;
		/* CHECK DIT */
		// compute the correct result
		y=x;
		MulDom.FFT_DIT_Harvey_mod4p_iterative2x2(y);
		// check 2x2 		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIT_Harvey_mod4p_iterative2x2<Vect>,x,y, "DIT_Harvey_mod4p_iterative2x2");
		// check 3x3 		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIT_Harvey_mod4p_iterative3x3<Vect>,x,y, "DIT_Harvey_mod4p_iterative3x3");
		// check 4x1 SSE		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIT_Harvey_mod4p_iterative4x1_SSE<Vect>,x,y, "DIT_Harvey_mod4p_iterative4x1_SSE");
#ifdef __AVX2__
		// check 8x1 AVX		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIT_Harvey_mod4p_iterative8x1_AVX<Vect>,x,y, "DIT_Harvey_mod4p_iterative8x1_AVX");
#endif
		// check Harvey SSE		
		DFT_sanity_check(MulDom,&FFT_t::template FFT_DIT_Harvey_SSE<Vect>,x,y, "DIT_Harvey_SSE");
		
		cout<<endl;
	}
}

/**************************************
 ****** DFT PERFORMANCE FUNCTION ******
 **************************************/
template<typename Funct, typename FFT, typename Vect>
void DFT_performance(FFT& FFTDom, Funct f, size_t lpts, const Vect& x, string msg){
	Vect z(x);
	auto Functor = bind(f, &FFTDom, ref(z));
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
		size_t pts = 1 << lpts;
		cout<<"********************************************************"<<endl;
		cout<<"*** Testing polynomials of size 2^" << lpts <<endl;
		cout<<"********************************************************"<<endl;
		vector<Element> x(pts);

		// Generate random inputs
		typename Field::RandIter Gen(fld,seed);
		randomVect(Gen,x);
		FFT_transform<Field> MulDom(fld,lpts);
		typedef FFT_transform<Field> FFT_t; 
		typedef vector<Element>  Vect;

		// check 2x2 		
		DFT_performance(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative2x2<Vect>,lpts, x,     "DIF_Harvey_mod2p_iterative2x2");
		// check 3x3 		
		DFT_performance(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative3x3<Vect>,lpts, x,     "DIF_Harvey_mod2p_iterative3x3");
		// check 4x1 SSE		
		DFT_performance(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative4x1_SSE<Vect>,lpts, x, "DIF_Harvey_mod2p_iterative4x1_SSE");
		// check 4x2 SSE		
		DFT_performance(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative4x2_SSE<Vect>,lpts, x, "DIF_Harvey_mod2p_iterative4x2_SSE");
#ifdef __AVX2__
		// check 8x1 AVX		
		DFT_performance(MulDom,&FFT_t::template FFT_DIF_Harvey_mod2p_iterative8x1_AVX<Vect>,lpts, x, "DIF_Harvey_mod2p_iterative8x1_AVX");
#endif
		// check Harvey SSE		
		DFT_performance(MulDom,&FFT_t::template FFT_DIF_Harvey_SSE<Vect>,lpts, x, "DIF_Harvey_SSE");
		cout<<"---------------------------------------------------------------"<<endl;
		// check 2x2 		
		DFT_performance(MulDom,&FFT_t::template FFT_DIT_Harvey_mod4p_iterative2x2<Vect>,lpts, x,     "DIT_Harvey_mod4p_iterative2x2");
		// check 3x3 		
		DFT_performance(MulDom,&FFT_t::template FFT_DIT_Harvey_mod4p_iterative3x3<Vect>,lpts, x,     "DIT_Harvey_mod4p_iterative3x3");
		// check 4x1 SSE		
		DFT_performance(MulDom,&FFT_t::template FFT_DIT_Harvey_mod4p_iterative4x1_SSE<Vect>,lpts, x, "DIT_Harvey_mod4p_iterative4x1_SSE");
#ifdef __AVX2__
		// check 8x1 AVX		
		DFT_performance(MulDom,&FFT_t::template FFT_DIT_Harvey_mod4p_iterative8x1_AVX<Vect>,lpts, x, "DIT_Harvey_mod4p_iterative8x1_AVX");
#endif
		// check Harvey SSE		
		DFT_performance(MulDom,&FFT_t::template FFT_DIT_Harvey_SSE<Vect>,lpts, x, "DIT_Harvey_SSE");


		cout<<endl;
	}
}


int main(int argc, char** argv){
	if (argc < 2 || argc >3){
		cerr<<"usage : prime_bitsize , (seed)"<<endl;
		exit(0);
	}
	size_t bits =atoi(argv[1]);
	long seed=((argc>2)?atoi(argv[2]):time(NULL));	
	RandomFFTPrime Rd(bits,seed);
	size_t p = Rd.randomPrime(5);
	size_t k = bits-4;
	cout<<"prime : "<<p<<endl;
	cout<<endl;
	
	Givaro::Modular<int32_t> F(p);
	check_DIF(F,k,seed);
	bench_DIF(F,k,seed);


	return 0;
}
 
 
