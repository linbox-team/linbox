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

#define __FFLASFFPACK_SEQUENTIAL


#include <linbox/linbox-config.h>
#ifdef VERB
#define KARA_TIMING
#define FFT_PROFILE
#endif
#ifdef VERBINT
#define INTFFT_TIMING
#endif


#ifdef HAVE_OPENMP
#include <omp.h>
#define GIVARO_USES_OMP
#include <givaro/givtimer.h>
#define gettime realtime
typedef Givaro::OMPTimer myTimer;
#else
#include <givaro/givtimer.h>
#define gettime usertime
typedef Givaro::Timer myTimer;
#endif


#include <functional>
#include <iostream>
#include <vector>
using namespace std;
#include <linbox/ring/modular.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/randiter/random-fftprime.h>
//#include <linbox/field/unparametric.h>
#include <givaro/zring.h>
#include <recint/rint.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/matrix/polynomial-matrix.h>
#include <linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h>
#include <linbox/algorithms/polynomial-matrix/polynomial-mul.h>



using namespace LinBox;

#ifdef BENCH_NTL
#include "eric-mul.h"
#endif

template <typename Rand, typename Vect>
void randomVect (Rand& r, Vect& v) {
	size_t s = v.size();				   
	for (size_t i = 0; i < s; ++i)
		r.random(v[i]); 
}


template<typename Field, typename RandIter>
void check_pol_mulfft(const Field& fld,  RandIter& Gen, size_t d) {
	typedef typename  Field::Element Element;

	std::vector<Element> A(d),B(d),C(2*d-1),D(2*d-1,0);
	// Generate random matrix of polynomial
	randomVect(Gen,A);
	randomVect(Gen,B);

	std::cout<<"degree= "<<d<<" "<<std::endl;
	typedef PolynomialFFTMulDomain<Field>    FFT;
	FFT  PMFFT(fld);
	myTimer chrono;
	//size_t mmul=2*n*n*n;
	//size_t madd=n*n;
	//size_t fft= 17 *d *log(2.*d)/log(2.);
	//size_t costFFT  = mmul*2*d + 3*madd*fft;
	size_t count=0;
	bool correct=true;
	do {
		std::vector<Element> D(2*d-1,0);
		PMFFT.mul(C,A,B);
		for(size_t i=0;i<d;i++)
			for(size_t j=0;j<d;j++)
				fld.axpyin(D[i+j],A[i],B[j]);		
		size_t k=0;
		while(k<2*d-1 && (D[k]==C[k])) {k++;}

		correct= (correct&&(k==(2*d-1)));
		if (!correct){
			cout<<"A= "; std::copy(A.begin(), A.end(), std::ostream_iterator<Element>(std::cout, " "));cout<<endl;
			cout<<"B= "; std::copy(B.begin(), B.end(), std::ostream_iterator<Element>(std::cout, " "));cout<<endl;
			cout<<"at : "<<k<<endl;
			cout<<"C= "; std::copy(C.begin(), C.end(), std::ostream_iterator<Element>(std::cout, " "));cout<<endl;
			cout<<"D= "; std::copy(D.begin(), D.end(), std::ostream_iterator<Element>(std::cout, " "));cout<<endl;
			
	}

		count++;
	} while (count<1 && correct);
	cout<<"polmul fft is "<<(correct?" OK":" KO")<<endl;
	cout<<endl;
}


template<typename Field, typename RandIter>
void profile_pol_mulfft(const Field& fld,  RandIter& Gen, size_t d,size_t b) {

	std::vector<typename Field::Element> A(d),B(d),C(2*d-1);
	// Generate random matrix of polynomial
	randomVect(Gen,A);
	randomVect(Gen,B);

	std::cout<<"degree= "<<d<<" "<<std::endl;
	typedef PolynomialFFTMulDomain<Field>    FFT;
	FFT  PMFFT(fld);
	myTimer chrono;
	//size_t mmul=2*n*n*n;
	//size_t madd=n*n;
	//size_t fft= 17 *d *log(2.*d)/log(2.);
	//size_t costFFT  = mmul*2*d + 3*madd*fft;
	size_t count=0;
	chrono.start();
	do {
		PMFFT.mul(C,A,B);
		count++;
	} while (false);//count<20 && chrono.userElapsedTime()<1);
	//cout<<"FFT Multiplication total: "<<chrono.userElapsedTime()<<" s, "
	//   <<costFFT/(1e6*chrono.userElapsedTime())<<" Miops"<<endl;
	chrono.stop();
	//cout<<chrono.userElapsedTime()/count<<" (";
	//cout<<chrono.realElapsedTime()/count<<") ";
	//cout<<chrono.gettime()/count<<" ";
	cout<<endl<<endl;
	cout<<"FFLAS code : "<<chrono<<endl;

#ifdef BENCH_NTL
	ZZ p;
	p = RandomBits_ZZ(b);
	p = 2*p + 1;
	ZZ_p::init(p);
	
	ZZ_pX f = random_ZZ_pX(d);
	ZZ_pX g = random_ZZ_pX(d);
	ZZ_pX h;
	
	chrono.clear();
	chrono.start();
	my_mul(h, f, g);
	chrono.stop();
	cout<<"Eric code  : "<<chrono<<endl;
	chrono.clear();
	chrono.start();	
	ZZ_pX h2 = f*g;
	chrono.stop();
	cout<<"NTL code   : "<<chrono<<endl;
#endif


	
}





template<typename Field>
void runTest(const Field& F, long b, long d, long seed, std::string test){
	
	typename Field::RandIter G(F,b,seed);
	//typename Field::RandIter G(F,seed);
	if (test == "all"){
		check_pol_mulfft(F,G,d);
		profile_pol_mulfft(F,G,d,b);
	}
	if (test == "check")
		check_pol_mulfft(F,G,d);
	if (test == "fft")
		profile_pol_mulfft(F,G,d,b);
	if (test == "longfft"){
		size_t D[15]={64,128,256,512,1024,2048,4096,8192};
		for (auto x:D)
			profile_pol_mulfft(F,G,x,b);
	}
}

int main(int argc, char** argv){
	static long    b = 20; // entries bitsize
	static uint64_t d = 32;  // polynomial degree
	static bool    z = false; // computer over  Z[x]
	static long    seed = time(NULL);
	static std::string  test ="all";
	static bool genprime=false;
	
	static Argument args[] = {
		{ 'd', "-d D", "Set degree of test matrices to D.", TYPE_INT,     &d },
		{ 'b', "-b B", "Set bitsize of the matrix entries", TYPE_INT, &b },
		{ 'z', "-z y", "Perform the computation over Z[x]", TYPE_BOOL, &z},
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		{ 't', "-t t", "Choose the targeted test {all,check,fft,longfft}", TYPE_STR, &test},
		{ 'g', "-g y", "Perform the computation with generic prime", TYPE_BOOL, &genprime},
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	if (z){
#ifdef FFT_PROFILER
		FFT_PROF_LEVEL=2;		
#endif
		cout<<"Computation over Z[x]  "<<endl;
		cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
		Givaro::ZRing<integer> F;
		runTest (F,b,d,seed,test);
	}
	else {
		if (b > 29){
#ifdef FFT_PROFILER		
			FFT_PROF_LEVEL=2;
#endif
			// RandomPrimeIter Rd(b,seed);
			// integer p= Rd.random();
			// Givaro::Modular<integer> F(p);
			// //Givaro::Modular<RecInt::ruint128,RecInt::ruint512> F(p);
			// cout<<"Computation over Fp[x] with p=  "<<p<<" (Generic prime)"<<endl;
			// cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
			// runTest (F,n,b,d,seed,test);
		}
		else {
#ifdef FFT_PROFILER 
			FFT_PROF_LEVEL=1; 
#endif
			integer p;
			if (genprime){
				RandomPrimeIter Rd(b,seed);
				p= Rd.random();
			}else {
				RandomFFTPrime Rd(1<<b,seed);
				p = Rd.randomPrime(integer(d).bitsize()+1);
			}
			//Givaro::Modular<int32_t,int64_t> F((int32_t)p);
			Givaro::Modular<double> F((int32_t)p);
			cout<<"Computation over Fp[x] with p=  "<<p<< (genprime?" (Generic prime)":" (FFT prime)")<<endl;
			cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
			runTest (F,b,d,seed,test);
		}
	}
	return 0;
} 
 



