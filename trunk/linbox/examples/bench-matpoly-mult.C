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
#include <linbox/field/modular.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/randiter/random-fftprime.h>
#include <linbox/field/unparametric.h>
//#include <givaro/unparametric.h>
#include <linbox/matrix/matrix-domain.h>
#include <linbox/util/commentator.h>
#include <linbox/util/timer.h>
#include <linbox/matrix/polynomial-matrix.h>
#include <linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h>


#ifdef BENCH_FLINT
#define __GMP_BITS_PER_MP_LIMB 64
extern "C" {
#include "flint/longlong.h"
#include "flint/ulong_extras.h"
#include "flint/nmod_poly_mat.h"
#include "flint/flint.h"
}
#endif

using namespace LinBox;



template <typename Rand, typename Vect>
void randomVect (Rand& r, Vect& v) {
	size_t s = v.size();
	for (size_t i = 0; i < s; ++i)
		r.random(v[i]);
}

template <typename Rand, typename Mat>
void randomMat (Rand& r, Mat& m) {
	for (size_t i = 0; i < m.rowdim(); ++i)
		for (size_t j = 0; j < m.coldim(); ++j)
			r.random(m.refEntry(i,j));
}

template<typename MatPol>
bool operator==(const MatPol& A, const MatPol& B){
	MatrixDomain<typename MatPol::Field> MD(A.field());
	if (A.size()!=B.size()|| A.rowdim()!= B.rowdim() || A.coldim()!=B.coldim())
		return false;
	size_t i=0;
	while (i<A.size() && MD.areEqual(A[i],B[i]))
		i++;

	if (i<A.size() && A.rowdim()<10 && A.coldim()<10){
		cout<<"first:"<<endl<<A<<endl;
		cout<<"second:"<<endl<<B<<endl;
	}

	return i==A.size();
}

/******************************************
 ****** MATPOLY MUL CHECKING FUNCTION *****
 ******************************************/
template<typename T>
struct congruent{
	T p;
	congruent(T _p): p(_p){}
	bool operator()(T a, T b) const { return (a%p) == (b%p);}
};

template<typename MULDOM, typename MatPol>
void MATPOLMUL_sanity_check(MULDOM& MulDom, const MatPol& C, const MatPol& A, const MatPol& B, string msg){
	MatPol CC(C.field(),C.rowdim(),C.coldim(),C.size());
	//auto Functor = bind(f, &MulDom, ref(CC),A,B);
	//Functor();
#ifdef FFT_PROFILER
	FFT_PROF_LEVEL=3;
#endif
	MulDom.mul(CC,A,B);
	msg+="  ";
	msg.resize(45,'.');
	cout<<"  Checking ... "<<msg
	    << ((C==CC)?" done":" error")<<endl;
}

template<typename Field, typename RandIter>
void check_matpol_mul(const Field& fld, const RandIter& Gen, size_t n, size_t d) {
	typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> MatrixP;
	MatrixP A(fld,n,n,d),B(fld,n,n,d),C(fld,n,n,2*d-1);

	// Generate random matrix of polynomial
	for (size_t i=0;i<d;i++){
		randomMat(Gen,A[i]);
		randomMat(Gen,B[i]);
	}

	typedef PolynomialMatrixNaiveMulDomain<Field>       Naive;
	typedef PolynomialMatrixKaraDomain<Field>            Kara;
	typedef PolynomialMatrixFFTMulDomain<Field>           FFT;

	Naive NMD(fld);
	Kara PMKD(fld);
	FFT  PMFFT(fld);

	// compute the correct result
	for (size_t r=0;r<n;r++)
		for (size_t c=0;c<n;c++)
			for (size_t k=0;k<n;k++)
				for (size_t i=0;i<A.size();i++)
					for (size_t j=0;j<B.size();j++)
						fld.axpyin(C.ref(r,c,i+j),A.get(r,k,i),B.get(k,c,j));


	// check naive
	MATPOLMUL_sanity_check(NMD,C,A,B, "Naive Multiplication");
	// check karatsuba
	MATPOLMUL_sanity_check(PMKD,C,A,B, "Karatsuba Multiplication");
	// check fft
	MATPOLMUL_sanity_check(PMFFT,C,A,B, "FFT Multiplication");

	cout<<endl;
}


/***********************************************
 ****** MATPOLY MUL  PERFORMANCE FUNCTION ******
 ***********************************************/
template<typename MULDOM, typename MatPol>
void MATPOLMUL_performance(MULDOM& MulDom,  const MatPol& A, const MatPol& B, double Miops, string msg){
	MatPol C(A.field(),A.rowdim(),A.coldim(),A.size()+B.size()-1);
	//auto Functor = bind(f, &MulDom, ref(C),A,B);
	Timer chrono;
	double time;
	size_t ct,minct=4;
	size_t prec=6;
	ct = 0;
	chrono.start();
	while (chrono.realElapsedTime() < 1){
		for (size_t k=0;k<minct;k++)
			//Functor();
			MulDom.mul(C,A,B);
		ct+=minct;
	}

	time = chrono.userElapsedTime()/ct;
	Miops/=(1e6*time);
	msg+="  ";
	msg.resize(45,'.');
	cout << "Timings ... " << msg <<" : ";
	cout.precision(prec);
	cout.width(10);
	cout<<fixed<<time << " s, ";
	cout.precision(2);
	cout.width(10);
	cout<<fixed<<Miops << " Miops\n";
}


template<typename Field, typename RandIter>
void bench_matpol_mul(const Field& fld, const RandIter& Gen, size_t n, size_t d) {
	typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> MatrixP;
	MatrixP A(fld,n,n,d),B(fld,n,n,d),C(fld,n,n,2*d-1);

	// Generate random matrix of polynomial
	for (size_t i=0;i<d;i++){
		randomMat(Gen,A[i]);
		randomMat(Gen,B[i]);
	}

	typedef PolynomialMatrixNaiveMulDomain<Field>       Naive;
	typedef PolynomialMatrixKaraDomain<Field>            Kara;
	typedef PolynomialMatrixFFTMulDomain<Field>           FFT;

	Naive NMD(fld);
	Kara PMKD(fld);
	FFT  PMFFT(fld);
	size_t mmul=2*n*n*n;
	size_t madd=n*n;
	size_t kara=pow((double)d, log(3.)/log(2.));
	size_t fft= 17 *d *log(2.*d)/log(2.);
	size_t costNaive= mmul*d*d + madd*(d-1)*(d-1);
	size_t costKara = mmul*kara+ 6*madd*kara;
	size_t costFFT  = mmul*2*d + 3*madd*fft;

#ifdef FFT_PROFILER
	FFT_PROF_LEVEL=3;
#endif
	// bench naive
	MATPOLMUL_performance(NMD,A,B,costNaive, "Naive Multiplication");
	// bench karatsuba
	MATPOLMUL_performance(PMKD,A,B,costKara, "Karatsuba Multiplication");
	// bench fft
	MATPOLMUL_performance(PMFFT,A,B,costFFT, "FFT Multiplication");



#ifdef BENCH_FLINT
	Timer chrono;
	nmod_poly_mat_t AA,BB,CC;
	nmod_poly_mat_init(AA,n,n,fld.cardinality());
	nmod_poly_mat_init(BB,n,n,fld.cardinality());
	nmod_poly_mat_init(CC,n,n,fld.cardinality());
	flint_rand_t state;
	flint_randinit(state);
	nmod_poly_mat_randtest(AA,state,d);
	nmod_poly_mat_randtest(BB,state,d);
	cout<<"-----------------------"<<endl;
	chrono.start();
 	nmod_poly_mat_mul(CC,AA,BB);
	cout.precision(6);
	chrono.stop();
	cout<<"FLINT Multiplication    : "<<chrono.usertime()<<" s"<<endl;
	chrono.start();
	nmod_poly_mat_mul_KS(CC,AA,BB);
	cout.precision(6);
	chrono.stop();
	cout<<"FLINT Multiplication KS : "<<chrono.usertime()<<" s"<<endl;
	chrono.start();
	nmod_poly_mat_mul_interpolate(CC,AA,BB);
	cout.precision(6);
	chrono.stop();
	cout<<"FLINT Multiplication E/I: "<<chrono.usertime()<<" s"<<endl;
#endif
	cout<<endl;
}



template<typename Field, typename RandIter>
void profile_matpol_mulfft(const Field& fld, const RandIter& Gen, size_t n, size_t d) {
	typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;
	MatrixP A(fld,n,n,d),B(fld,n,n,d),C(fld,n,n,2*d-1);

	// Generate random matrix of polynomial
	for (size_t i=0;i<n*n;i++){
		randomVect(Gen,A(i));
		randomVect(Gen,B(i));
	}
	cout<<n<<" "<<d<<" ";
	typedef PolynomialMatrixFFTMulDomain<Field>    FFT;
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
	cout<<chrono.gettime()/count<<" ";

#ifdef BENCH_FLINT
	nmod_poly_mat_t AA,BB,CC;
	nmod_poly_mat_init(AA,n,n,fld.cardinality());
	nmod_poly_mat_init(BB,n,n,fld.cardinality());
	nmod_poly_mat_init(CC,n,n,fld.cardinality());
	flint_rand_t state;
	flint_randinit(state);
	nmod_poly_mat_randtest(AA,state,d);
	nmod_poly_mat_randtest(BB,state,d);
	count=0;
	chrono.start();
	do {
		nmod_poly_mat_mul(CC,AA,BB);
		count++;
	} while (count<20 && chrono.userElapsedTime()<1);
	cout<<chrono.userElapsedTime()/count<<" ";
	//cout<<"FLINT Multiplication: "<<chrono.userElapsedTime()<<" s, "
	//    <<costFFT/(1e6*chrono.userElapsedTime())<<" Miops"<<endl;

#endif
	cout<<endl;
}



template<typename Field, typename RandIter>
void profile_matpol_mulkara(const Field& fld, const RandIter& Gen, size_t n, size_t d) {
	typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> MatrixP;
	MatrixP A(fld,n,n,d),B(fld,n,n,d),C(fld,n,n,2*d-1);

	// Generate random matrix of polynomial
	for (size_t i=0;i<d;i++){
		randomMat(Gen,A[i]);
		randomMat(Gen,B[i]);
	}

	typedef PolynomialMatrixKaraDomain<Field>    Kara;
	Kara  PM(fld);
	Timer chrono;
	chrono.start();
	PM.mul(C,A,B);
	size_t mmul=2*n*n*n;
	size_t madd=n*n;
	size_t kara=pow((double)d, log(3.)/log(2.));
	size_t costKara = mmul*kara+ 6*madd*kara;
	cout<<"Kara Multiplication total: "<<chrono.userElapsedTime()<<" s, "
	    <<costKara/(1e6*chrono.userElapsedTime())<<" Miops"<<endl;
}


template<typename Field>
void runTest(const Field& F, size_t n, long b, long d, long seed, string test){
	//typename Field::RandIter G(F,b,seed);
	typename Field::RandIter G(F,seed);
	if (test == "check"|| test == "all")
		check_matpol_mul(F,G,n,d);
	if (test == "bench" || test == "all")
		bench_matpol_mul(F,G,n,d);
	if (test == "fft")
		profile_matpol_mulfft(F,G,n,d);
	if (test == "longfft"){
		size_t N[15]={16,16,16,16,16,16,16,16,             64,128,256,512,512,1024,2048};
		size_t D[15]={64,128,256,512,1024,2048,4096,8192,1024,512,256,512,128,  64,32 };
		for (size_t i=0;i<15;i++)
			profile_matpol_mulfft(F,G,N[i],D[i]);
	}
	if (test == "kara")
		profile_matpol_mulkara(F,G,n,d);
}

int main(int argc, char** argv){
	static size_t  n = 32; // matrix dimension
	static long    b = 20; // entries bitsize
	static long    d = 32;  // matrix degree
	static bool    z = false; // computer over  Z[x]
	static long    seed = time(NULL);
	static string  test ="all";

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'd', "-d D", "Set degree of test matrices to D.", TYPE_INT,     &d },
		{ 'b', "-b B", "Set bitsize of the matrix entries", TYPE_INT, &b },
		{ 'z', "-z y", "Perform the computation over Z[x]", TYPE_BOOL, &z},
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		{ 't', "-t t", "Choose the targeted test {all,check,bench,fft,kara,longfft}", TYPE_STR, &test},
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	if (z){
#ifdef FFT_PROFILER
		FFT_PROF_LEVEL=2;
#endif
		cout<<"Computation over Z[x]  "<<endl;
		cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
		Givaro::UnparametricRing<integer> F;
		runTest (F,n,b,d,seed,test);
	}
	else {
		if (b > 29){
#ifdef FFT_PROFILER
			FFT_PROF_LEVEL=2;
#endif
			RandomPrimeIter Rd(b,seed);
			integer p= Rd.random();
			Givaro::Modular<integer> F(p);
			cout<<"Computation over Fp[x] with p=  "<<p<<" (Generic prime)"<<endl;
			cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
			runTest (F,n,b,d,seed,test);
		}
		else {
#ifdef FFT_PROFILER
			FFT_PROF_LEVEL=1;
#endif
			RandomFFTPrime Rd(b,seed);
			integer p = Rd.randomPrime(integer(d).bitsize()+1);
			Givaro::Modular<int32_t> F((int32_t)p);
			cout<<"Computation over Fp[x] with p=  "<<p<<" (FFT prime)"<<endl;
			cout<<"++++++++++++++++++++++++++++++++++++"<<endl;
			runTest (F,n,b,d,seed,test);
		}
	}
	return 0;
} 




