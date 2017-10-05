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


template<typename Field, typename Rand>
void randomMatPol(Rand& r,  PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field>& A){
	for(size_t i=0;i<A.size();i++)
		randomMat(r, A[i]);       	
}

template<typename Field, typename Rand>
void randomMatPol(Rand& r,  PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field>& A){
	for(size_t i=0;i<A.rowdim()*A.coldim();i++)
		randomVect(r, A(i));
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


template<typename MatrixP, typename Field, typename RandIter>
bool check_matpol_mul(const Field& fld,  RandIter& Gen, size_t n, size_t d) {       
	MatrixP A(fld,n,n,d),B(fld,n,n,d),C(fld,n,n,2*d-1);

	// Generate random matrix of polynomial
	randomMatPol(Gen,A);
	randomMatPol(Gen,B);
	typedef PolynomialMatrixDomain<Field>    PolMatDom;
	PolMatDom  PMD(fld);
	PMD.mul(C,A,B);
	return check_mul(C,A,B,C.size());
}


template<typename MatrixP, typename Field, typename RandIter>
bool check_matpol_midp(const Field& fld,  RandIter& Gen, size_t n, size_t d) {
	MatrixP A(fld,n,n,d),C(fld,n,n,2*d-1);	
	MatrixP B(fld,n,n,d);
	// Generate random matrix of polynomial
	randomMatPol(Gen,A);
	randomMatPol(Gen,C);
	typedef PolynomialMatrixDomain<Field>    PolMatDom;
	PolMatDom  PMD(fld);
	PMD.midproduct(B,A,C) ;
	return check_midproduct(B,A,C);
}

template<typename MatrixP, typename Field, typename RandIter>
bool check_matpol_midpgen(const Field& fld,  RandIter& Gen, size_t n, size_t d) {
	size_t d0,d1; // n0 must be in [d, 2d-1[

	d1=d/2;
	d0=d-d1;
	
	MatrixP A(fld,n,n,d0+1),B(fld,n,n,d1),C(fld,n,n,d);

	// Generate random matrix of polynomial
	randomMatPol(Gen,A);
	randomMatPol(Gen,C);
	typedef PolynomialMatrixDomain<Field>    PolMatDom;
	PolMatDom  PMD(fld);
	PMD.midproductgen(B,A,C,true,d0+1,d) ;
	return check_midproduct(B,A,C,true,d0+1,d);
}


template<typename MatrixP, typename Field, typename RandIter>
bool debug_midpgen_dlp(const Field& fld,  RandIter& Gen) {
	size_t d0,d1;
	size_t n0,n1;

	n0=22;
	n1=42;
	
		
	MatrixP A(fld,48,48,22),B(fld,48,32,21),C(fld,48,32,42);

	// Generate random matrix of polynomial
	randomMatPol(Gen,A);
	randomMatPol(Gen,C);
	typedef PolynomialMatrixDomain<Field>    PolMatDom;
	PolMatDom  PMD(fld);
	PMD.midproductgen(B,A,C,true,n0,n1) ;
	return check_midproduct(B,A,C,true,n0,n1);
}



template<typename Field>
bool launchTest(const Field& F, size_t n, long b, long d, long seed){
	bool ok=true;
	typename Field::RandIter G(F,b,seed);
	typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;
	ostream& report = LinBox::commentator().report();
	report<<"Polynomial matrix (polfirst) testing over ";F.write(report)<<std::endl;
	ok&=check_matpol_mul<MatrixP> (F,G,n,d);
	ok&=check_matpol_midp<MatrixP> (F,G,n,d);
	ok&=check_matpol_midpgen<MatrixP> (F,G,n,d); 

	//typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> PMatrix;
	// std::cerr<<"Polynomial matrix (matfirst) testing:\n";F.write(std::cerr)<<std::endl;
	// check_matpol_mul<PMatrix> (F,G,n,d);
	// check_matpol_midp<PMatrix> (F,G,n,d);

	report<<"Debugging midpgen for DLP over ";F.write(report)<<std::endl;
	debug_midpgen_dlp<MatrixP>(F,G);
	return ok;
}


bool runTest(uint64_t n, uint64_t d, long seed){

	bool ok=true;

	// fourier prime < 2^(53--log(n))/2
	{
		size_t bits= (53-integer(n).bitsize())/2;
		RandomFFTPrime Rd(1<<bits,seed);
		integer p = Rd.randomPrime(integer(d).bitsize()+1);
		
		Givaro::Modular<double> F((int32_t)p);
		ok&=launchTest (F,n,bits,d,seed);
		
	}
	// normal prime < 2^(53--log(n))/2
	{
		size_t bits= (53-integer(n).bitsize())/2;;
		RandomPrimeIter Rd(bits,seed);
		integer p;
		Rd.random(p);

		Givaro::Modular<double> F((int32_t)p);
		ok&=launchTest (F,n,bits,d,seed);
	}

	// multi-precision prime
	 {
	 	size_t bits=114;
	 	RandomPrimeIter Rd(bits,seed);
	 	integer p= Rd.random();

	 	Givaro::Modular<integer> F1(p);			
	 	ok&=launchTest (F1,n,bits,d,seed);
	 	Givaro::Modular<RecInt::ruint128,RecInt::ruint256> F2(p);
	 	ok&=launchTest (F2,n,bits,d,seed);
	
	 }
	 // over the integer
	{
		Givaro::ZRing<integer> F;
		ok&=launchTest (F,n,128,d,seed);
	 }	
	return ok;
}

int main(int argc, char** argv){
	static size_t  n = 16; // matrix dimension
	static size_t  d = 512; // polynomial size
	static long    seed = time(NULL);
	
	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'd', "-d D", "Set degree of test matrices to D.", TYPE_INT,     &d },
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	return (runTest(n,d,seed)? 0: -1);
} 
 



