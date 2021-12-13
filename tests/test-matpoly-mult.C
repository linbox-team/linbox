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
    A.random(Gen);
    B.random(Gen);
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
    A.random(Gen);
    C.random(Gen);

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
    A.random(Gen);
    C.random(Gen);


	typedef PolynomialMatrixDomain<Field>    PolMatDom;
	PolMatDom  PMD(fld);
	PMD.midproductgen(B,A,C,true,d0+1,d) ;
	return check_midproduct(B,A,C,true,d0+1,d);
}


template<typename MatrixP, typename Field, typename RandIter>
bool debug_midpgen_dlp(const Field& fld,  RandIter& Gen) {
	size_t n0,n1;

	n0=22;
	n1=42;


	MatrixP A(fld,48,48,22),B(fld,48,32,21),C(fld,48,32,42);

	// Generate random matrix of polynomial
    A.random(Gen);
    C.random(Gen);

	typedef PolynomialMatrixDomain<Field>    PolMatDom;
	PolMatDom  PMD(fld);
	PMD.midproductgen(B,A,C,true,n0,n1) ;
	return check_midproduct(B,A,C,true,n0,n1);
}



template<typename Field>
bool launchTest(const Field& F, size_t n, uint64_t b, long d, long seed){
    bool ok=true;
    typename Field::RandIter::Residu_t samplesize(1); samplesize<<=b;
    typename Field::RandIter G(F,seed,samplesize);
	typedef PolynomialMatrix<Field,PMType::polfirst> MatrixP;
	ostream& report = LinBox::commentator().report();
	report<<"Polynomial matrix (polfirst) testing over ";F.write(report)<<std::endl;
	ok&=check_matpol_mul<MatrixP> (F,G,n,d);
	ok&=check_matpol_midp<MatrixP> (F,G,n,d);
	ok&=check_matpol_midpgen<MatrixP> (F,G,n,d);

	//typedef PolynomialMatrix<Field,PMType::matfirst> PMatrix;
	// std::cerr<<"Polynomial matrix (matfirst) testing:\n";F.write(std::cerr)<<std::endl;
	// check_matpol_mul<PMatrix> (F,G,n,d);
	// check_matpol_midp<PMatrix> (F,G,n,d);

	report<<"Debugging midpgen for DLP over ";F.write(report)<<std::endl;
	debug_midpgen_dlp<MatrixP>(F,G);
	return ok;
}


bool runTest(uint64_t n, uint64_t d, long seed){

	bool ok=true;
	size_t bits= (53-integer(n).bitsize())/2;
	ostream &report = commentator().report (Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION);
	// fourier prime < 2^(53--log(n))/2
	{
        commentator().start("Half wordsize Fourrier prime");
		integer p;
		RandomFFTPrime::seeding (seed);
		if (!RandomFFTPrime::randomPrime (p, 1<<bits, integer(d).bitsize()+1))
			throw LinboxError ("RandomFFTPrime::randomPrime failed");

		Givaro::Modular<double> F((int32_t)p);
		ok&=launchTest (F,n,bits,d,seed);
        commentator().stop(MSG_STATUS (ok), (const char *) 0,"Half wordsize Fourrier prime");

	}
	// normal prime < 2^(53--log(n))/2
	{
        commentator().start("Half wordsize normal prime");
		typedef Givaro::Modular<double> Field;
		PrimeIterator<IteratorCategories::HeuristicTag> Rd(FieldTraits<Field>::bestBitSize(n),seed);
		integer p;
		p=*Rd;
        report<<"prime bits : "<<p.bitsize()<<std::endl;
		Field F((int32_t)p);
		ok&=launchTest (F,n,bits,d,seed);
        commentator().stop(MSG_STATUS (ok), (const char *) 0,"Half wordsize normal prime");
	}

	// multi-precision prime
	 {

         commentator().start("Multiprecision generic prime");
	 	size_t bits=114;
	 	PrimeIterator<IteratorCategories::HeuristicTag> Rd(bits,seed);
	 	integer p= *Rd;
        report<<"prime bits : "<<p.bitsize()<<std::endl;
	 	Givaro::Modular<integer> F1(p);			
	 	ok&=launchTest (F1,n,bits,d,seed);
	 	Givaro::Modular<RecInt::ruint128,RecInt::ruint256> F2(p);
	 	ok&=launchTest (F2,n,bits,d,seed);
	    commentator().stop(MSG_STATUS (ok), (const char *) 0,"Multiprecision generic prime");
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

	commentator().start ("Testing polynomial matrix multiplication", "testMatpolyMult", 1);
    bool pass=    runTest(n,d,seed);
    commentator().stop(MSG_STATUS(pass),(const char *) 0,"testMatpolyMult");
    
    return (pass? 0: -1);
} 

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
