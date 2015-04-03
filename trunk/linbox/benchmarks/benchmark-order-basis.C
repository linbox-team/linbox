/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#include <iostream>
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/field/modular.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"
#include "linbox/algorithms/block-coppersmith-domain.h"

using namespace LinBox;
using namespace std;

template<typename Field, typename Mat>
string check_sigma(const Field& F, const Mat& sigma,  Mat& serie, size_t ord){
	Mat T(F,sigma.rowdim(),serie.coldim(),sigma.size()+serie.size()-1);
	PolynomialMatrixMulDomain<Field> PMD(F);
	PMD.mul(T,sigma,serie);
	MatrixDomain<Field> MD(F);
	size_t i=0;
	string msg(".....");
	bool nul_sigma=true;
	while(i<ord && MD.isZero(T[i])){
		if (!MD.isZero(sigma[i])) nul_sigma=false;		
		i++;
	}
	if (i<ord){
		cout<<"error at degree="<<i<<endl;
		T[i].write(std::cout, Tag::FileFormat::Plain);
		cout<<"***"<<endl;
		cout<<serie<<endl;
		cout<<sigma<<endl;	
	}
	
	
	if (i==ord && !nul_sigma)
		msg+="done";
	else
		msg+="error";
	return msg;
}

template<typename MatPol>
bool operator==(const MatPol& A, const MatPol& B){
	MatrixDomain<typename MatPol::Field> MD(A.field());
	if (A.size()!=B.size()|| A.rowdim()!= B.rowdim() || A.coldim()!=B.coldim()){
		cout<<A.size()<<"("<<A.rowdim()<<"x"<<A.coldim()<<") <> "
		    <<B.size()<<"("<<B.rowdim()<<"x"<<B.coldim()<<") <> "<<endl;
		return false;
	}
	size_t i=0;
	while (i<A.size() && MD.areEqual(A[i],B[i]))
		i++;

	if (i<A.size() && A.rowdim()<10 && A.coldim()<10){
		cout<<"first:"<<endl<<A<<endl;
		cout<<"second:"<<endl<<B<<endl;
	}

	return i==A.size();
}
 

template<typename Field, typename RandIter>
void bench_sigma(const Field& F, const RandIter& Gen, size_t m, size_t n, size_t d, string target) {
	//typedef typename Field::Element Element;
	//typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> MatrixP;
	typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;

	MatrixP Serie(F, m, n,  d);
	MatrixP Sigma2(F, m, m, d+1);
	

	
	// set the Serie at random
	for (size_t k=0;k<d;++k)
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<n;++j)
				Gen.random(Serie.ref(i,j,k));

	// define the shift
	vector<size_t> shift1(m,0);
	vector<size_t> shift2(shift1);
	vector<size_t> shift3(shift1);

	OrderBasis<Field> SB(F);
	Timer chrono;
#ifndef __MINMEMORY
	if (target=="ALL"){
		MatrixP Sigma1(F, m, m, d+1);
		chrono.start();
		SB.M_Basis(Sigma1, Serie, d, shift1);
		chrono.stop();
		std::cout << "M-Basis       : " <<chrono.usertime()<<" s"<<std::endl;
	}
#endif
	chrono.clear();		
	chrono.start();
	SB.PM_Basis(Sigma2, Serie, d, shift2);
	chrono.stop();
	std::cout << "PM-Basis      : " <<chrono.usertime()<<" s"<<std::endl;
	chrono.clear();

	// MatrixP Sigma3(F, m, m, d+1);
	// chrono.start();
	// SB.oPM_Basis(Sigma3, Serie, d, shift3);
	// chrono.stop();
	// std::cout << "PM-Basis iter : " <<chrono.usertime()<<" s"<<std::endl;
	std::cout<<endl;
}

int main(int argc, char** argv){
	static size_t  m = 64; // matrix dimension
	static size_t  n = 32; // matrix dimension
	static size_t  b = 20; // entries bitsize
	static size_t  d = 32;  // matrix degree
	static long    seed = time(NULL);
	static string target="BEST";

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of matrix series to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set column dimension of matrix series to N.", TYPE_INT,     &n },
		{ 'd', "-d D", "Set degree of  matrix series to D.", TYPE_INT,     &d },
		{ 'b', "-b B", "Set bitsize of the matrix entries", TYPE_INT, &b },
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		{ 't', "-t T", "Set the targeted benchmark {ALL, BEST}.",            TYPE_STR , &target },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	
	typedef Givaro::Modular<double>              SmallField;	
	typedef Givaro::Modular<Givaro::Integer>      LargeField;

	size_t logd=integer(d).bitsize();
	
	std::cout<<"###  matrix series is of size "<<m<<" x "<<n<<" of degree "<<d<<std::endl;
	if (b < 26){
		if (logd>b-2){
			std::cout<<"degree is to large for field bitsize: "<<b<<std::endl;
			exit(0);
		}
		RandomFFTPrime Rd(b,seed);	
		integer p = Rd.randomPrime(logd+1);
		std::cout<<"# starting sigma basis computation over Fp[x] with p="<<p<<endl;;		
		SmallField F(p);
		typename SmallField::RandIter G(F,0,seed);
		bench_sigma(F,G,m,n,d,target);
	}
	else {
		RandomPrimeIterator Rd(b,seed);	
		integer p = Rd.randomPrime();
		std::cout<<"# starting sigma basis computation over Fp[x] with p="<<p<<endl;;		

		LargeField F(p);
		typename LargeField::RandIter G(F,0,seed);
		bench_sigma(F,G,m,n,d,target);
	}
	
	
	return 0;
}

