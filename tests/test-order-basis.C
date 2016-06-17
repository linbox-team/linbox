/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#include <iostream>
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/randiter/random-fftprime.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/ring/modular.h"

#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"
#include "linbox/algorithms/block-coppersmith-domain.h"

using namespace LinBox;
using namespace std;


//ostream& report = commentator().report();
//ostream& report = std::cout;

template<typename Field, typename Mat>
string check_sigma(const Field& F, const Mat& sigma,  Mat& serie, size_t ord){
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
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
		report<<"error at degree="<<i<<endl;
		T[i].write(report, Tag::FileFormat::Plain);
		report<<"***"<<endl;
		report<<serie<<endl;
		report<<sigma<<endl;	
	}
	
	
	if (i==ord && !nul_sigma)
		msg+="done";
	else
		msg+="error";
	return msg;
}

template<typename MatPol>
bool operator==(const MatPol& A, const MatPol& B){
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	MatrixDomain<typename MatPol::Field> MD(A.field());
	if (A.real_degree()!=B.real_degree()|| A.rowdim()!= B.rowdim() || A.coldim()!=B.coldim()){
		report<<A.size()<<"("<<A.rowdim()<<"x"<<A.coldim()<<") <> "
		    <<B.size()<<"("<<B.rowdim()<<"x"<<B.coldim()<<") <> "<<endl;
		return false;
	}
	size_t i=0;
	while (i<=A.real_degree() && MD.areEqual(A[i],B[i]))
		i++;

	if (i<=A.real_degree() && A.rowdim()<10 && A.coldim()<10){
		report<<"first:"<<endl<<A<<endl;
		report<<"second:"<<endl<<B<<endl;
	}

	return i>A.real_degree();
}
 

template<typename Field, typename RandIter>
void check_sigma(const Field& F, RandIter& Gen, size_t m, size_t n, size_t d) {
	ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	//typedef typename Field::Element Element;
	typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> MatrixP;
	//typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> MatrixP;
	MatrixP Serie(F, m, n,  d);
	MatrixP Sigma1(F, m, m, d+1),Sigma2(F, m, m, d+1),Sigma3(F, m, m, d+1);

	// set the Serie at random
	for (size_t k=0;k<d;++k)
		for (size_t i=0;i<m;++i)
			for (size_t j=0;j<n;++j)
				Gen.random(Serie.ref(i,j,k));

	//report<<"Serie:="<<Serie<<std::endl;
	
	// define the shift
	vector<size_t> shift(m,0);
	vector<size_t> shift2(shift),shift3(shift);

	OrderBasis<Field> SB(F);

	SB.M_Basis(Sigma3, Serie, d, shift3);
	report << "M-Basis       : " <<check_sigma(F,Sigma3,Serie,d)<<endl;
	SB.PM_Basis(Sigma1,Serie, d, shift);
	report << "PM-Basis      : " <<check_sigma(F,Sigma1,Serie,d)<<endl;
	//SB.oPM_Basis(Sigma2, Serie, d, shift2);
	//report << "PM-Basis iter : " <<check_sigma(F,Sigma2,Serie,d)<<endl;

	// if (!(Sigma1==Sigma2)){
	// report<<"---> different basis for PM-Basis and PM-Basis iter"<<endl;
	// report<<Sigma1<<endl;
	// report<<Sigma2<<endl;
	// }
	report<<endl;
}

int main(int argc, char** argv){
	static size_t  m = 64; // matrix dimension
	static size_t  n = 32; // matrix dimension
	static size_t  b = 20; // entries bitsize
	static size_t  d = 32;  // matrix degree
	static long    seed = time(NULL);

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of matrix series to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set column dimension of matrix series to N.", TYPE_INT,     &n },
		{ 'd', "-d D", "Set degree of  matrix series to D.", TYPE_INT,     &d },
		{ 'b', "-b B", "Set bitsize of the matrix entries", TYPE_INT, &b },
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};
	
	parseArguments (argc, argv, args);

	typedef Givaro::Modular<double>              SmallField;	
	typedef Givaro::Modular<Givaro::Integer>      LargeField;

	size_t logd=integer((uint64_t)d).bitsize();
	commentator().start ("Testing order basis computation", "testOrderBasis", 1);

	
	ostream &report = commentator().report (Commentator::LEVEL_ALWAYS, INTERNAL_DESCRIPTION);
	report<<"###  matrix series is of size "<<m<<" x "<<n<<" of degree "<<d<<std::endl;
	if (b < 26){
		if (logd>b-2){
			report<<"degree is to large for field bitsize: "<<b<<std::endl;
			exit(0);
		}
		RandomFFTPrime Rd(1<<b,seed);	
		integer p = Rd.randomPrime(logd+1);
		report<<"# starting sigma basis computation over SmallField [x] with p="<<p<<endl;;		
		SmallField F(p);
		typename SmallField::RandIter G(F,0,seed);
		check_sigma(F,G,m,n,d);
	}
	else {
		RandomPrimeIterator Rd(b,seed);	
		integer p = Rd.randomPrime();
		report<<"# starting sigma basis computation over LargeField Fp[x] with p="<<p<<endl;;		

		LargeField F(p);
		typename LargeField::RandIter G(F,0,seed);
		check_sigma(F,G,m,n,d);
	}

	commentator().stop (MSG_STATUS (true), (const char *) 0, "testOrderBasis"); 
	return 0;
}

 
