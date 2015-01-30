/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

#include <iostream>

#include "linbox/randiter/random-fftprime.h"
#include "linbox/field/modular.h"
#include "linbox/matrix/polynomial-matrix.h"
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
void check_sigma(const Field& F, const RandIter& Gen, size_t n, size_t d) {
	//typedef typename Field::Element Element;
	typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> MatrixP;

	MatrixP Serie(F, 2*n, n,  d);
	MatrixP Sigma1(F, 2*n, 2*n, d+1),Sigma2(F, 2*n, 2*n, d+1),Sigma3(F, 2*n, 2*n, d+1);


	// set the first top n*n part of serie at random and lower part with identity
	for (size_t k=0;k<d;++k)
		for (size_t i=0;i<n;++i)
			for (size_t j=0;j<n;++j)
				Gen.random(Serie.ref(i,j,k));

	for (size_t i=0;i<n;++i)
		Serie.ref(i+n,i,0)=1;

	// define the shift
	vector<size_t> shift(2*n,0);
	for (size_t i=n;i<2*n;++i)
		shift[i]=1;
	vector<size_t> shift2(shift),shift3(shift);

	OrderBasis<Field> SB(F);


	SB.M_Basis(Sigma3, Serie, d, shift3);
	SB.PM_Basis(Sigma1, Serie, d, shift);
	SB.oPM_Basis(Sigma2, Serie, d, shift2);

	std::cout << "M-Basis       : " <<check_sigma(F,Sigma3,Serie,d)<<endl;
	std::cout << "PM-Basis      : " <<check_sigma(F,Sigma1,Serie,d)<<endl;
	std::cout << "PM-Basis iter : " <<check_sigma(F,Sigma2,Serie,d)<<endl;
	if (!(Sigma1==Sigma2))
	cout<<"---> different basis for PM-Basis and PM-Basis iter"<<endl;
	cout<<endl;
}

int main(int argc, char** argv){
	static size_t  n = 32; // matrix dimension
	static size_t  b = 20; // entries bitsize
	static size_t  d = 32;  // matrix degree
	static long    seed = time(NULL);

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'd', "-d D", "Set degree of test matrices to D.", TYPE_INT,     &d },
		{ 'b', "-b B", "Set bitsize of the matrix entries", TYPE_INT, &b },
		{ 's', "-s s", "Set the random seed to a specific value", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	typedef Givaro::Modular<int32_t>           Field;
	//typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> MatrixP;

	if (b> 26){
		cerr<<"bitsize is too large ... exiting"<<endl;
		exit(1);
	}

	RandomFFTPrime Rd(b,seed);
	integer p = Rd.randomPrime(integer(d).bitsize()+1);
	Field F(p);
	std::cout<<"# starting sigma basis computation over Fp[x] with p="<<p<<endl;;

	typename Field::RandIter G(F,0,seed);
	check_sigma(F,G,n,d);

	return 0;
}

