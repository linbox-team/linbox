// examples/cholesky.C
/*
 * from examples/local2.C  -bds, 2016Dec
 * This file is part of the library LinBox.  See COPYING for license info.
 */

/** \file examples/cholesky.C
 * @example  examples/cholesky.C
  \brief mod p^e Smith form by Cholesky style elmination
  \ingroup examples

  \author bds & ml

  */

#include <linbox/linbox-config.h>

#include <iostream>
#include <string>
#include <vector>
#include <list>

using namespace std;

#include <linbox/util/timer.h>

#include <linbox/matrix/dense-matrix.h>
//#include <linbox/matrix/symmetric-embedded.h>
//#include <linbox/matrix/symmetric-vspacked.h>
//#include <linbox/matrix/symmetric-cspacked.h>

#include <linbox/blackbox/fibb-product.h>
#include <linbox/blackbox/transpose-bb.h>
#include <linbox/ring/modular.h>
#include <linbox/matrix/matrix-domain.h>

#include "linbox/algorithms/cholesky.h"

using namespace LinBox;

template <class Ring, class Symmetric>
bool testCholesky(Symmetric & S, const typename Ring::Element & p)
{
  S.write(cout << "S is " << endl,Tag::FileFormat::Plain) << endl;
	const Ring & R = S.field();
	// do it
	UserTimer T;
	typedef typename FIBB<Ring>::MotherMatrix MotherMatrix;
	typedef typename FIBB<Ring>::Matrix Matrix;

	FIBBProduct<Ring> F(R);
	T.start();
	size_t r = cholesky(F, S, p);
	T.stop();
	cout << "time is " << T << endl;
	// check it
  cout << "rank is " << r << endl;
	//DenseMatrix<Ring> 
  	MotherMatrix XX(R, S.coldim(), 1); Matrix X(XX);
	X.random();
  	MotherMatrix YY(R, S.rowdim(), 1); Matrix Y(YY);
  	MotherMatrix ZZ(R, S.rowdim(), 1); Matrix Z(ZZ);
	BlasMatrixDomain<Ring> MD(R);
	MD.mul(Y,S,X), 
	Y.write(cout << "Y is " << endl) << endl;
	F.applyRight(Z,X); // want MD.mul(Z,F,X)

	Z.write(cout << "Z is " << endl) << endl;
	bool pass = MD.areEqual(Y,Z);
  cout << "equality check is " << pass << endl;
  if (not pass and S.rowdim() <= 20) 
  {
	X.write(cout << "X is " << endl) << endl;
	Y.write(cout << "Y = SX is " << endl) << endl;
	Z.write(cout << "Z = FX is " << endl) << endl;
	S.write(cout << "S became " << endl,Tag::FileFormat::Plain) << endl;
	cout << F.bbTag() << " " << F.rowdim() << " " << F.coldim() << endl;
	F.write(cout << "fibb prod is " << endl) << endl;
  }
	return pass;
}
int main(int argc, char* argv[])
{
	if (argc != 2) {

		cout << "usage: " << argv[0] << " n, where "; 

		cout << "Positive n is cube dimension, 2^n is matrix order, \n";
		cout << "n = 0 selects a 5x5 example\n";
		cout << "n < 0 selects a |n|x|n| random dense matrix\n";
		cout << "Currently not checking local SNF\n";

		return 0;
	}

	typedef Givaro::Modular<double> Ring;
	Ring R(7);

	typedef DenseMatrix<Ring> MotherMatrix;
	typedef DenseSubmatrix<Ring> Matrix;

	int n = atoi(argv[1]);
	size_t m;
	bool pass = true;
	if (n > 0)
	{
		m = 1<<n;
		MotherMatrix A_mother(R,m,m);
		Matrix A(A_mother);

		Ring::Element nn; R.init(nn, n);
		for (size_t i = 0; i < m; ++i) 
		{	A.setEntry(i,i,nn);
			for (size_t k = 1; k < m; k<<=1) 
			{	size_t j = (i & k) ? i - k : i + k;
				A.setEntry(i,j,R.mOne);
				//A.setEntry(j,i,R.mOne);
			}
		}
		pass = testCholesky<Ring>(A,R.zero);
	} else if (n == 0) {
		m = 5;
		size_t a = 1;
		MotherMatrix A_mother(R,m,m);
		Matrix A(A_mother);
		A.setEntry(a+0,a+1,A.field().one);
		A.setEntry(a+1,a+0,A.field().one);
		A.setEntry(a+2,a+2,A.field().one);
		MotherMatrix L(R,m,m);
		MotherMatrix LT(R,m,m);
		Ring::RandIter r(R);
		for (size_t i = 0; i < m; ++i) 
		{	L.setEntry(i,i,R.one);
			LT.setEntry(i,i,R.one);
			for (size_t j = 0; j < i; ++j)
			{	typename Ring::Element x; R.init(x);
				r.random(x);
				L.setEntry(i,j,x);
				LT.setEntry(j,i,x);
			}
		}
		MatrixDomain<Ring> FC(R);
		MotherMatrix B(R,m,m);
		FC.mul(B,L,A);
		FC.mul(A,B,LT);

		pass = testCholesky<Ring>(A,R.zero);
	} else {
		m = -n;
		MotherMatrix A_mother(R,m,m);
		Matrix A(A_mother);
		A.random();
		pass = testCholesky<Ring>(A,R.zero);
	}

	
	return pass ? 0 : -1;
} // main

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

