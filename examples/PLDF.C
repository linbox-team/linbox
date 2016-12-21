/*
 * from examples/local2.C  -bds, 2016Dec
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

/** \file examples/smith.C
 * @example  examples/smith.C
  \brief mod m Smith form by elmination
  \ingroup examples

  \author bds & zw

  Various Smith form algorithms may be used for matrices over the
  integers or over Z_m.  Moduli greater than 2^32 are not supported.
  Several types of example matrices may be constructed or matrix read from file.
  Run the program with no arguments for a synopsis of the
  command line parameters.

  For the "adaptive" method, the matrix must be over the integers.
  This is expected to work best for large matrices.

  For the "2local" method, the computation is done mod 2^32.

  For the "local" method, the modulus must be a prime power.

  For the "ilio" method, the modulus may be arbitrary composite.
  If the modulus is a multiple of the integer determinant, the intege Smith form is obtained.  Determinant plus ilio may be best for smaller matrices.

  This example was used during the design process of the adaptive algorithm.
  */

#include <linbox/linbox-config.h>

#include <iostream>
#include <string>
#include <vector>
#include <list>

using namespace std;

#include <linbox/util/timer.h>

#include <linbox/matrix/dense-matrix.h>
#include <linbox/matrix/permutation-matrix.h>
#include <linbox/ring/modular.h>
#include <linbox/matrix/matrix-domain.h>

#include "linbox/algorithms/PLDF.h"

using namespace LinBox;

// #ifndef BIG

template<class I1, class Lp> void distinct (I1 a, I1 b, Lp& c);
template <class I> void display(I b, I e);

int main(int argc, char* argv[])
{
	if (argc < 2) {

		cout << "usage: " << argv[0] << " n, where "; 

		cout << "n is cube dimension, 2^n is matrix order, \n";

		return 0;
	}

	UserTimer T;

	typedef Givaro::Modular<double> Ring;
	Ring R(7);
	int n = atoi(argv[1]);
	size_t m;
	if (n != 0)
	{
	m = 1<<n;
	DenseMatrix<Ring> A_mother(R,m,m);
	DenseSubmatrix<Ring> A(A_mother);

	Ring::Element nn; R.init(nn, n);
	for (size_t i = 0; i < m; ++i) {
		A.setEntry(i,i,nn);
		for (size_t k = 1; k < m; k<<=1) {
			size_t j = (i & k) ? i - k : i + k;
			A.setEntry(i,j,R.mOne);
			//A.setEntry(j,i,R.mOne);
		}
	}

	A.write(cout,Tag::FileFormat::Plain) << endl;
	LinBox::MatrixPermutation<uint32_t> P; P.resize((uint32_t)m);
	//P.write(cout) << endl;
	size_t k = PLD(P,A);
	cout << "k is " << k << endl;
	A.write(cout << "A became " << endl,Tag::FileFormat::Plain) << endl;
	//P.write(cout << "P became " << endl) << endl;
	} else {
		m = 5;
		size_t a = 1;
	DenseMatrix<Ring> A_mother(R,m,m);
	DenseSubmatrix<Ring> A(A_mother);
		A.setEntry(a+0,a+1,A.field().one);
		A.setEntry(a+1,a+0,A.field().one);
		A.setEntry(a+2,a+2,A.field().one);
	DenseMatrix<Ring> L(R,m,m);
	DenseMatrix<Ring> LT(R,m,m);
	Ring::RandIter r(R);
	for (size_t i = 0; i < m; ++i) {
		L.setEntry(i,i,R.one);
		LT.setEntry(i,i,R.one);
		for (size_t j = 0; j < i; ++j){
			typename Ring::Element x; R.init(x);
			r.random(x);
			L.setEntry(i,j,x);
			LT.setEntry(j,i,x);
		}
	}

	MatrixDomain<Ring> FC(R);
	DenseMatrix<Ring> B(R,m,m);
	FC.mul(B,L,A);
	FC.mul(A,B,LT);
		
	A.write(cout,Tag::FileFormat::Plain) << endl;
	LinBox::MatrixPermutation<uint32_t> P; P.resize((uint32_t)m);
	//P.write(cout) << endl;
	size_t k = PLD(P,A);
	cout << "k is " << k << endl;
	A.write(cout << "A became " << endl,Tag::FileFormat::Plain) << endl;
	//P.write(cout << "P became " << endl) << endl;
	}
#if 0	
	typedef list< Ring::Element > List;

	SmithFormLocal<Ring> SmithForm;
#if 0
	{
	List L;

	T.start();
	SmithForm( L, A, R );
	T.stop();

	list<pair<Ring::Element, size_t> > p;

	distinct(L.begin(), L.end(), p);

//	for (List::iterator i = L.begin(); i != L.end(); ++i) cout << *i << ",";
//	cout << endl;

	cout << "# Local2_32, n = " << n << ", m = " << m << endl;
	cout << "smith diag: "; display(p.begin(), p.end());
	T.print(cout << "# time = ") << endl;
	}
#endif
// handle next case as 2A + A^2
	{
	SparseMatrix<Ring> B(R,m,m), C(R,m,m);

	Ring::Element nn; R.init(nn, n);
	for (size_t i = 0; i < m; ++i) {
		B.setEntry(i,i,n);
		for (size_t k = 1; k < m; k<<=1) {
			size_t j = (i & k) ? i - k : i + k;
			B.setEntry(i,j,R.mOne);
		}
	}
	A.zero();
	MatrixDomain<Ring> MD(R);
	MD.mul(A, B, B);
	MD.addin(A,B);
	MD.addin(A,B); //A = 2B + B^2

	List L;

	T.start();
	SmithForm( L, A, R );
	T.stop();

	list<pair<Ring::Element, size_t> > p;

	distinct(L.begin(), L.end(), p);

//	for (List::iterator i = L.begin(); i != L.end(); ++i) cout << *i << ",";
//	cout << endl;

	cout << "# Local2_32, n = " << n+1 << ", m = " << 2*m << ", faked" << endl;
	cout << "smith diag: "; display(p.begin(), p.end());
	T.print(cout << "# time = ") << endl;
	}

#endif
	return 0 ;
}


///// end krat ////////////////////////////
//! @bug this already exists elsewhere
template<class I1, class Lp>
void distinct (I1 a, I1 b, Lp& c)
{
	typename iterator_traits<I1>::value_type e;
	size_t count = 0;
	if (a != b) {e = *a; ++a; count = 1;}
	else return;
	while (a != b)
	{  if (*a == e) ++count;
		else
		{ c.push_back(typename Lp::value_type(e, count));
			e = *a; count = 1;
		}
		++a;
	}
	c.push_back(typename Lp::value_type(e, count));
	return;
}

template <class I>
void display(I b, I e)
{ cout << "(";
	for (I p = b; p != e; ++p) cout << p->first << " " << p->second << ", ";
	cout << ")" << endl;
}

//@}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

