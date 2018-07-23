/*
 * examples/mats.C
 *
 * Copyright (C) 2017  D. Saunders, Z. Wang, J-G Dumas
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

/** \file examples/mats.C
 * @example  examples/mats.C
 \brief example matrices that were chosen for Smith form testing.
 \ingroup examples

 \author bds & zw

 Various Smith form algorithms may be used for matrices over the
 integers or over Z_m.  Moduli greater than 2^32 are not supported here.
 Several types of example matrices may be constructed or the matrix be read from a file.
 Run the program with no arguments for a synopsis of the command line parameters.

 For the "adaptive" method, the matrix must be over the integers.
 This is expected to work best for large matrices.

 For the "2local" method, the computation is done mod 2^32.

 For the "local" method, the modulus must be a prime power.

 For the "ilio" method, the modulus may be arbitrary composite.
 If the modulus is a multiple of the integer determinant, the integer Smith form is obtained.  
 Determinant plus ilio may be best for smaller matrices.

 This example was used during the design process of the adaptive algorithm.
*/

#include <linbox/linbox-config.h>

#include <iostream>
#include <string>

using namespace std;


//#include <linbox/ring/modular.h>
//#include "givaro/zring.h"
//#include "linbox/integer.h"

#include <linbox/util/timer.h>
#include <linbox/matrix/dense-matrix.h>

// place A: Edit here and at place B for ring change
//#include <linbox/ring/pir-modular-int32.h>
//#include <linbox/ring/pir-ntl-zz_p.h>
#include <linbox/ring/pir-ntl-zz_p.h>

using namespace LinBox;

template <class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int & n, string src) ;

int main(int argc, char* argv[])
{
	if (argc < 3 or argc > 4) {

		cout << "usage: " << argv[0] << " type n [filename]"  << endl;

		cout << " type = `random', `random-rough', `tref', or `fib',"
			 << " and n is the dimension" << endl;
		cout  << " If filename is present, matrix is written there, else to cout." << endl;

		return 0;
	}

	string type = argv[1];

	int n = atoi(argv[2]);

// place B: Edit here and at place A for ring change
	//typedef PIRModular<int32_t> PIR;
	typedef Givaro::ZRing<Integer> PIR;
	PIR R;
	DenseMatrix<PIR> M(R);

	Mat(M,R,n,type);
	if (argc == 4) {
		ofstream out(argv[3]);
		M.write(out) << endl;
	} else {
		M.write(cout) << endl;
	}
}// main

template < class Ring >
void scramble(DenseMatrix<Ring>& M)
{

	Ring R = M.field();

	int N,n = (int)M.rowdim(); // number of random basic row and col ops.
	N = n;

	for (int k = 0; k < N; ++k) {

		int i = rand()%(int)M.rowdim();

		int j = rand()%(int)M.coldim();

		if (i == j) continue;

            // M*i += alpha M*j and Mi* += beta Mj

            //int a = rand()%2;
		int a = 0;

		for (size_t l = 0; l < M.rowdim(); ++l) {

			if (a)

				R.subin(M[(size_t)l][(size_t)i], M[(size_t)l][(size_t)j]);

			else

				R.addin(M[(size_t)l][(size_t)i], M[(size_t)l][(size_t)j]);

                //K.axpy(c, M.getEntry(l, i), x, M.getEntry(l, j));
                //M.setEntry(l, i, c);
        }

            //a = rand()%2;

		for (size_t l = 0; l < M.coldim(); ++l) {

			if (a)

				R.subin(M[(size_t)i][l], M[(size_t)j][l]);
			else

				R.addin(M[(size_t)i][l], M[(size_t)j][l]);
		}
	}
}

// This mat will have s, near sqrt(n), distinct invariant factors,
// each repeated twice), involving the s primes 101, 103, ...
template <class PIR>
void RandomRoughMat(DenseMatrix<PIR>& M, PIR& R, int n) {
	M.resize((size_t)n, (size_t)n, R.zero);
	if (n > 10000) {cerr << "n too big" << endl; exit(-1);}
	int jth_factor[130] =
        {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67,
         71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149,
         151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199, 211, 223, 227, 229,
         233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283, 293, 307, 311, 313,
         317, 331, 337, 347, 349, 353, 359, 367, 373, 379, 383, 389, 397, 401, 409,
         419, 421, 431, 433, 439, 443, 449, 457, 461, 463, 467, 479, 487, 491, 499,
         503, 509, 521, 523, 541, 547, 557, 563, 569, 571, 577, 587, 593, 599, 601,
         607, 613, 617, 619, 631, 641, 643, 647, 653, 659, 661, 673, 677, 683, 691,
         701, 709, 719, 727, 733};

	for (int j= 0, i = 0 ; i < n; ++j)
	{
		typename PIR::Element v; R.init(v, jth_factor[25+j]);
		for (int k = j ; k > 0 && i < n ; --k)
		{   M[(size_t)i][(size_t)i] = v; ++i;
        if (i < n) {M[(size_t)i][(size_t)i] = v; ++i;}
		}
	}
	scramble(M);
}

// This mat will have the same nontrivial invariant factors as
// diag(1,2,3,5,8, ... 999, 0, 1, 2, ...).
template <class PIR>
void RandomFromDiagMat(DenseMatrix<PIR>& M, PIR& R, int n) {
	M.resize((size_t)n,(size_t) n, R.zero);

	for (int i= 0 ; i < n; ++i)

		R.init(M[(size_t)i][(size_t)i], i % 1000 + 1);
	scramble(M);

}

// This mat will have the same nontrivial invariant factors as
// diag(1,2,3,5,8, ... fib(k)), where k is about sqrt(n).
// The basic matrix is block diagonal with i-th block of order i and
// being a tridiagonal {-1,0,1} matrix whose snf = diag(i-1 1's, fib(i)),
// where fib(1) = 1, fib(2) = 2.  But note that, depending on n,
// the last block may be truncated, thus repeating an earlier fibonacci number.
template <class PIR>
void RandomFibMat(DenseMatrix<PIR>& M, PIR& R, int n) {
	M.resize((size_t)n,(size_t) n, R.zero);

	typename PIR::Element pmone; R.assign(pmone, R.one);

	for (int i= 0 ; i < n; ++i) M[(size_t)i][(size_t)i] = R.one;

	int j = 1, k = 0;

	for (int i= 0 ; i < n-1; ++i) {

		if ( i == k) {

			M[(size_t)i][(size_t)i+1] = R.zero;

			k += ++j;
		}

		else {

			M[(size_t)i][(size_t)i+1] = pmone;

			R.negin(pmone);
		}
		R.neg(M[(size_t)i+1][(size_t)i], M[(size_t)i][(size_t)i+1]);
	}
	scramble(M);
}


//////////////////////////////////
// special mats tref and krat

// Trefethen's challenge #7 mat (primes on diag, 1's on 2^e bands).
template <class PIR>
void TrefMat(DenseMatrix<PIR>& M, PIR& R, int n) {
	M.resize((size_t)n, (size_t)n, R.zero);

	std::vector<int> power2;

	int i = 1;

	do {

		power2. push_back(i);

		i *= 2;
	} while (i < n);

	std::ifstream in ("prime", std::ios::in);

	for ( i = 0; i < n; ++ i)

		in >> M[(size_t)i][(size_t)i];

	std::vector<int>::iterator p;

	for ( i = 0; i < n; ++ i) {

		for ( p = power2. begin(); (p != power2. end()) && (*p <= i); ++ p)
			M[(size_t)i][(size_t)(i - *p)] = 1;

		for ( p = power2. begin(); (p != power2. end()) && (*p < n - i); ++ p)
			M[(size_t)i][(size_t)(i + *p)] = 1;
	}

}
//// end tref ///////  begin krat /////////////////////////////

struct pwrlist
{
	vector<integer> m;
	pwrlist(integer q)
        { m.push_back(1); m.push_back(q); //cout << "pwrlist " << m[0] << " " << m[1] << endl;
        }
	integer operator[](int e)
        {
            for (int i = (int)m.size(); i <= e; ++i) m.push_back(m[1]*m[(size_t)i-1]);
            return m[(size_t)e];
        }
};

// Read "1" or "q" or "q^e", for some (small) exponent e.
// Return value of the power of q at q = _q.
template <class num>
num& qread(num& Val, pwrlist& M, istream& in)
{
	char c;
	in >> c; // next nonwhitespace
	if (c == '0') return Val = 0;
	if (c == '1') return Val = 1;
	if (c != 'p' && c != 'q') { cout << "exiting due to unknown char " << c << endl; exit(-1);}
	in.get(c);
	if (c !='^') {in.putback(c); return Val = M[1];}
	else
	{ int expt; in >> expt;
    return Val = M[expt];
	};
}

template <class PIR>
void KratMat(DenseMatrix<PIR>& M, PIR& R, int q)
{
	pwrlist pwrs(q);
	for (unsigned int i = 0; i < M.rowdim(); ++ i)

		for ( unsigned int j = 0; j < M.coldim(); ++ j) {
			int Val;
			qread(Val, pwrs, cin);
			R. init (M[(size_t)i][(size_t)j], Val);
		}
}

///// end krat ////////////////////////////

/** Output matrix is determined by src which may be:
  "random-rough"
  This mat will have s, near sqrt(n), distinct invariant factors,
  each repeated twice), involving the s primes 101, 103, ...
  "random"
  This mat will have the same nontrivial invariant factors as
  diag(1,2,3,5,8, ... 999, 0, 1, 2, ...).
  "fib"
  This mat will have the same nontrivial invariant factors as
  diag(1,2,3,5,8, ... fib(k)), where k is about sqrt(n).
  The basic matrix is block diagonal with i-th block of order i and
  being a tridiagonal {-1,0,1} matrix whose snf = diag(i-1 1's, fib(i)),
  where fib(1) = 1, fib(2) = 2.  But note that, depending on n,
  the last block may be truncated, thus repeating an earlier fibonacci number.
  "file" (or any other string)
  Also "tref" and file with format "kdense"
  */
template <class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int & n,
         string src) {

	if (src == "random-rough") RandomRoughMat(M, R, n);

	else if (src == "random") RandomFromDiagMat(M, R, n);

	else if (src == "fib") RandomFibMat(M, R, n);

	else if (src == "tref") TrefMat(M, R, n);

	else if (src == "krat") KratMat(M, R, n);

	else { // from cin, mostly pointless, but may effect a file format change.

			M.read(cin);
		n = M.rowdim();
	}

} // Mat
