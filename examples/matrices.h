/*
 * examples/matrices.h
 *
 * Copyright (C) 2017-2019  D. Saunders, Z. Wang, J-G Dumas
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
#ifndef __LinBox_Matrices_H_
#define __LinBox_Matrices_H_

#include <linbox/matrix/dense-matrix.h>

template < class Ring >
void scramble(LinBox::DenseMatrix<Ring>& M)
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
void RandomRoughMat(LinBox::DenseMatrix<PIR>& M, PIR& R, int n) {
	if (n > 10000) {std::cerr << "n too big" << std::endl; exit(-1);}
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
void RandomFromDiagMat(LinBox::DenseMatrix<PIR>& M, PIR& R, int n) {

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
void RandomFibMat(LinBox::DenseMatrix<PIR>& M, PIR& R, int n) {

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
// special matrices tref and krat

#include <givaro/givintprime.h>

// Trefethen's challenge #7 mat (primes on diag, 1's on 2^e bands).
template <class PIR>
void TrefMat(LinBox::DenseMatrix<PIR>& M, PIR& R, int n) {

	std::vector<int> power2;

	int i = 1;

	do {

		power2. push_back(i);

		i *= 2;
	} while (i < n);

    Givaro::IntPrimeDom IPD; Givaro::Integer prime(1);

    for ( i = 0; i < n; ++ i) {
        R.init( M[(size_t)i][(size_t)i], IPD.nextprimein(prime) );
    }

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
	std::vector<Givaro::Integer> m;
	pwrlist(Givaro::Integer q)
        { m.push_back(1); m.push_back(q);
//cout << "pwrlist " << m[0] << " " << m[1] << endl;
        }
	Givaro::Integer operator[](int e)
        {
            for (int i = (int)m.size(); i <= e; ++i) m.push_back(m[1]*m[(size_t)i-1]);
            return m[(size_t)e];
        }
};

// Read "1" or "q" or "q^e", for some (small) exponent e.
// Return value of the power of q at q = _q.
template <class num>
num& qread(num& Val, pwrlist& M, std::istream& in)
{
	char c;
	in >> c; // next nonwhitespace
	if (c == '0') return Val = 0;
	if (c == '1') return Val = 1;
	if (c != 'p' && c != 'q') { std::cerr << "exiting due to unknown char " << c << std::endl; exit(-1);}
	in.get(c);
	if (c !='^') {in.putback(c); return Val = M[1];}
	else
	{ int expt; in >> expt;
    return Val = M[expt];
	};
}

template <class PIR>
void KratMat(LinBox::DenseMatrix<PIR>& M, PIR& R, int q)
{
	M.resize((size_t)q, (size_t)q, R.zero);
	pwrlist pwrs(q);
	for (unsigned int i = 0; i < M.rowdim(); ++ i)

		for ( unsigned int j = 0; j < M.coldim(); ++ j) {
			int Val;
			qread(Val, pwrs, std::cin);
			R. init (M[(size_t)i][(size_t)j], Val);
		}
}

///// end krat ////////////////////////////



template <class PIR>
void MolerMat(LinBox::DenseMatrix<PIR>& A, PIR& R, int n)
{
    A.resize((size_t)n, (size_t)n, R.zero);
    typename PIR::Element tmp; R.init(tmp);
std::cout << A.rowdim() << 'x' << A.coldim() << std::endl;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < i; ++j)
            A.setEntry( i, j, R.init(tmp, j-1) );
        A.setEntry(i,i, R.init(tmp, i+1) );
        for (int j = i+1; j < n; ++j)
            A.setEntry( i, j, R.init(tmp, i-1) );
    }
}

//  n-by-n matrix of 0's and 1's defined by
//  A(i,j) = 1, if j = 1 or if i divides j,
//  and A(i,j) = 0 otherwise.
template <class PIR>
void RedhefferMat(LinBox::DenseMatrix<PIR>& A, PIR& R, int n)
{
    A.resize((size_t)n, (size_t)n, R.zero);
    for (int i = 0; i < n; ++i) {
        A.setEntry(i,0, R.one);
        for (int j = 1; j < n; ++j)
            if ( ((j+1)%(i+1)) == 0 )
                A.setEntry( i, j, R.one );
    }
}


#endif
