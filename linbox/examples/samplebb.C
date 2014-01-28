
/*
 * examples/samplebb.C
 *
 * Copyright (C) 2005, 2010 D Saunders
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

/** \file examples/samplebb.C
 * @example  examples/samplebb.C
 * \ingroup examples
 * \brief generate an example matrix with specified frobenius form.
 *
 * samplebb takes options and any number of argument triples denoting companion
 * matrix blocks.
 * For example, the call "samplebb -r 7 2 3 a3 1 1" generates a sparsely
 * randomized matrix (because of the '-r' option) matrix which is similar
 * (because of the two triples '7 2 3' and 'a2 1 1') to a direct sum of 3
 * companion matrices for (x-7)^2 plus one companion matrix for x^3 + x + 2, the
 * polynomial denoted by 'a3'.
 *
 *  In general, in the first position of each triple 'aK' denotes the polynomial
 *  x^k + x + K-1 and a number n denotes the polynomial x-n.  The second number
 *  in the triple specifies a power of the polynomial and the third specifies how
 *  many companion matrix blocks for that power of that polynomial.
 *
 *  Possible options are
 *  -r lightly randomized similarity transform, matrix remains sparse.
 *  -R fully randomized similarity transform, matrix becomes dense.
 *
 *  The matrix is written to standard out in SMS format (triples).
 *
 *  For some other examples:
 *  "samplebb  1 1 2  2 1 2  4 1 2  12 1 1  0 1 1" is a 8 by 8 diagonal matrix in smith form,
 *  diag(1,1,2,2,4,4,12,0)
 *
 */

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <linbox/blackbox/direct-sum.h>
#include <linbox/blackbox/companion.h>
#include <linbox/algorithms/matrix-hom.h>
#include <linbox/field/ntl-ZZ.h>
#include <NTL/ZZX.h>

using std::string;
using std::list;
using std::vector;
using LinBox::Companion;
using LinBox::DirectSum;
using LinBox::BlasMatrix;
using LinBox::NTL_ZZ;
using NTL::ZZX;


void stripOptions(int& acp, char* avp[], string& opts, const int ac, char** av)
{
	acp = 0;
	for (int i = 1; i < ac; ++i)
	{
		//std::cout << av[i] << " ";
		if (av[i][0] == '-') {
			for (const char* j = av[i]+1; *j != 0; ++j)
				opts.push_back(*j);
		}

		else {
			avp[acp] = av[i];
			++acp;
		}
	}
}

template <class List, class Ring>
void augmentBB(List& L, char* code, int e, int k, const Ring& R)
{
	typedef typename Ring::Element Int;
	Int a, one, zero;
	R.init(one, 1);
	R.init(zero, 0);
	ZZX p;

	// build poly p

	if ( *code != 'a')  // build linear poly
	{
		R.init(a, -atoi(code));
		p += ZZX(0, a);
		p += ZZX(1, one);
	}
	else // build long poly
	{
		int n = atoi(code+1);
		R.init(a, n-1);
		p += ZZX(n, one);
		p += ZZX(1, one);
		p += ZZX(0, a);
	}

	//std::cout << "(code, e, k) =(" << code << ", " << e << ", " << k << ")" << std::endl;
	//std::cout << "Correspoding poly: " << p << std::endl;
	// compute q =  p^e
	ZZX q(0, one);
	for(int i = 0; i < e; ++i) q *= p;
	//std::cout <<"Polynomial: " << q << std::endl;

	vector<Int>  v(deg(q)+1);
	for (int i = 0; i < v.size(); ++i) v[i] = coeff(q, i);

	// companion matrix of q
	Companion<Ring>* C = new Companion<Ring>(R, v);
	for(int i = 0; i < k; ++i) L.push_back(C);

}

template < class Ring >
void scramble(BlasMatrix<Ring>& M)
{

	Ring R = M.field();

	int N,n = M.rowdim(); // number of random basic row and col ops.
	N = 2*n;

	for (int k = 0; k < N; ++k) {

		int i = rand()%M.rowdim();

		int j = rand()%M.coldim();

		if (i == j) continue;

		// M*i += alpha M*j and Mj* -= alpha Mi*

		typename Ring::Element alpha, beta, x;
		R.init(alpha, rand()%5 - 2);
		R.neg(beta, alpha);

		for (size_t l = 0; l < M.rowdim(); ++l) if (!R.isZero(alpha)) {
			R.mul(x, alpha, M[l][j]);
			R.addin(M[l][i], x);
		}


		for (size_t l = 0; l < M.rowdim(); ++l) if (!R.isZero(alpha)) {
			R.mul(x, beta, M[i][l]);
			R.addin(M[j][l], x);
		}
	}

	/*
	   std::ofstream out("matrix", std::ios::out);

	//M. write(std::cout);

	out << n << " " << n << "\n";

	for (int i = 0; i < n; ++ i) {

	for ( int j = 0; j < n; ++ j) {

	R. write(out, M[i][j]);

	out << " ";
	}

	out << "\n";

	}
	*/

}

template <class Matrix>
void printMatrix (const Matrix& A)
{
	int m = A. rowdim();
	int n = A. coldim();
	typedef typename Matrix::Field Ring;
	typedef typename Ring::Element Element;
	std::vector<Element> x(m), y(n);
	Ring r = A. field();
	Element one, zero;
	r. init (one, 1);
	r. init (zero, 0);

	std::cout << m << " " << n <<  " M" << std::endl;
	typename std::vector<Element>::iterator y_p;
	for (int i = 0; i < m; ++ i) {
		r. assign (x[i], one);
		A. applyTranspose(y, x);
		for(y_p = y. begin(); y_p != y. end(); ++ y_p)
			if (! r.isZero(*y_p))
				std::cout << i+1 << " " << y_p - y.begin() + 1 << " " << *y_p << std::endl;
		r. assign (x[i], zero);
	}
	std::cout << "0 0 0" << std::endl;
}


int main(int ac, char* av[])
{
	if (ac < 2)
	{	std::cout << "usage: " << av[0] <<
		" options block-groups." << std::endl;
		std::cout << av[0] << " -r 1 2 3 a4 1 1" << std::endl;
		std::cout <<
		"for lightly randomized matrix similar to direct sum of 3 copies of companion " << std::endl
		<< "matrix of (x-1)^2 and one copy of companion matrix of (x^4 + x + 3)^1." << std::endl;
	}

	typedef NTL_ZZ Ring;
	Ring Z;
	typedef Companion<Ring> BB;
	int acp; char* avp[ac];
	string opts;
	stripOptions(acp, avp, opts, ac, av);
	//std::cout << "number of triples: " << acp << std::endl;
	//for (int i = 0; i < acp; ++ i)
	//	std::cout << avp[i];
	//std::cout << std::endl;
	//std::cout << "Begin to ....\n";
	list<BB*> L;

	for (int i = 0; i < acp; i += 3)
		augmentBB(L, avp[i], atoi(avp[i+1]), atoi(avp[i+2]), Z);

	DirectSum<BB> A(L);
	//std::cout <<"Option: " << opts.c_str() << std::endl;

	if (opts.size() >= 1)
	{	if (opts[0] == 'r')
		{
			// into sparse matrix, then 3n row ops with corresponding col ops
			BlasMatrix<Ring>* B;//(Z,A.rowdim(), A.coldim());
			//MatrixDomain<Ring> MD(Z);
			LinBox::MatrixHom::map (B, A, Z);

			scramble(*B);
			printMatrix(*B);
			delete B;
		}

		if (opts[0] == 'R') ;
		// into dense matrix, then many row ops
		//...

	}
	else {
		printMatrix (A);
	}

	return 0 ;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

