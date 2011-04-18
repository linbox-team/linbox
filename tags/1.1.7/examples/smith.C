/** 
 * examples/smith.C
 *
 * Copyright (C) 2005, 2010  D. Saunders, Z. Wang, J-G Dumas
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */

/** \file examples/smith.C
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

For the "2local" method, the computaattion is done mod 2^32.

For the "local" method, the modulus must be a prime power.

For the "ilio" method, the modulus may be arbitrary composite.
If the modulus is a multiple of the integer determinant, the intege Smith form is obtained.  Determinant plus ilio may be best for smaller matrices.

This example was used during the design process of the adaptive algorithm.
*/

#include <iostream>
#include <string>
#include <vector>
#include <list>

using namespace std;


#include "linbox/field/modular-int32.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/algorithms/smith-form-sparseelim-local.h"

#include "linbox/util/timer.h"

#include "linbox/field/unparametric.h"
#include "linbox/field/local2_32.h"
#include "linbox/field/PIR-modular-int32.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/algorithms/smith-form-local2.h"
#include <linbox/algorithms/smith-form-iliopoulos.h>
#include "linbox/algorithms/smith-form-adaptive.h"
#include "linbox/blackbox/dense.h"

using namespace LinBox;

// #ifndef BIG


template<class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int n, 
		 string src, string file, string format);

template<class I1, class Lp> void distinct (I1 a, I1 b, Lp& c);
template <class I> void display(I b, I e);

int main(int argc, char* argv[]) {
	typedef PIRModular<LinBox::int32> PIR;

	if (argc < 5) {
	
		cout << "usage: " << argv[0] << " alg m n source format \n"  << endl;

	    cout << "alg = `adaptive', `ilio', `local', or `2local', \n"
		     << "m is modulus (ignored by 2local, adaptive), "
			 << "n is matrix order, \n" 
			 << "source is `random', `random-rough', `fib', `tref', or a filename \n"
			 << "format is `dense' or `sparse' (if matrix from a file)\n"
			 << "compile with -DBIG if you want big integers used.\n";

	     return 0;
	}

	string algo = argv[1];

	int m = atoi(argv[2]);

    int n = atoi(argv[3]);

	string src = argv[4];

	string file = src;

	string format = (argc >= 6 ? argv[5] : "");

	UserTimer T;

	if (algo == "adaptive")
	{   
        typedef PID_integer Ints;
		Ints Z;
	    DenseMatrix<Ints> M(Z);
	    Mat(M, Z, n, src, file, format);
		vector<integer> v(n);
	    T.start();
		SmithFormAdaptive::smithForm(v, M);
	    T.stop();
	    list<pair<integer, size_t> > p;

	    distinct(v.begin(), v.end(), p);

	    cout << "#";

	    display(p.begin(), p.end());

	    cout << "# adaptive, Ints, n = " << n << endl;

	    cout << "T" << n << "adaptive" << m << " := ";

	}
	else if (algo == "ilio") { 

		PIR R(m);

	    DenseMatrix<PIR> M(R);

	    Mat(M, R, n, src, file, format);

	    T.start();

	    SmithFormIliopoulos::smithFormIn (M);

	    T.stop();

	    typedef list< PIR::Element > List;

	    List L;

	    for (size_t i = 0; i < M.rowdim(); ++i) L.push_back(M[i][i]);

	    list<pair<PIR::Element, size_t> > p;

	    distinct(L.begin(), L.end(), p);

	    cout << "#";

	    display(p.begin(), p.end());

	    cout << "# ilio, PIR-Modular-int32(" << m << "), n = " << n << endl;

	    cout << "T" << n << "ilio" << m << " := ";
	} 

	else if (algo == "local") { // m must be a prime power
	
           if (format == "sparse" ) {
                typedef Modular<LinBox::int32> Field;
                Field F(m);
                std::ifstream input (argv[4]);
                if (!input) { std::cerr << "Error opening matrix file: " << argv[1] << std::endl; return -1; }

                MatrixStream<Field> ms( F, input );
                SparseMatrix<Field, Vector<Field>::SparseSeq > B (ms);
                std::cout << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;
                if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(std::cout) << std::endl;


                
                Integer p(m), im(m);
                    // Should better ask user to give the prime !!!
                for(unsigned int k = 2; ( ( ! ::probab_prime(p) ) && (p > 1) ); ++k)
                    ::root( p, im, k );

                    // using Sparse Elimination
                LinBox::PowerGaussDomain< Field > PGD( F );
                std::vector<std::pair<size_t,size_t> > local;

                PGD(local, B, m, (int)p);    

                typedef list< Field::Element > List;
				List L;
				for ( std::vector<std::pair<size_t,size_t> >::iterator 
					  p = local.begin(); p != local.end(); ++p) {
					for(size_t i = 0; i < p->first; ++i) L.push_back(p->second);
				}
				size_t m = (B.rowdim() > B.coldim() ? B.coldim() : B.rowdim());
				for (size_t i = L.size(); i < m; ++i) L.push_back(0);

                list<pair<Field::Element, size_t> > pl;
                
                distinct(L.begin(), L.end(), pl);
                
                std::cout << "#";
                
                //display(local.begin(), local.end());
                display(pl.begin(), pl.end());
           		cout << "# local, PowerGaussDomain<int32>(" << m << "), n = " << n << endl;
           
            } else {

                PIR R(m);
		
                DenseMatrix<PIR> M(R);
		
                Mat(M, R, n, src, file, format);
                
                typedef list< PIR::Element > List;
                
                List L;
                
                SmithFormLocal<PIR> SmithForm;
                
                T.start();
                
                SmithForm( L, M, R );
                
                T.stop();
                
                list<pair<PIR::Element, size_t> > p;
                
                distinct(L.begin(), L.end(), p);
                
                cout << "#";
                
                display(p.begin(), p.end());
                
           		cout << "# local, PIR-Modular-int32(" << m << "), n = " << n << endl;
           
            }
           	cout << "T" << n << "local" << m << " := ";
	}

	else if (algo == "2local") { 

		Local2_32 R;

	    DenseMatrix<Local2_32> M(R);

	    Mat(M, R, n, src, file, format);

	    typedef list< Local2_32::Element > List;

	    List L;

	    SmithFormLocal<Local2_32> SmithForm;

	    T.start();

	    SmithForm( L, M, R );

	    T.stop();

	    list<pair<Local2_32::Element, size_t> > p;

	    distinct(L.begin(), L.end(), p);

	    cout << "#";

	    display(p.begin(), p.end());

	    cout << "# 2local, Local2_32, n = " << n << endl;

	    cout << "T" << n << "local2_32 := ";
	}

	else {

		printf ("Unknown algorithms\n");

		exit (-1);

	}

	T.print(cout); cout << ";" << endl;

}

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
   mat read from named file with format "sparse" or "dense".
  Also "tref" and file with format "kdense"
*/
template <class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int n, 
			string src, string file, string format) {

	typename PIR::Element zero; R.init(zero, 0);
    if (src == "random-rough") RandomRoughMat(M, R, n);

    else if (src == "random") RandomFromDiagMat(M, R, n);

    else if (src == "fib") RandomFibMat(M, R, n);

    else if (src == "tref") TrefMat(M, R, n);

    else // from file
	{

		int rdim, cdim;

		std::ifstream in (file.c_str(), std::ios::in);
		if (! in) { cerr << "error: unable to open file" << endl; exit(-1); }

		in >> rdim >> cdim;

		M. resize (rdim, cdim);

		integer val;

		if (format == "dense" ) {

			for (int i = 0; i < rdim; ++ i)

				for ( int j = 0; j < cdim; ++ j) {

					in >> val;

					R. init (M[i][j], val);

				}
		}

		else if (format == "sparse") {

			int i, j;

			char mark;

			in >> mark;

			LinBox::integer val;

			do {

				in >> i >> j;
				in. ignore (1);
				in >> val;

				if ( i == 0) break;

				R. init (M[i-1][j-1], val);

			} while (true);

		}
		  //Krattenthaler's q^e matrices, given by exponent
		else if (format == "kdense") KratMat(M, R, n, in);

		else {

			cout << "Format: " << format << " Unknown\n";

			exit (-1);

		}
	}

    /*show some entries
	for (int k = 0; k < 10; ++k)
	cout << M.getEntry(0,k) <<  " " << M.getEntry(M.rowdim()-1, M.coldim()-1 - k) << endl;
	cout << endl << M.rowdim() << " " << M.coldim() << endl;
	*/

	/* some row ops and some col ops */
} // Mat

// This mat will have s, near sqrt(n), distinct invariant factors, 
// each repeated twice), involving the s primes 101, 103, ...
template <class PIR>
void RandomRoughMat(DenseMatrix<PIR>& M, PIR& R, int n) {
	typename PIR::Element zero; R.init(zero, 0);
	M.resize(n, n, zero);
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
		{   M[i][i] = v; ++i;
		    if (i < n) {M[i][i] = v; ++i;}
		}
	}
    scramble(M);
}

// This mat will have the same nontrivial invariant factors as
// diag(1,2,3,5,8, ... 999, 0, 1, 2, ...).
template <class PIR>
void RandomFromDiagMat(DenseMatrix<PIR>& M, PIR& R, int n) {
	typename PIR::Element zero; R.init(zero, 0);
	M.resize(n, n, zero);

	for (int i= 0 ; i < n; ++i) 
	
		R.init(M[i][i], i % 1000 + 1);
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
	typename PIR::Element zero; R.init(zero, 0);
	M.resize(n, n, zero);

	typename PIR::Element one; R.init(one, 1);

	for (int i= 0 ; i < n; ++i) M[i][i] = one;

	int j = 1, k = 0;

	for (int i= 0 ; i < n-1; ++i) { 

		if ( i == k) {
		
			M[i][i+1] = zero;
			
			k += ++j;
		}

    	else { 
		
			M[i][i+1] = one; 
			
			R.negin(one);
		}
    	R.neg(M[i+1][i], M[i][i+1]);
	}
    scramble(M);
}
	
template < class Ring >
void scramble(DenseMatrix<Ring>& M)
{
	
	    Ring R = M.field();

		int N,n = M.rowdim(); // number of random basic row and col ops.
		N = n;
	
		for (int k = 0; k < N; ++k) {

	    	int i = rand()%M.rowdim(); 
			
	    	int j = rand()%M.coldim(); 
			
	    	if (i == j) continue;

		    // M*i += alpha M*j and Mi* += beta Mj

	   		//int a = rand()%2;
			int a = 0;

	   	 	for (size_t l = 0; l < M.rowdim(); ++l) {

				if (a)

					R.subin(M[l][i], M[l][j]);

				else 

					R.addin(M[l][i], M[l][j]);

				//K.axpy(c, M.getEntry(l, i), x, M.getEntry(l, j));
				//M.setEntry(l, i, c);
   	    	}

	    	//a = rand()%2;

	    	for (size_t l = 0; l < M.coldim(); ++l) {

				if (a)

					R.subin(M[i][l], M[j][l]);
				else 

					R.addin(M[i][l], M[j][l]);
   	    	}
		}

		std::ofstream out("matrix", std::ios::out);

		out << n << " " << n << "\n";

		for (int i = 0; i < n; ++ i) {

			for ( int j = 0; j < n; ++ j) {

				R. write(out, M[i][j]);

				out << " ";
			}

			out << "\n";

		}

	//}
}

//////////////////////////////////
// special mats tref and krat

// Trefethen's challenge #7 mat (primes on diag, 1's on 2^e bands).
template <class PIR>
void TrefMat(DenseMatrix<PIR>& M, PIR& R, int n) {
	typename PIR::Element zero; R.init(zero, 0);
	M.resize(n, n, zero);

	std::vector<int> power2;

	int i = 1;

	do {

		power2. push_back(i);

		i *= 2;
	} while (i < n);

	std::ifstream in ("prime", std::ios::in);

	for ( i = 0; i < n; ++ i)

		in >> M[i][i];

	std::vector<int>::iterator p;

	for ( i = 0; i < n; ++ i) {

		for ( p = power2. begin(); (p != power2. end()) && (*p <= i); ++ p)
			M[i][i - *p] = 1;

		for ( p = power2. begin(); (p != power2. end()) && (*p < n - i); ++ p)
			M[i][i + *p] = 1;
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
        for (int i = m.size(); i <= e; ++i) m.push_back(m[1]*m[i-1]);
		return m[e];
	}
};

// Read "1" or "q" or "q^e", for some (small) exponent e.
// Return value of the power of q at q = _q.
template <class num>
num& qread(num& val, pwrlist& M, istream& in)
{	
	char c;
	in >> c; // next nonwhitespace
	if (c == '0') return val = 0;
	if (c == '1') return val = 1;
	if (c != 'p' && c != 'q') { cout << "exiting due to unknown char " << c << endl; exit(-1);}
	in.get(c); 
	if (c !='^') {in.putback(c); return val = M[1];}
	else 
	{ int expt; in >> expt; 
	  return val = M[expt];
	};
}	

template <class PIR>
void KratMat(DenseMatrix<PIR>& M, PIR& R, int q, istream& in) 
{
	pwrlist pwrs(q); 
	for (unsigned int i = 0; i < M.rowdim(); ++ i)

		for ( unsigned int j = 0; j < M.coldim(); ++ j) {
			int val; 
			qread(val, pwrs, in);
			R. init (M[i][j], val);
		}
}

///// end krat ////////////////////////////
template<class I1, class Lp>
void distinct (I1 a, I1 b, Lp& c)
{ typename I1::value_type e;
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
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
