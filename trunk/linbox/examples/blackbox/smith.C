/** \file smith.C examples/blackbox/smith.C
\brief mod m Smith form by elmination

 \author bds & zw

*/
#include <iostream>
#include <string>
#include <list>

using namespace std;
#include "linbox/util/timer.h"
#include "linbox/field/unparametric.h"
#include "linbox/field/local2_32.h"
#include "linbox/field/PIR-modular-int32.h"
#include "linbox/algorithms/2local-smith.h"
#include "linbox/algorithms/local-smith.h"
#include "linbox/algorithms/iliopoulos-elimination.h"
#include "linbox/blackbox/dense.h"

using namespace LinBox;
typedef PIRModular<int32> PIR;

template<class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int n, 
		 string src, string file, string format);

template<class I1, class Lp> void distinct (I1 a, I1 b, Lp& c);
template <class I> void display(I b, I e);

int main(int argc, char* argv[]) {

	if (argc < 5) {
	
		cout << "usage: " << argv[0] << " alg m n source format \n"  << endl;

	    cout << "alg = `ilio' or `local' or `2local', \n"
		     << "m is modulus (ignored by 2local), n is matrix order, \n" 
			 << "source is `random' or `fib' or a filename \n"
			 << "format is `dense' or `sparse' (ignored for random or fib) \n";

	     return 0;
	}

	string algo = argv[1];

	int m = atoi(argv[2]);

    int n = atoi(argv[3]);

	string src = argv[4];

	string file = src;

	string format = (argc >= 6 ? argv[5] : "");

	UserTimer T;

	if (algo == "ilio") { 

		PIR R(m);

	    DenseMatrix<PIR> M(R);

	    Mat(M, R, n, src, file, format);

	    T.start();

	    IliopoulosElimination::smithIn (M);

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
	
		PIR R(m);
		
	    DenseMatrix<PIR> M(R);
		
	    Mat(M, R, n, src, file, format);

	    typedef list< PIR::Element > List;

	    List L;

	    LocalSmith<PIR> SmithForm;

	    T.start();

	    SmithForm( L, M, R );

	    T.stop();

	    list<pair<PIR::Element, size_t> > p;

	    distinct(L.begin(), L.end(), p);

	    cout << "#";

	    display(p.begin(), p.end());

	    cout << "# local, PIR-Modular-int32(" << m << "), n = " << n << endl;

	    cout << "T" << n << "local" << m << " := ";
	}

	else if (algo == "2local") { 

		Local2_32 R;

	    DenseMatrix<Local2_32> M(R);

	    Mat(M, R, n, src, file, format);

	    typedef list< Local2_32::Element > List;

	    List L;

	    LocalSmith<Local2_32> SmithForm;

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

/** 
  random-rough:
   This mat will have s (near sqrt(n)) distinct invariant factors, each repeated
   twice), involving the s primes 101, 103, ...
  random:
   This mat will have the same nontrivial invariant factors as
   diag(1,2,3,5,8, ... 999, 0, 1, 2, ...).
  fib:
   This mat will have the same nontrivial invariant factors as
   diag(1,2,3,5,8, ... fib(k)), where k is about sqrt(n). 
   The basic matrix is block diagonal with i-th block of order i and
   being a tridiagonal {-1,0,1} matrix whose snf = diag(i-1 1's, fib(i)),
   where fib(1) = 1, fib(2) = 2.  But note that, depending on n, 
   the last block may be truncated, thus repeating an earlier fibonacci number.
  file
   mat read from file file with format `sparse or `dense
*/
template <class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int n, 
			string src, string file, string format) {

	M.resize(n, n);

	typename PIR::Element one; 
	
	R.init(one, 1);

	typename PIR::Element zero; 
	
	R.init(zero, 0);

    if (src == "random-rough") {
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
		//cerr << "M is built" << endl;

	}
	
    else if (src == "random") {
	
		for (int i= 0 ; i < n; ++i) 
		
			R.init(M[i][i], i % 1000 + 1);
	    scramble(M);

	}
	
    else if (src == "fib") {

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
    else 
	{

		int n, m;

		char mark;

		std::ifstream in (file.c_str(), std::ios::in);
		if (! in) { cerr << "error: unable to open file" << endl; exit(-1); }

		in >> n;

		in >> m;

		M. resize (n, m);

		int val;

		if (format == "dense") {

			for (int i = 0; i < n; ++ i)

				for ( int j = 0; j < m; ++ j) {

					in >> val;

					R. init (M[i][j], val);

				}
		}

		else if (format == "sparse") {


			int i, j;

			in >> mark;

			do {

				in >> i >> j >> val;

				if ( i == 0) break;

				R. init (M[i-1][j-1], val);

			} while (true);

		}

		else {

			cout << "Format: " << format << " Unknown\n";

			exit (-1);

		}
	}
	/*
	else {

		cout << "Source: " << src << " Unknow choice\n";

		exit (-1);

	}
	*/

    /*show some entries*/
	//for (k = 0; k < n-1; ++k)
	//cout << M.getEntry(k,k) <<  " " << M.getEntry(k, k+1) << " ";
	//cout << M.getEntry(n-1,n-1) <<  endl;

	/* some row ops and some col ops */
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

		//M. write(std::cout);

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
