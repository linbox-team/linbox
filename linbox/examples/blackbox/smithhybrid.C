/** \file smith.C examples/blackbox/smith.C
\brief mod m Smith form by elmination

 \author bds & zw

*/
// Consider rank = 0, 1 case.
#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <linbox/integer.h>
#include "linbox/util/timer.h"
#include <linbox/util/commentator.h>
#include "linbox/blackbox/dense.h"
#include "linbox/blackbox/sparse.h"
#include <linbox/algorithms/smith-form-adaptive.h>
#include <linbox/field/ntl-ZZ.h>

using namespace LinBox;
using std::string;

template<class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int n, int m,
		 string src, string file, string format);

template<class I1, class Lp> void distinct (I1 a, I1 b, Lp& c);
template <class I> void display(I b, I e);

template <class Ring>
void compute_vf_dense (std::vector<integer>& L, const DenseMatrix<Ring>& M, const string& alg);

template <class Ring>
void compute_vf_sparse (std::vector<integer>& L, const SparseMatrix<Ring>& M, const string& alg);

int main(int argc, char* argv[]) {

	srand(time(0));
	using namespace std;
	if (argc < 6) {
		cout << "usage: " << argv[0] << " alg n m source format \n";
	    cout << "alg (hybrid, ilio, val, bis, backwardbis) "
			 << "n is the matrix row dimension\n"
			 << "m is matrix column dimension\n"
			 << "source is `random' or `fib' or a filename or random-rough\n"
			 << "format is `dense' or `sparse' (random, random-rough, fib just support dense) \n";
	     return 0;
	}

	string alg = argv[1];
	int n = atoi (argv[2]); int m = atoi (argv[3]);
	string src = argv[4]; string format = argv[5];
	UserTimer T; 
	typedef NTL_ZZ PIR;

	commentator.setBriefReportStream(std::cout);
	commentator.setReportStream (std::cout);
    T.start();
	PIR R;
	// default integer = 0
	std::vector<integer> L;

	if (format == "dense" || alg == "ilio" || alg == "bis" || alg ==  "backwardbis") {
		DenseMatrix<PIR> M(R, n, m);
		Mat(M, R, n, m, src, src, format);
		int order = M. rowdim() <= M.coldim() ? M. rowdim() : M.coldim();
		L. resize (order);
		compute_vf_dense (L, M, alg);
	}
	else if (format == "sparse") {
		SparseMatrix<PIR> M(R);
		std::ifstream in(src.c_str(), std::ios::in);
		M. read (in); in. close ();
		int order = M. rowdim() <= M.coldim() ? M. rowdim() : M.coldim();
		L. resize (order);
		compute_vf_sparse (L, M, alg);
	}
	else {
		std::cout <<"Wrong format\n";
		exit (1);
	}
	T. stop ();
	
	list<pair<integer, size_t> > p;
	distinct(L.begin(), L.end(), p);
	cout << "#";
	display(p.begin(), p.end());
	cout << "T_" << alg << n << "x" << m << src << (format == "dense" ? 'D' : 'S') << " := " ;
	T. print (cout);
	std::cout << "sec\n";

	return 0;
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
 "biglfif" This mat will have 3 or 4 huge largest few invariant factors.
 "file" (or any other string)
   mat read from named file with format "sparse" or "dense".
*/
template <class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int n, int m,
			string src, string file, string format) {

	typename PIR::Element one; 
	R.init(one, 1);
	typename PIR::Element zero; 
	R.init(zero, 0);
    if (src == "random-rough") {
		M.resize(n, m);
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
		M.resize(n, m);
		for (int i= 0 ; i < n; ++i) 
			R.init(M[i][i], i % 1000 + 1);
	    scramble(M);

	}
	
    else if (src == "biglfif") { 
	// big last two invariant factors. Element must be bignums.
		M.resize(n, m);
        int k = (n < m) ? n : m;
		for (int i= 0 ; i < k; ++i) 
			R.init(M[i][i], 1);
	    typename PIR::Element p, t, p4, p5; R.init(t, 10000000000);
		R.mul(p, t, t); R.mulin(p, p); R.mulin(p, t); 
		R.addin(p, R.init(t, 151));  // first prime after 10^50.
		// p = 100000000000000000000000000000000000000000000000151;
		R.mul(t, p, p); R.mul(p4, t, t); R.mul(p5, p, p4);
		R.init(M[k-1][k-1], p5); // p^5
		R.init(M[k-2][k-2], p4); // p^4
	    scramble(M);

	}
	
    else if (src == "fib") {
		M.resize(n, m);
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
		if (!in) { cerr << "error: unable to open file" << endl; exit(-1); }
		in >> n; in >> m; 
		M. resize (n, m);

		if (format == "dense") {

			for (int i = 0; i < n; ++ i)

				for ( int j = 0; j < m; ++ j) {

					in. ignore(1);
					R. read(in, M[i][j]);
					//in >> val;
					//R. init (M[i][j], val);
				}
		}

		else if (format == "sparse") {

			int i, j;
			in >> mark;
			do {

				in >> i >> j;
				if ( i == 0) break;

				in. ignore (1);
				//R. init (M[i-1][j-1], val);
				R. read (in, M[i-1][j-1]);

			} while (true);

		}

		else {

			cout << "Format: " << format << " Unknown\n";

			exit (-1);

		}
	}
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
template<class Matrix>
void compute_hybrid (std::vector<integer>& L, const Matrix& M);
template<class Matrix>
void compute_ilio (std::vector<integer>& L, const Matrix& M);
template<class Matrix>
void compute_val (std::vector<integer>& L, const Matrix& M);
template<class Matrix>
void compute_bis (std::vector<integer>& L, const Matrix& M);

template <class Ring>
void compute_vf_dense (std::vector<integer>& L, const DenseMatrix<Ring>& M, const string& alg) {

	if (alg == "hybrid") {
		compute_hybrid (L, M);
	}
	else if (alg == "ilio") {
		compute_ilio (L, M);
	}
	else if (alg == "bis") {
		compute_bis (L, M);
	}
	else if (alg == "backwardbis") {
		compute_backwardbis (L, M);
	}
	else if (alg == "val") {
		compute_val (L, M);
	}
	else if (alg == "print") {
	    M.write(std::cout);
	}
	else {
		std::cerr << "Unknown algorithm.\n";
		exit (-1);
	}
}

template <class Ring>
void compute_vf_sparse (std::vector<integer>& L, const SparseMatrix<Ring>& M, const string& alg) {

	if (alg == "hybrid") {
		compute_hybrid (L, M);
	}
	else if (alg == "ilio") {
		std::cout << "This should not be called.\n";
		exit (1);
	}
	else if (alg == "bis" || alg == "backwardbis") {
		std::cout << "This should not be called.\n";
		exit (1);
	}
	else if (alg == "val") {
		compute_val (L, M);
	}
	else if (alg == "print") {
	    M.write(std::cout, SparseMatrix<Ring>::FORMAT_GUILLAUME );
	}
	else {
		std::cerr << "Unknown algorithm.\n";
		exit (-1);
	}
}


template <class Matrix>
void compute_hybrid (std::vector<integer>& L, const Matrix& M){
	SmithFormAdaptive::smithForm (L, M);
}

template <class Matrix>
void compute_ilio (std::vector<integer>& L, const Matrix& M) {

		std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		std::vector<integer>::iterator L_p;
		report << "Computation of the rank starts:\n";
		typedef typename Matrix::Field Ring;
		unsigned long r;
		MatrixRank<Ring, Modular<int> > MR;
		r = MR. rank (M);
		report << "   Matrix rank over a random prime field: " << r << '\n';
		report << "Computation of the rank finished.\n";
		
		report << "Computation of the largest invariant factor with bonus starts:\n";
		typedef RationalSolverAdaptive Solver;
		typedef LastInvariantFactor<Ring, Solver> LIF;
		typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
		OIF oif; oif. setThreshold  (10); oif.getLastInvariantFactor().setThreshold (10);
		typename Ring::Element _lif, _bonus; integer lif, bonus;
		oif. oneInvariantFactor_Bonus (_lif, _bonus, M, (int)r);
		M. field(). convert (lif, _lif); M. field(). convert (bonus, _bonus);
		//report << "   The largest invariant factor: " << lif << std::endl;
		//report << "   Bonus (previous one): " << bonus << std::endl;
		//std::cout << "   The largest invariant factor: " << lif << std::endl;
		//std::cout << "   Bonus (previous one): " << bonus << std::endl;
		if (bonus == 1) {
			for (L_p = L. begin(); L_p != L. end(); ++ L_p)
				*L_p = 1;
		}
		else if ( bonus <  PIRModular<int>::getMaxModulus() ) {
			report << "    Elimination starts:\n";
			PIRModular<int> R (bonus);
			DenseMatrix<PIRModular<int> >* M_ilio;
			MatrixMod::mod (M_ilio, M, R);
			IliopoulosElimination::smithIn (*M_ilio);
			int i; 
			for (i = 0, L_p = L. begin(); L_p != L. end(); ++ i, ++ L_p)
			R. convert(*L_p, (*M_ilio) [i][i]);
			delete M_ilio;
			report << "    Elimination ends.\n";
		}
	 	else {
			report << "    Elimination start:\n";
			PIR_ntl_ZZ_p R (bonus);
			DenseMatrix<PIR_ntl_ZZ_p>* M_ilio;
			MatrixMod::mod (M_ilio, M, R);
			IliopoulosElimination::smithIn (*M_ilio);
			int i; 
			for (i = 0, L_p = L. begin(); L_p != L. end(); ++ i, ++ L_p)
				R. convert(*L_p, (*M_ilio) [i][i]);
			delete M_ilio;
			report << "    Elimination ends.\n";
		}
		for (L_p = L. begin(); L_p != L. begin() + r; ++ L_p)
			if (*L_p == 0) *L_p = bonus;
		L[r-1] = lif;
}

template <class Matrix>
void compute_val (std::vector<integer>& L, const Matrix& M) {

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	report << "Computation of the rank starts:\n";
	typedef typename Matrix::Field Ring;
	unsigned long r;
	MatrixRank<Ring, Modular<int> > MR;
	r = MR. rank (M);
	report << "   Matrix rank over a random prime field: " << r << '\n';
	report << "Computation of the rank finished.\n";
	std::vector<long> e(SmithFormAdaptive::NPrime); std::vector<long>::iterator e_p;
	report <<"   Compute the degree of min poly of AA^T: \n"; 
	typedef Modular<int> Field; 
	integer Val; Field::Element v; unsigned long degree; 
	typename MatrixModTrait<Matrix, Field>::value_type* Mp; 
	RandomPrime rg ((int)(log( (double)(Field::getMaxModulus()) ) / M_LN2 - 2)); 
	Field F (rg. randomPrime()); MatrixMod::mod (Mp, M, F); 
	Valence::one_valence (v, degree, *Mp); delete Mp; 
	report <<"   Degree of minial polynomial of AA^T = " << degree << '\n'; 
	// if degree is small if (degree < sqrt(double(order))) {
	report << "   Computation of the valence starts:\n"; 
	Valence::valence (Val, degree, M); 
	report << "      Valence = " << Val << std::endl; 
	report << "   Computation of the valence of ends.\n"; 
	Val = abs (Val); 
	//Factor the k-smooth part of Val 
	const long* prime_p;
	for (prime_p = SmithFormAdaptive::prime, e_p = e. begin(); e_p != e. end(); ++ prime_p, ++ e_p) { 
		*e_p = 0; 
		while (Val % *prime_p == 0) { 
			++ *e_p; Val = Val / *prime_p; 
		}
	}
	if (Val == 1) { 
		SmithFormAdaptive::smithFormVal (L, M, r, e); 
		report << "Computation of the invariant factors ends." << std::endl; return; } 
	else {std::cerr << "   Valence is rough.\n"; exit (-1);}
}

template <class Matrix>
void compute_bis (std::vector<integer>& L, const Matrix& M) {

	typedef typename Matrix::Field Ring;
	typedef Modular<int> Field;
	typedef RationalSolverAdaptive Solver;
	typedef LastInvariantFactor<Ring, Solver> LIF;
	typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
	SmithForm<Ring, OIF, MatrixRank<Ring, Field > > sf;;
	sf. setOIFThreshold (6); sf. setLIFThreshold (10);
	int order = M. rowdim() <= M. coldim() ? M. rowdim() : M. coldim();
	std::vector<typename Ring::Element> out(order);
	sf. smithForm (out, M);
   	typename std::vector<typename Ring::Element>::iterator out_p; std::vector<integer>::iterator L_p;
	for (L_p = L. begin(), out_p = out. begin(); out_p != out. end(); ++ out_p, ++ L_p)
		M. field(). convert (*L_p, *out_p);
}

template <class Matrix>
void compute_backwardbis (std::vector<integer>& L, const Matrix& M) {

	typedef typename Matrix::Field Ring;
	typedef Modular<int> Field;
	typedef RationalSolverAdaptive Solver;
	typedef LastInvariantFactor<Ring, Solver> LIF;
	typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
	SmithForm<Ring, OIF, MatrixRank<Ring, Field > > sf;;
	sf. setOIFThreshold (6); sf. setLIFThreshold (10);
	int order = M. rowdim() <= M. coldim() ? M. rowdim() : M. coldim();
	std::vector<typename Ring::Element> out(order);
	sf. smithFormBackward (out, M);
   	typename std::vector<typename Ring::Element>::iterator out_p; std::vector<integer>::iterator L_p;
	for (L_p = L. begin(), out_p = out. begin(); out_p != out. end(); ++ out_p, ++ L_p)
		M. field(). convert (*L_p, *out_p);
}

		

//@}
