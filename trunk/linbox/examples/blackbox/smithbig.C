/** @name examples/blackbox/ju.C
 * @author bds
 *
 * @doc
 */
//@{

#include <iostream>
#include <string>
#include <list>

using namespace std;
#include "linbox/util/timer.h"
#include "linbox/field/unparametric.h"
#include "linbox/field/local2_32.h"
//#include "linbox/field/PIR-modular-int32.h"
#include "linbox/field/PIR-ntl-ZZ_p.h"
#include "linbox/algorithms/2local-smith.h"
#include "linbox/algorithms/local-smith.h"
#include "linbox/algorithms/iliopoulos-elimination.h"
#include "linbox/blackbox/dense.h"

using namespace LinBox;
//typedef PIRModular<int32> PIR;
typedef PIR_ntl_ZZ_p PIR;

template<class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int n);

template<class I1, class Lp> void distinct (I1 a, I1 b, Lp& c);
template <class I> void display(I b, I e);

int main(int argc, char* argv[])
{
        if (argc < 3)
	{    cout << "usage: " << argv[0] << " n alg m"  << endl;
	     cout << "for n x n matrix, alg = ilio or local or 2local, modulus m" << endl;
	     return 0;
	}
        int n = atoi(argv[1]);
	//int m = (argc == 4) ? atoi(argv[3]) : 0;
	NTL::ZZ m;
	if (argc == 4) NTL::conv (m, argv[3]);

	UserTimer T;
	string algo = argv[2];

	if (algo == "ilio")
	{   PIR R(m);
	    DenseMatrix<PIR> M(R);
	    Mat(M, R, n);
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
	else if (algo == "local") // m must be a prime power
	{   PIR R(m);
	    DenseMatrix<PIR> M(R);
	    Mat(M, R, n);
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
	else if (algo == "2local")
	{   Local2_32 R;
	    DenseMatrix<Local2_32> M(R);
	    Mat(M, R, n);
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
	T.print(cout); cout << ";" << endl;

}

template <class PIR>
void Mat(DenseMatrix<PIR>& M, PIR& R, int n)
{
	M.resize(n, n);

	/* give v about 50 each of 50 primes */
	/*
	int32 p[50] = 
	{2,3,5,7,11,13,17, 19, 23, 29, 
	 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 
	 73, 79, 83, 89, 97, 101, 103, 107, 109, 113,
	 127, 131, 137, 139, 149, 151, 157, 163, 167, 173,
	 179, 181, 191, 193, 197, 199, 211, 223, 227, 229};

	int k = 0;
	for (int i= 0 ; i < 50; ++i)
	  for (int j= 0; j < 50-i ; ++j)
	  {
	    R.init(M[k][k], p[i]); ++k;
	  }
	*/
	for (int i= 0 ; i < n; ++i) M[i][i]=(i+1)%1000;

        /*show some entries*/
//      for (k = 0; k < n-1; ++k)
//	  cout << M.getEntry(k,k) <<  " " << M.getEntry(k, k+1) << " ";
//	cout << M.getEntry(n-1,n-1) <<  endl;

	/* some row ops and some col ops */
	int N = 7*n; // number of random basic row and col ops.
	for (int k = 0; k < N; ++k)
	{
	    int i = rand()%M.rowdim(); 
	    int j = rand()%M.coldim(); 
	    if (i == j) continue;

	    // M*i += alpha M*j and Mi* += beta Mj
	    int a = rand()%2;
	    for (size_t l = 0; l < M.rowdim(); ++l) 
	    {
		if (a)
		R.subin(M[l][i], M[l][j]);
		else 
		R.addin(M[l][i], M[l][j]);
		//K.axpy(c, M.getEntry(l, i), x, M.getEntry(l, j));
		//M.setEntry(l, i, c);
   	    }
	    a = rand()%2;
	    for (size_t l = 0; l < M.coldim(); ++l) 
	    {
		if (a)
		R.subin(M[i][l], M[j][l]);
		else 
		R.addin(M[i][l], M[j][l]);
   	    }
	}

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
