#include <algorithm>

typedef unsigned int size_t;

/*
Argument meanings:

m = number of rows.
n = number of cols.
s = number of entries specified in input data file (any other matrix entries are implicitly zero).
d = number of bits in largest entry, approximately max_v(1 + log_2(abs(v))).
r = rank, 0 <= r <= min(m, n).  Generally unknown to user, but occasionally known a priori.
*/

/* File formats:
SMS format example:
4 6 M
1 1 3
1 3 7
2 2 9
3 1 27
4 3 -100
0 0 0
m = 4, n = 6, s = 5, d = 1 + lg(100), r = 4.
s is from number of lines strictly between '4 6 M' and '0 0 0'.
d from 3 7 9 27 -100

Other basic sparse format:
4 6 5
1 1 3
1 3 7
2 2 9
3 1 27
4 6 -100
parameters as before, s is from third number on line 1.

Dense format:
4 6
3  0 7 0 0  0
0  9 0 0 0  0
27 0 0 0 0  0
0  0 0 0 0 -100
m = 4, n = 6, s = 24, d = 1 + lg(100),
s is from m*n.
*/

// machine specific costs.
bool GoodMachineCosts = false; // true iff machine cost factors have been obtained.
double C; // cost of wordsize modular axpy.
double D; // cost per entry of wordsize modular dot product.
		  // Thus D*n is cost of length n vector dot product, wordsize modular.
double E; // average cost of wordsize modular axpy, when implicit within blas.


void initMachineCosts (/* url Service */) 
// to become a service call getting these values.
{   C = 1/25.42e6; D = 1/20.46e6; E = 0.5e-8; GoodMachineCosts = true;   }

bool isDense (size_t m, size_t n, long s) 
// a kludge for now.
{   return m*n < 3*s;   }

double traceCost (size_t m, size_t n, long s, size_t d, size_t r = 0)
{	// Trace is very fast.  Perhaps just return 1 second, ignoring input.
	double F = 0.5; // trace fudge factor.
	if (m > n) std::swap(m,n); 
	return F*C*d*m;
}

double minpolyCost (size_t m, size_t n, long s, size_t d, size_t r = 0)
// cost estimate in seconds to compute minimum polynomial of mxn integer matrix 
// of rank r and having s nonzero entries of size bounded by d bits.
{
	if ( ! GoodMachineCosts ) initMachineCosts();
	if (m > n) std::swap(m,n); if (r == 0) r = m;
	double F = 1; // minpoly fudge factor.

	return isDense(m,n,s) ? F*E*d*m*n*r : F*d*(C*s + D*n)*r;
}

double detCost (size_t m, size_t n, long s, size_t d, size_t r = 0)
{
	if ( ! GoodMachineCosts ) initMachineCosts();
	if (m > n) std::swap(m,n); if (r == 0) r = m;
	double F = 1; // det fudge factor.

	return isDense(m,n,s) ? F*d*E*m*n*r : F*d*minpolyCost(m, n, s, d, m);
}

double rankCost (size_t m, size_t n, long s, size_t d, size_t r = 0)
{ 
	double F = 1; // rank fudge factor.
	return F*minpolyCost(m,n,s,d,r); 
}

double valenceCost (size_t m, size_t n, long s, size_t d, size_t r = 0)
{ 
	double F = 1; // valence fudge factor.
	return F*minpolyCost(m,n,s,d,r); 
}

double smithformCost (size_t m, size_t n, long s, size_t d, size_t r = 0)
{
	double F = 1; // smithform fudge factor.
	if (m > n) std::swap(m,n); if (r == 0) r = m;
	return F*E*d*m*m*n*r;  // future tuning...
}
/*
double charpolyCost (size_t m, size_t n, long s, size_t d, size_t r = 0)
*/
double solveCost (size_t m, size_t n, long s, size_t d, size_t r = 0)
{
	double F = 1; // solve fudge factor.
	double BBfactor = 6;
	if (m > n) std::swap(m,n); if (r == 0) r = m;
	return isDense(m,n,s) ? F*d*E*m*n*r : BBfactor*F*d*minpolyCost(m, n, s, d, m);
}


#include <iostream>
#include <sstream>
using namespace std;

string prettyIter (double t, int k)
{
	const int m = 60, h = 60*m, d = 24*h;
	ostringstream foo;
	if (t < 1) foo << "1s";
	else if (k > 1)
	{
		if (t < m) foo << floor(t) << "s";
		else if (t < h) foo << floor(t/m) << "m " << prettyIter(t - m*floor(t/m),1);
		else if (t < d) foo << floor(t/h) << "h " << prettyIter(t - h*floor(t/h),1);
		else foo << floor(t/d) << "d " << prettyIter(t - d*floor(t/d),1);
	}
	else 
	{
		if (t < m) foo << ceil(t) << "s";
		else if (t < h) foo << ceil(t/m) << "m";
		else if (t < d) foo << ceil(t/h) << "h";
		else foo << ceil(t/d) << "d";
	}
	return foo.str();
}
string pretty(double t) { return prettyIter (t, 2); }

int main()
{
	cout << "minpoly: "   << pretty (minpolyCost   (2000,1000,10000,4)) << endl;
	cout << "det: "       << pretty (detCost       (2000,1000,10000,4)) << endl;
	cout << "rank: "      << pretty (rankCost      (2000,1000,10000,4)) << endl;
	cout << "valence: "   << pretty (valenceCost   (2000,1000,10000,4)) << endl;
	cout << "trace: "     << pretty (traceCost     (2000,1000,10000,4)) << endl;
	cout << "smithform: " << pretty (smithformCost (2000,1000,10000,4)) << endl;
	cout << "solve: "     << pretty (solveCost     (2000,1000,10000,4)) << endl;
}
