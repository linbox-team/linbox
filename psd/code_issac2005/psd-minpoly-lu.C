/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* lieatlas/codes/psd-minpoly.C

Read in Lie group model number and compute minpoly of symmetric integer matrix 
specified by a model for E8 and facet e8test.m.

04dec28:
This code is apparently not working currently, since we expect 
nonsingular pos def operators for this facet..  The minpoly for Q2 is
correct, but for others it is singular and/or with all or most 
coefficients positive (expected to alternate).  most probable
bug(s) in cra code ??

*/
// Assuming all principal minor are non-singular
#include <linbox-config.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <linbox/field/gmp-rational.h>
#include <linbox/field/modular-double.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/blackbox/dense.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/fflapack/fflapack.h>
#include <linbox/field/gmp-integers.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/util/timer.h>
#include "cra.h"
#include "read.h"
#include "util.h"

using namespace std;
using namespace LinBox;

int main(int argc, char* argv[]) {   
	std::ostream cnull(0);
	commentator. setBriefReportStream(cnull);
	if (argc < 2 || argc > 3) {	
		cerr << "Usage: " << argv[0] << " <n> [check]" << endl
		     << "returns to stdout the minpoly of the operator defined by " 
		     << "facets/e8test.m, rootData/e8data.m, and models/Q<n>Cpt."
		     << "If `check' is given, the minpoly is checked for correctness"
		     << endl;
		return -1;
	}
	typedef GMPRationalField Rationals;
	typedef vector<Rationals::Element> Q_Vector;
	typedef SparseMatrix<Rationals> Blackbox;
	typedef LieMatrix<Rationals, Blackbox> LBlackbox;

	Rationals Q; 
	LinBox::VectorDomain<Rationals> QVD(Q);
	int rk = 8;
	vector<Q_Vector> root(rk);
	Q_Vector rho(rk);
	Q_Vector nu(rk);
	vector<Blackbox*> op(rk);
	string model = "Q"; model += argv[1]; model += "Cpt";
	Blackbox* F;
	read_lie (nu, root, rho, op, F, model.c_str(), "e8data.m", "e8test.m", Q, rk);
	int dim = op[0]->rowdim();
	clog << "model is " << model << ", dimension " << dim << endl;
	clog << "facet is ";
	QVD.write(clog, nu) << endl;
	Q_Vector r; vector<int> i_l;
	compute_r (r, i_l, nu, root, rho, Q);
	LBlackbox* Lbb;

	buildLieMatrixGMP(Lbb, op, i_l, r, F);
	integer lcm_y = 1, tmp;
	lcmDenEntries (lcm_y, *Lbb);
	clog << "Denominator is " << lcm_y << endl; 
	integer rgcd_M;
	gcdNumEntries (rgcd_M, *Lbb);
	clog << "Gcd of Nnumerator: " << rgcd_M;
	clog << std::endl;
	Rationals::Element q_lcm, q_gcd;
	Q. init (q_lcm, lcm_y); Q. init (q_gcd, rgcd_M);
	Q. invin (q_gcd);
	Lbb -> scalarMulIn (q_lcm);
	Lbb -> scalarMulIn (q_gcd);

	GMP_Integers Z;
	int n = Lbb -> rowdim();

	std::cout << "Matrix order: " << Lbb -> rowdim() << std::endl;
	UserTimer utimer;
	utimer.start();
	DenseMatrix<GMP_Integers> M (Z, n, n);
	DenseMatrix<GMP_Integers>::RawIterator raw_p;
	Rationals::Element q_one, q_zero;
	integer q_tmp;
	Q. init (q_one, 1); Q. init (q_zero, 0);
	Q_Vector y(Lbb -> rowdim()), x (Lbb -> coldim(), q_zero);
	Q_Vector::iterator y_p;

	raw_p = M. rawBegin();
	for (int i = 0; i < Lbb -> rowdim(); ++ i) {
		Q. assign (x[i], q_one);
		Lbb -> apply (y, x);
		for (y_p = y. begin(); y_p != y. end(); ++ y_p, ++ raw_p) {
			Q. get_num (*raw_p, *y_p);
			Q. get_den (q_tmp, *y_p);
			if (q_tmp != 1) {
				std::cout << "Error:\n";
				return 1;
			}
		}
		Q. assign (x[i], q_zero);
	}

	utimer.stop();
	std::cout <<" User time building the matrix: "
			  << utimer.time() << std::endl;

	delete Lbb;

	typedef Modular<double> Field;
	typedef vector<Field::Element> K_Vector;
	integer mmodulus;
	FieldTraits<Field>::maxModulus(mmodulus);
	long bit1 = (long) floor (log((double)mmodulus)/M_LN2);
	long bit2 = (long) floor (log(sqrt(double(4503599627370496LL/n)))/M_LN2);
	RandomPrime primeg(bit1 < bit2 ? bit1 : bit2); long prime;

	utimer.start();
	// modular images
	bool first_time = true;
	Field::Element* FA = new Field::Element[n*n];
	size_t* P= new size_t[n], *PQ = new size_t[n];
	size_t* P_p, * PQ_p;
	vector<integer> v(n);
	Field::Element* p; Field::Element f_tmp;
	int i = 0, j = 0; bool goodPrime;
	CRA<integer> cra; integer m = 1;
	while (! cra.terminated() ){
		// get a prime. Get minpoly mod that prime. Accumulate into v with CRA. 
		primeg. randomPrime(prime);
		if (m % prime == 0) {
			std::clog << "Repeated prime.\n";
			continue;
		}
		m *= prime;
		Field K(prime); 
		//clog << "Computing blackbox matrix mod " << prime;
		for (p = FA, raw_p = M. rawBegin(); p != FA + (n*n); ++ p, ++ raw_p)
			K. init (*p, *raw_p);

		//clog << "\rComputing min poly mod " << prime << ". ";
		FFLAPACK::LUdivine( K, FFLAS::FflasNonUnit, n, n, FA, n, P, FFLAPACK::FflapackLQUP, PQ);

		++i; //clog << '\r' << i << ' ' << prime;

		goodPrime = true;
		for ( j = 0, P_p = P, PQ_p = PQ; j < n; ++ j, ++ P_p, ++ PQ_p) 
			if ((*P_p != j) || (*PQ_p != j))  {
				std::clog << "Bad prime.\n";
				goodPrime = false;
			}

		if(!goodPrime) continue;
		
		//KVD. write (clog, p) << endl;
		vector<Field::Element>::iterator pp; vector<integer>::iterator vp;
		K. init (f_tmp, 1);
		for (j = 0, vp = v.begin(); vp != v.end(); ++j, ++vp) {
			K. mulin (f_tmp, *(FA + (j * n + j)));
			*vp = (long int) f_tmp;
		}
		integer Prime = prime;
		cra.step(Prime, v); 
	}
	cra.result(v);
	delete FA; 
	utimer. stop();

	/*
	VectorDomain<GMP_Integers> ZVD(Z);
	clog << endl;
	ZVD.write(clog, v) << endl;
	*/
	std::cout << "Number of prime needed: " << cra. steps() << std::endl;
	std::cout << "User time computing the minimal polynomial: ";
	std::cout << utimer.time() << std::endl;

	std::cout << "Check if it is positive definite: ";
	if (isAllPositive(v)) std::cout << "Yes\n";
	else std::cout << "No:\n";

	return 0;
}
