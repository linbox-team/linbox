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
#include <linbox-config.h>

#include <iostream>
#include <fstream>
#include <vector>

#include <linbox/field/gmp-rational.h>
#include <linbox/field/gmp-integers.h>
#include <linbox/field/modular-int.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/vector/vector-domain.h>
#include <linbox/util/timer.h>
#include "read.h"
#include "util.h"
#include "lie-matrix.h"
#include "minpoly.h"

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

	Rationals Q; 
	LinBox::VectorDomain<Rationals> QVD(Q);
	int rk = 8;
	vector<Q_Vector> root(rk);
	Q_Vector rho(rk);
	Q_Vector nu(rk);
	vector<Blackbox*> op(rk);
	string model = "Q"; model += argv[1]; model += "Cpt";
	Blackbox* F;
	read_lie (nu, root, rho, op, F, 
			  model.c_str(), "e8data.m", "e8test.m", Q, rk);
	int dim = op[0]->rowdim();
	clog << "model is " << model << ", dimension " << dim << endl;
	clog << "facet is ";
	QVD.write(clog, nu) << endl;
	Q_Vector r; vector<int> i_l;
	compute_r (r, i_l, nu, root, rho, Q);

	LieMatrix<Rationals, SparseMatrix<Rationals> >* M;
	buildLieMatrixGMP(M, op, i_l, r, F);
	typedef std::vector<integer> Polynomial;
	typedef Modular<int> Field;
	Polynomial P;

	integer lcm_M, gcd_M;
	lcmDenEntries (lcm_M, *M);
	std::clog << "Denominator: " << lcm_M << std::endl;
	Rationals::Element rlcm_M, rgcd_M; Q. init (rlcm_M, lcm_M);
	M -> scalarMulIn (rlcm_M);
	gcdNumEntries (gcd_M, *M);
	std::clog << "GCD of numerators: " << gcd_M << std::endl;
	Q. init (rgcd_M, gcd_M);
	Q. invin (rgcd_M);
	M -> scalarMulIn (rgcd_M);

	std::cout << "Matrix order: " << M -> rowdim() << std::endl;

	UserTimer timer;
	timer. start ();
	MinPoly<integer, Field>::minPoly (P, *M);
	timer. stop ();

	GMP_Integers Z; VectorDomain<GMP_Integers> ZVD(Z);
	/*
	std::clog << "MinPoly:\n";
	ZVD.write(std::clog, P);
	std::clog << endl;
	*/

	std::cout << "User time: " << timer. time() << std::endl;

	std::cout << "Check if it is positive definite: ";
	if (isAlternativeSign(P)) std::cout << "Yes\n";
	else std::cout << "No:\n";

	if (M -> rowdim() < 100) {
		std::clog << "Check if it is correct. ";
		bool isMinPoly = check_minpoly(*M, P);
		delete M;
		if (isMinPoly) std::clog << "Yes\n";
		else std::clog << "no\n";
	}

	return 0;
}
