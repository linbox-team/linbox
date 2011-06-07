/* -*- mode:C++ -*- */
/* File: solve-nonsigular.C
	This file was used to generate solver comparison data for the paper "Symbolic-Numeric Exact Rational Linear System Solver" submitted to ISSAC'11
*/

#include <iostream>
#include <fstream>
#include <cstdlib>

#include <linbox/linbox-config.h>

#include <linbox/algorithms/rational-solver.h>
#include <linbox/randiter/random-prime.h>

#include <linbox/field/PID-integer.h>
#include <linbox/field/param-fuzzy.h>
#include <linbox/blackbox/blas-blackbox.h>
#include "tests/test-common.h"
#include "linbox/vector/stream.h"
#include "linbox/util/commentator.h"
#include "linbox/util/timer.h"

#ifdef __LINBOX_HAVE_LAPACK
#include <linbox/algorithms/numeric-solver-lapack.h>
#endif

#ifdef __LINBOX_HAVE_MATLAB
	#include <linbox/algorithms/numeric-solver-matlab.h>
#endif

//or #include "other-numeric-solver.h"
/* a numeric solver is a FAIBB (fast approximate inverse blackbox). It provides
 * 1. constructor from whatever parameters
 * 2. init(A) // init from a matrix of double.  evolve this a bit...
 *    LU or other initial prep may happen at this moment.
 * 3. solve(x, b) // x <-- A^{-1}b approximately, for vector of double x, b.
 * 4. apply(y, x) // y <-- Ax, approximately, for vector of double y, x.
 */

#include <linbox/algorithms/rational-solver-sn.h>
/* rational-solver provides
 * 1. constructor with a numerical solver as argument (call it NS).
 * 2. solve(num, den, A, b, NS)
 *    In our impl, solve prepares the double versions of A, b, initializes the NS,
 *    and calls rsol().
 */

//  matrix types
#include "matrix/coin.h"
#include "matrix/invhilb.h"
#include "matrix/randommat.h"
#include "matrix/randomans.h"
#include "matrix/hadamard.h"
#include "matrix/minmax.h"
#include "matrix/jmat.h"

using namespace LinBox;

enum MatType {diag=0, tref=1, hilb=2, zo=3, rand_sp=4,
				I=5, jordan2=6, rand_near_sing=7, Hadamard=8, minIJ=9,
				maxIJ=10, DlehmerD=11, je1=12, je2=13 };
enum SolverType {diagonal, lapack, matlab, superlu, dixon};

size_t nextPower2(size_t n){
    size_t p = 1;
    while(p < n) p <<= 1;
    return p;
}

template<class Ring, class Matrix, class Vector>
void generateProblem(const Ring& R, Matrix &D, Vector &b,
		LinBox::VectorStream<Vector>& stream1,
		LinBox::VectorStream<Vector>& stream2,
		MatType mt,
		int k = 10)
{

	Vector d, x, y;
	VectorWrapper::ensureDim (d, stream1.n ());
	VectorWrapper::ensureDim (b, stream1.n ());
	VectorWrapper::ensureDim (x, stream1.n ());
	VectorWrapper::ensureDim (y, stream1.n ());

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	int n =(int) d.size();

	bool zeroEntry;
	do {
		stream1.next (d);
		zeroEntry = false;
		for (size_t i=0; i<stream1.n(); i++)
			zeroEntry |= R.isZero(d[i]);
	} while (zeroEntry);

	//  set up RHS
	report << "Setting up RHS... ";
	int randLim = R_CEILING;
 	switch (mt) {
		//  random RHSs
		case zo: //randLim = 10000;
		//case random:
		case I:
		case diag: stream2.next (b);
			//  special case?
			if (n == 4) for (size_t i = 0; i < b.size(); ++i) b[i] = 2*(i+1);
			for (size_t i = 0; i < b.size(); ++i) b[i] %= randLim;
			break;
		//  RHS with just first element 1
		//case zo:
		case rand_sp:
		case rand_near_sing:
		case jordan2:
		case Hadamard:
		case tref:
		case DlehmerD:
		case minIJ:
		case maxIJ:
		case je1:
		case hilb: b[0] = 1; break;
		case je2: b[n-1] = 1; break;
	}
	report << "Done." << endl;

	report << "Setting up matrix order " << n << "... ";
	//  set up Matrix
	typename Ring::Element tmp;
 	switch (mt) {
		case rand_near_sing: randomAns(R, D, n, n); break;
		case hilb: invhilb(R, D, n); break;
		case Hadamard: hadamard(R, D, n); break;
		case minIJ: minmat(R, D, n); break;
		case maxIJ: maxmat(R, D, n); break;
		case DlehmerD: qlehmer(R, D, n); break;
		case je1:
		case je2: jordanform(R, D, n); break;
		case rand_sp:
			randomMat(R, D, n, k); break;
			//  modified for steffy's random model
			/*
			randomMat(R, D, n, n);
			R.init(tmp, 10000);
			for(int i=0; i<n; ++i)
				D.setEntry(i, i, tmp);
			break;
			*/
		case diag:
		  {
		  //typename Ring::Element product;
		  //R.init(product, 1);
		  randLim = 100000;
		  for(int i = 0; i < n; ++i) {
		    int xx = d[i]%randLim;
			if (xx == 0) xx = 1;
			R.init (tmp,  xx);
			//R.mulin(product, tmp);
			//if (n == 4) tmp = i+1;
			//if (tmp == 4) tmp = -4;
			D.setEntry(i, i, tmp);
		  }
		  }
		  break;
		case tref: //trefethen(R, D, n); break;
		case I:
			R.init(tmp, 1);
			for(int i = 0; i < n; ++i)
				D.setEntry(i, i, tmp);
			break;
		case jordan2:
			//randomMat(R, D, n, n);
			for(int i = 0; i < n; ++i){
				R.init(tmp, 1);
				D.setEntry(i, i, tmp);
				R.init(tmp, 0);
				for(int j = i+1; j < n; ++j)
					D.setEntry(i, j, tmp);
				R.init(tmp, 2);
				if (i > 0) D.setEntry(i, i-1, tmp);
			}
			break;
		case zo:
			for(int i = 0; i < n; ++i)
				for(int j = 0; j < n; ++j){
					R.init(tmp, rand()%2);
					D.setEntry(i, j, tmp);
				}
			break;
	}
	report << "Done." << endl;

	stream1.reset ();
   	stream2.reset ();
}


template <class Ring, class RSolver, class Matrix, class Vector>
bool testRandomSolve (const Ring& R, RSolver& rsolver, Matrix& D, Vector &b) {

	int n = (int) b.size();
	Vector d, tmpb, x, y;
	VectorWrapper::ensureDim (d, n);
	VectorWrapper::ensureDim (x, n);
	VectorWrapper::ensureDim (y, n);
	VectorWrapper::ensureDim (tmpb, n);
	VectorDomain<Ring> VD (R);

	for(int i=0; i<n; ++i) tmpb[i] = b[i];

	//std::ostringstream str;
	//std::ostream &report = cerr;
	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	//  print small mats
	if(n <= 20){
		report << "Matrix: " << endl; D.write(report);
		VD.write (report << "Right-hand side:  ", b) << endl;
	}

	std::vector<typename Ring::Element> num(n);
	typename Ring::Element den;
	Timer timer;

	timer.clear(); timer.start();
	int solveResult = rsolver.solve(num, den, D, b);
	timer.stop();

	report << "Total time: " << timer << endl;

	if(n <= 20){
		VD.write (report << "solution numerator: ", num) << endl;
		report << "solution denominator: " << den << endl;
	}

#ifdef WRITE_MATRICES
	/*  write matrices to file */
	ofstream mat;
	stringstream ms;
	ms << n;
	string file = "matrix." + ms.str();
	mat.open(file.c_str());
	D.write(mat);
#endif
#ifdef WRITE_RESULTS
	ofstream out;
	stringstream ss;
	ss << n;
	string res = "output." + ss.str();
	out.open(res.c_str());
	rsolver.writeVec(num, "first value in numerator", 0, 1, out);
	out << endl << endl << "denominator: " << endl << den << endl;
#endif

	if ( solveResult != 0 ) {
	    report << "ERROR: Did not return OK solving status" << endl;
		return false;
	}
	if ( solveResult == 0 && R.isZero(den) ) {
		report << "ERROR: Solver set denominator to zero" << endl;
		return false;
	}
	if ( solveResult == 0 && !VD.areEqual(D.apply(y, num), VD.mulin(tmpb, den)) ) {
		report << "ERROR: Computed solution is incorrect" << endl;
		return false;
	}
	else {
	   return true;
	}
}

int main(int argc, char** argv) {
	bool pass = true;
	int run = 1;
    static size_t n = 10;
	static size_t k = 10;
	bool e = false;
	MatType mt = rand_sp;
	string file;
	SolverType st = lapack;

   static Argument args[] = {
		{ 't', NULL, "(0)diag, (1)tref, (2)ihilb (3)zo, (4)rand_sp\n\t\t(5)I (6)jordan2 (7)rand_near_sing (8)Hadamard\n\t\t(9)minIJ (10)maxIJ (11)DLehmerD (12)Je1 (13)Je2.\n\t\tFor benchmarking use: 2,3,4,6,8,9,10,11, with 6 a special case", TYPE_INT, &mt },
		{ 's', NULL, "Set numerical solver to (0)diagonal, (1)lapack, (2)matlab, (3)superlu", TYPE_INT, &st },
      { 'n', "-n N", "Set order of test matrices to N.", TYPE_INT, &n},
		{ 'k', "-k K", "Set # entries per row to K (for rand_sp case).", TYPE_INT, &k},
		{ 'e', NULL, "Use exact apply", TYPE_BOOL, &e},
		{ 'r', "-r R", "Run solvers with corresponding bit on: numsym(1), zw(2), dixon(4)", TYPE_INT, &run},
		//{ 'f', "-f FILE", "Set input file to FILE.", TYPE_STRING, &file },
		END_OF_ARGUMENTS
		//{ '\0' }
   };
	parseArguments (argc, argv, args);

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	std::ostream &report = commentator.report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

	typedef PID_integer	Ring;  		Ring R;

	typedef ParamFuzzy Field;
	typedef Modular<int32_t> ZField;
	typedef Modular<double> DField;

	typedef BlasBlackbox<Field> Matrix;
	typedef BlasBlackbox<Ring> CommonMatrix;
	typedef vector<Ring::Element> Vector;

	if(mt == Hadamard)
		n = nextPower2(n);

	RandomDenseStream<Ring> s1 (R, n, 1), s2 (R, n, 1);

	CommonMatrix A(R, n, n);
	Vector b(n);
	generateProblem(R, A, b, s1, s2, mt, (int) k);

	if(run & 1){
		/*  choose your numerical solver */
		switch (st){
#ifdef __LINBOX_HAVE_LAPACK
			case lapack: {
				 report << "Using lapack numeric solver." << endl;
				 typedef LPS<Matrix> NumSolver;	NumSolver numSolver;
				 RationalSolverSN<Ring, NumSolver > rsolver(R, numSolver, e);
				 pass = pass && testRandomSolve(R, rsolver, A, b);
				}
				break;
#endif
#ifdef __LINBOX_HAVE_MATLAB
			case matlab: {
				 report << "Using matlab numeric solver." << endl;
				 typedef MLS<Matrix> NumSolver;	NumSolver numSolver;
				 RationalSolverSN<Ring, NumSolver > rsolver(R, numSolver, e);
				 pass = pass && testRandomSolve(R, rsolver, A, b);
				}
				break;
#endif
						 /*
			case superlu: {
				report << "Using SuperLU numeric solver." << endl;
				typedef SLU<Matrix> NumSolver;	NumSolver numSolver(file);
				SNRationalSolver<Ring, NumSolver > rsolver(R, numSolver);
				pass = pass && testRandomSolve(R, rsolver, s1, s2, mt, 1, e, k);
				}
				break;
						  */
			default:
				 break;
		}
		report << "numsym: " << (pass ? "pass" : "fail") << std::endl << std::endl;
	}
	if(run & 2){
		RationalSolver<Ring, ZField, RandomPrimeIterator, WanTraits> rsolver(R);
		pass = testRandomSolve(R, rsolver, A, b);
		report << "zw: " << (pass ? "pass" : "fail") << std::endl << std::endl;
	}
	if(run & 4){
		RandomPrimeIterator genprime( 26-(int)ceil(log((double)n)*0.7213475205) );
		RationalSolver<Ring, DField, RandomPrimeIterator, DixonTraits> rsolver(R, genprime);
		pass = testRandomSolve(R, rsolver, A, b);
		report << "dixon: " << (pass ? "pass" : "fail") << std::endl << std::endl;
	}

	return pass ? 0 : -1;
}
