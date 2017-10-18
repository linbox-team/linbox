#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
#include <list>

#include "givaro/givtimer.h"

// Polynomial Matrix / Order Basis
#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"
#include "linbox/ring/givaro-poly.h"
#include "linbox/ring/givaro-poly-quotient.h"
#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/algorithms/poly-dixon.h"
#include "linbox/matrix/factorized-matrix.h"
#include "linbox/algorithms/invert-tb.h"
#include "linbox/algorithms/det-tb.h"

//#define LINBOX_USES_OMP 1
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;

#include "test-poly-smith-form.h"

typedef NTL_zz_pX PolyRing;

typedef typename PolyRing::Element Polynomial;

typedef TestPolySmithFormUtil<PolyRing> Util;

typedef MatrixDomain<PolyRing> MatDom;

typedef BlasMatrix<PolyRing> Matrix;

typedef SmithFormKannanBachemDomain<PolyRing> SmithDom;
typedef PolyDixonDomain<PolyRing> DixonDom;

// typedef GivaroPolyQuotient<PolyRing> QuotRing;
// typedef BlasMatrix<QuotRing> QuotMatrix;
// typedef SmithFormKannanBachemDomain<QuotRing> QSmithDom;

Polynomial makeLump(PolyRing &PD, size_t p, size_t d) {
	std::vector<integer> coefs;
	for (size_t i = 0; i <= d; i++) {
		coefs.push_back(rand() % p);
	}
	
	Polynomial tmp;
	PD.init(tmp, coefs);
	return tmp;
}

void solveTextBook(const PolyRing &PD, const Matrix &A) {
	SmithDom SFD(PD);
	
	Matrix B(A);
	
	std::vector<Polynomial> result;
	
	Givaro::Timer TW;
	
	TW.clear();
	TW.start();
	SFD.solveTextBook(result, B);
	TW.stop();
	double tb_time = TW.usertime() + TW.systime();
	std::cout << " " << tb_time;
	
	/*
	std::cout << "Result Text Book:" << std::endl;
	for (size_t i = 0; i < result.size(); i++) {
		PD.write(std::cout, result[i]) << std::endl;
	}
	// */
}

void solveKannanBachem(const PolyRing &PD, const Matrix &A) {
	SmithDom SFD(PD);
	
	Matrix B(A);
	
	std::vector<Polynomial> result;
	
	Givaro::Timer TW;
	
	TW.clear();
	TW.start();
	SFD.solve(result, B);
	TW.stop();
	double kb_time = TW.usertime() + TW.systime();
	std::cout << " " << kb_time;
	
	/*
	std::cout << "Result Kannan-Bachem:" << std::endl;
	for (size_t i = 0; i < result.size(); i++) {
		PD.write(std::cout, result[i]) << std::endl;
	}
	//*/
}

void solveIliopoulos(const PolyRing &PD, const Matrix &A, const Polynomial &det) {
	SmithDom SFD(PD);
	
	Matrix B(A);
	std::vector<Polynomial> result;
	
	Givaro::Timer TW;
	
	TW.clear();
	TW.start();
	SFD.solveIliopoulos(result, B, det);
	TW.stop();
	double ilio_time = TW.usertime() + TW.systime();
	std::cout << " " << ilio_time;
	
	/*
	std::cout << "Result Iliopoulos:" << std::endl;
	for (size_t i = 0; i < result.size(); i++) {
		PD.write(std::cout, result[i]) << std::endl;
	}
	//*/
}

/*
void solveQuotTextBook(const PolyRing &PD, const Matrix &A, const Polynomial &det) {
	QuotRing QD(PD, det);
	QSmithDom QSFD(QD);
	Util util(PD);
	
	QuotMatrix QA(A, QD);
	
	// util.printMatrix(QA);
	
	std::vector<Polynomial> result;
	QSFD.solveTextBook(result, QA);
	
	std::cout << "Quot Text Book Result:" << std::endl;
	for (size_t i = 0; i < result.size(); i++) {
		PD.write(std::cout, result[i]) << std::endl;
	}
}

void solveLocal(const PolyRing &PD, const Matrix &A, const Polynomial &det) {
	QuotRing QD(PD, det);
	SmithFormLocal<QuotRing> SFD;
	Util util(PD);
	
	QuotMatrix QA(A, QD);
	
	// util.printMatrix(QA);
	
	std::list<Polynomial> result;
	SFD(result, QA, QD);
	
	std::cout << "Local Result:" << std::endl;
	for (std::list<Polynomial>::const_iterator iterator = result.begin(), end = result.end(); iterator != end; ++iterator) {
		PD.write(std::cout, *iterator) << std::endl;
	}
}

void factorizeMatrix(const PolyRing &PD, const Matrix &A) {
	std::vector<integer> v;
	v = {1, 1};
	//v = {1, 4, 0, 1};
	Polynomial f;
	PD.init(f, v);
	PD.write(std::cout << "f: ", f) << std::endl;
	
	QuotRing QD(PD, f);
	Util util(QD);
	
	QuotMatrix QA(A, QD);
	
	util.printMatrix(QA);
	
	InvertTextbookDomain<QuotRing> ID(QD);
	
	QuotMatrix Ainv(QD, A.rowdim(), A.coldim());
	
	ID.invert(Ainv, QA);
	
	util.printMatrix(Ainv);
	
	QuotMatrix C(QD, A.rowdim(), A.coldim());
	
	MatrixDomain<QuotRing> QMD(QD);
	QMD.mul(C, QA, Ainv);
	
	util.printMatrix(C);
}
*/

void solveDet(const PolyRing &PD, const Matrix &A) {
	DetTextbookDomain<PolyRing> DD(PD);
	Polynomial det;
	
	Givaro::Timer TW;
	
	TW.clear();
	TW.start();
	DD.solve(det, A);
	TW.stop();
	double kb_time = TW.usertime() + TW.systime();
	std::cout << " " << kb_time;
	
	//PD.write(std::cout << "det: ", det) << std::endl;
}
 
int main(int argc, char** argv)
{
	bool pass = true;
	
	size_t p = 7;
	size_t n = 20;
	size_t d = 5;
	size_t iterations = 2;
	size_t seed = time(NULL);
	std::string bumpsFilename;
	bool runLocal = false;
	
	static Argument args[] = {
		{ 'p', "-p P", "Set size of base field", TYPE_INT, &p },
		{ 'n', "-n N", "Set dim of test matrices to N.", TYPE_INT,  &n },
		{ 'd', "-d D", "Set degree of lumps", TYPE_INT, &d },
		{ 's', "-s S", "Random seed", TYPE_INT, &seed },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		{ 'b', "-b B", "Set the file containing bumps", TYPE_STR, &bumpsFilename },
		{ 'l', "-l L", "Run local smith form", TYPE_BOOL, &runLocal },
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);
	
	PolyRing PD(p);
	Util util(PD);
	
	std::vector<Polynomial> lumps;
	for (size_t i = 0; i < n * n; i++) {
		lumps.push_back(makeLump(PD, p, d));
	}
	
	std::vector<size_t> bump_idx;
	std::vector<Polynomial> bumps;
	
	{
		std::ifstream file(bumpsFilename);
		
		size_t nbumps;
		file >> nbumps;
		for (size_t i = 0; i < nbumps; i++) {			
			size_t idx;
			file >> idx;
			bump_idx.push_back(idx);
			
			Polynomial bump;
			PD.read(file, bump);
			bumps.push_back(bump);
		}
		
		file.close();
	}
		
	std::vector<Polynomial> diag;
	util.makeDiag(diag, bump_idx, bumps, n);
	
	Polynomial det;
	PD.assign(det, PD.one);
	for (size_t i = 0; i < diag.size(); i++) {
		// PD.write(std::cout, diag[i]) << std::endl;
		PD.mulin(det, diag[i]);
	}
	
	Matrix A(PD, n, n);
	util.makeExample(A, diag, lumps);
	
	// util.printMatrix(A);
	// PD.write(std::cout << "det: ", det) << std::endl << std::endl;
	
	std::cout << 2*d;
	solveTextBook(PD, A);
	solveKannanBachem(PD, A);
	solveIliopoulos(PD, A, det);
	// solveDet(PD, A);
	
	std::cout << std::endl;
	
	/*
	solveQuotTextBook(PD, A, det);
	if (runLocal) {
		solveLocal(PD, A, det);
	}
	// Compute inverse of matrix
	factorizeMatrix(PD, A);
	*/
	
	return pass ? 0 : -1;
}