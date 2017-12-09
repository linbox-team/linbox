#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

#include "givaro/givtimer.h"

#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"
#include "linbox/ring/polynomial-local-x.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/block-coppersmith-domain.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-massey-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/algorithms/poly-smith-form-local-x.h"

#include "sparse-matrix-generator.h"
#include "test-poly-smith-form.h"

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

typedef Field::RandIter RandIter;
typedef MatrixDomain<Field> MatrixDom;
typedef typename MatrixDom::OwnMatrix Matrix;
typedef RandomDenseMatrix<RandIter, Field> RandomMatrix;
typedef BlackboxBlockContainer<Field, SparseMat> Sequence;
typedef BlockMasseyDomain<Field, Sequence> MasseyDom;
typedef BlockCoppersmithDomain<MatrixDom, Sequence> CoppersmithDom;

typedef NTL_zz_pX PolynomialRing;
typedef typename PolynomialRing::Element Polynomial;

typedef MatrixDomain<PolynomialRing> PolyMatrixDom;
typedef typename PolyMatrixDom::OwnMatrix PolyMatrix;
typedef SmithFormKannanBachemDomain<PolynomialRing> SmithFormDom;

typedef PolynomialLocalX<PolynomialRing> LocalRing;
typedef typename LocalRing::Element LocalPolynomial;
typedef MatrixDomain<LocalRing> LocalMatrixDom;
typedef typename LocalMatrixDom::OwnMatrix LocalMatrix;
typedef PolySmithFormLocalXDomain<PolynomialRing> LocalSmithFormDom;

Givaro::Timer TW;

void writeInvariantFactors(std::ostream &os, PolynomialRing &R, std::vector<Polynomial> &factors) {
	for (size_t i = 0; i < factors.size(); i++) {
		Polynomial f;
		R.monic(f, factors[i]);
		R.write(os, f) << std::endl;
	}
}

void computeDet(const PolynomialRing &R, Polynomial &det, std::vector<Polynomial> &factors) {
	R.assign(det, R.one);
	for (size_t i = 0; i < factors.size(); i++) {
		R.mulin(det, factors[i]);
	}
	
	Polynomial monic_det;
	R.monic(monic_det, det);
	R.assign(det, monic_det);
}

void timeTextbook(const PolynomialRing &R, std::vector<Polynomial> &result, const PolyMatrix &M) {
	SmithFormDom SFD(R);
	result.clear();
	PolyMatrix G(M);
	
	TW.clear();
	TW.start();
	
	SFD.solveTextBook(result, G);
	
	TW.stop();
	double sf_time = TW.usertime();
	std::cout << sf_time << " " << std::flush;
}

void timeKannanBachem(const PolynomialRing &R, std::vector<Polynomial> &result, const PolyMatrix &M) {
	SmithFormDom SFD(R);
	result.clear();
	PolyMatrix G(M);
	
	TW.clear();
	TW.start();
	
	SFD.solve(result, G);
	
	TW.stop();
	double sf_time = TW.usertime();
	std::cout << sf_time << " " << std::flush;
}

void timeHybrid(const PolynomialRing &R, std::vector<Polynomial> &result, const PolyMatrix &M) {
	SmithFormDom SFD(R);
	result.clear();
	PolyMatrix G(M);
	
	TW.clear();
	TW.start();
	
	SFD.solveAdaptive(result, G);
	
	TW.stop();
	double sf_time = TW.usertime();
	std::cout << sf_time << " " << std::flush;
}

void timeHybrid2(const PolynomialRing &R, std::vector<Polynomial> &result, const PolyMatrix &M) {
	SmithFormDom SFD(R);
	result.clear();
	PolyMatrix G(M);
	
	TW.clear();
	TW.start();
	
	SFD.solveAdaptive2(result, G);
	
	TW.stop();
	double sf_time = TW.usertime();
	std::cout << sf_time << " " << std::flush;
}

void timeIliopoulos(const PolynomialRing &R, std::vector<Polynomial> &result, const PolyMatrix &M, const Polynomial &det) {
	SmithFormDom SFD(R);
	result.clear();
	PolyMatrix G(M);
	
	TW.clear();
	TW.start();
	
	SFD.solveIliopoulos(result, G, det);
	
	TW.stop();
	double sf_time = TW.usertime();
	std::cout << sf_time << " " << std::flush;
}

void timeLocalX(const PolynomialRing &R, std::vector<Polynomial> &result, Polynomial &det, const PolyMatrix &M, size_t exponent) {
	LocalSmithFormDom SFD(R, exponent);
	LocalRing L(R, exponent);
	
	result.clear();
	LocalMatrix G(M, L);
	
	TW.clear();
	TW.start();
	
	SFD.solve(result, G);
	
	TW.stop();
	double sf_time = TW.usertime();
	std::cout << sf_time << " " << std::endl;
	
	R.assign(det, R.one);
	for (size_t i = 0; i < result.size(); i++) {
		R.mulin(det, result[i]);
	}
	
	det = NTL::trunc(det, exponent);
	NTL::MakeMonic(det);
}

int main(int argc, char** argv) {
	size_t p = 7;
	size_t n = 10;
	size_t b = 5;
	double sparsity = 0.05;
	int seed = time(NULL);
	size_t t = 2;
	
	std::string bumpFile;
	std::string matrixFile;
	std::string outFile;
	std::string algo = "hyb";

	static Argument args[] = {
		{ 'm', "-m M", "Name of file for bumps", TYPE_STR, &bumpFile},
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		{ 'n', "-n N", "Dimension of matrix", TYPE_INT, &n},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 's', "-s S", "Target sparsity of matrix", TYPE_DOUBLE, &sparsity},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
		{ 'a', "-a A", "Smith form algorithm to use", TYPE_STR, &algo},
		{ 't', "-t T", "Run iliopoulos with t-th largest invariant factor", TYPE_INT, &t},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	srand(seed);

	Field F(p);
	MatrixDomain<Field> MD(F);
	PolynomialRing R(p);
	SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
	TestPolySmithFormUtil<Field> util(F);
	
	SparseMat M(F);
	Polynomial det;
	
	if (matrixFile == "" && bumpFile == "") {
		std::cout << "Must provide either matrix or bumps input" << std::endl;
		return -1;
	} else if (matrixFile == "") {
		// create sparse matrix from bumps and compute determinant
		M.resize(n, n);
		Gen.generate(M, det, bumpFile, sparsity);
	} else {
		std::ifstream iF(matrixFile);
		M.read(iF);
		M.finalize();
		iF.close();
		
		assert(M.rowdim() == M.coldim());
		n = M.rowdim();
	}
	
	std::cout << n << " " << b << " " << Gen.sparsity(M) << " " << std::flush;
	
	// Generate random left and right projectors
	RandIter RI(F);
	RandomMatrix RM(F, RI);
	
	Matrix U(F, b, n);
	Matrix V(F, n, b);
	
	RM.random(U);
	RM.random(V);
	
	// Construct block sequence to input to BM
	Sequence seq(&M, F, U, V);
	
	// Compute minimal generating polynomial matrix
	MasseyDom BMD(&seq); // pascal
	CoppersmithDom BCD(MD, &seq, 10); // george
	
	std::vector<Matrix> minpoly;
	std::vector<size_t> degree;
	
	TW.clear();
	TW.start();
	
	//BMD.left_minpoly_rec(minpoly, degree);
	BCD.right_minpoly(minpoly);
	
	TW.stop();
	double bm_time = TW.usertime();
	std::cout << bm_time << " " << std::flush;
	
	TW.clear();
	TW.start();
	
	// Convert to matrix with polynomial entries
	PolyMatrix G(R, b, b);
	for (size_t i = 0; i < b; i++) {
		for (size_t j = 0; j < b; j++) {
			std::vector<long> coeffs;
			for (size_t k = 0; k < minpoly.size(); k++) {
				long coeff;
				F.convert(coeff, minpoly[k].getEntry(i, j));
				coeffs.push_back(coeff);
			}
			Polynomial tmp;
			R.init(tmp, coeffs);
			
			G.setEntry(i, j, tmp);
		}
	}
	
	TW.stop();
	double cv_time = TW.usertime();
	std::cout << cv_time << " " << std::flush;
	
	TestPolySmithFormUtil<PolynomialRing> putil(R);
	//putil.printMatrix(G);
	//std::cout << std::endl;
	
	// Compute smith form of generator
	std::vector<Polynomial> result;
	
	timeKannanBachem(R, result, G);
	timeHybrid(R, result, G);
	computeDet(R, det, result);
	
	std::vector<Polynomial> result2;
	timeIliopoulos(R, result2, G, det);
	timeIliopoulos(R, result2, G, result[result.size() - t]);
	
	Polynomial det2;
	timeLocalX(R, result2, det2, G, n + 1);
	
	//R.write(std::cout << "det: ", det) << std::endl;
	//R.write(std::cout << "det2: ", det2) << std::endl;
	std::cout << "Pass? " << (R.areEqual(det, det2) ? "True" : "False") << std::endl;
	std::cout << std::endl;
	
	if (outFile == "") {
		writeInvariantFactors(std::cout, R, result);
	} else {
		std::ofstream out(outFile);
		writeInvariantFactors(out, R, result);
		out.close();
	}
	
	return 0;
}


