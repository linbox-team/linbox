#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

#include "givaro/givtimer.h"

#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/block-coppersmith-domain.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-massey-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"
#include "linbox/algorithms/smith-form-local.h"

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

void writeInvariantFactors(std::ostream &os, PolynomialRing &R, std::vector<Polynomial> &factors) {
	for (size_t i = 0; i < factors.size(); i++) {
		Polynomial f;
		R.monic(f, factors[i]);
		R.write(os, f) << std::endl;
	}
}

int main(int argc, char** argv)
{
	size_t p = 7;
	size_t n = 10;
	size_t b = 5;
	double sparsity = 0.05;
	int seed = time(NULL);
	
	std::string bumpFile;
	std::string matrixFile;
	std::string outFile;

	static Argument args[] = {
		{ 'm', "-m M", "Name of file for bumps", TYPE_STR, &bumpFile},
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		{ 'n', "-n N", "Dimension of matrix", TYPE_INT, &n},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 's', "-s S", "Target sparsity of matrix", TYPE_DOUBLE, &sparsity},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
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
	MasseyDom BMD(&seq);
	CoppersmithDom BCD(MD, &seq, 10);
	
	std::vector<Matrix> minpoly;
	std::vector<size_t> degree;
	
	Givaro::Timer TW;
	
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
	SmithFormDom SFD(R);
	std::vector<Polynomial> result;
	
	TW.clear();
	TW.start();
	
	// SFD.solve(result, G);
	// SFD.solveTextbook(result, G);
	SFD.solveAdaptive(result, G); // half tb then ilio w/ computed det
	// SFD.solveIliopoulos(result, G, det);
	
	TW.stop();
	double sf_time = TW.usertime();
	std::cout << sf_time << " " << std::endl;
	
	if (outFile == "") {
		writeInvariantFactors(std::cout, R, result);
	} else {
		std::ofstream out(outFile);
		writeInvariantFactors(out, R, result);
		out.close();
	}
	
	return 0;
}


