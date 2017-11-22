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
typedef SparseMatrix<Field, SparseMatrixFormat::ELL> SparseMat;

typedef Field::RandIter RandIter;
typedef MatrixDomain<Field> MatrixDom;
typedef typename MatrixDom::OwnMatrix Matrix;
typedef RandomDenseMatrix<RandIter, Field> RandomMatrix;
typedef BlackboxBlockContainer<Field, SparseMat> Sequence;
typedef BlockMasseyDomain<Field, Sequence> MasseyDom;

typedef NTL_zz_pX PolynomialRing;
typedef typename PolynomialRing::Element Polynomial;

typedef MatrixDomain<PolynomialRing> PolyMatrixDom;
typedef typename PolyMatrixDom::OwnMatrix PolyMatrix;
typedef SmithFormKannanBachemDomain<PolynomialRing> SmithFormDom;

int main(int argc, char** argv)
{
	size_t p = 7;
	size_t n = 10;
	size_t b = 5;
	std::string bumpFile;
	double sparsity = 0.05;
	int seed = time(NULL);

	static Argument args[] = {
		{ 'n', "-n N", "Dimension of matrix", TYPE_INT, &n},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'm', "-m M", "Name of file for bumps", TYPE_STR, &bumpFile},
		{ 's', "-s S", "Target sparsity of matrix", TYPE_DOUBLE, &sparsity},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	srand(seed);

	Field F(p);
	PolynomialRing R(p);
	SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
	TestPolySmithFormUtil<Field> util(F);
	
	// create sparse matrix from bumps and compute determinant
	SparseMat M(F, n, n);
	Polynomial det;
	Gen.generate(M, det, bumpFile, sparsity);
	util.printMatrix(M);
	
	R.write(std::cout << "det: ", det) << std::endl;
	
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
	
	std::vector<Matrix> minpoly;
	std::vector<size_t> degree;
	
	Givaro::Timer TW;
	
	TW.clear();
	TW.start();
	
	BMD.left_minpoly_rec(minpoly, degree);
	
	TW.stop();
	double bm_time = TW.usertime();
	
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
	
	TestPolySmithFormUtil<PolynomialRing> putil(R);
	putil.printMatrix(G);
	std::cout << std::endl;
	
	// Compute smith form of generator
	SmithFormDom SFD(R);
	std::vector<Polynomial> result;
	
	TW.clear();
	TW.start();
	
	SFD.solve(result, G);
	// SFD.solveTextbook(result, G);
	// SFD.solveAdaptive(result, G); // half tb then ilio w/ computed det
	// SFD.solveIliopoulos(result, G, det);
	
	TW.stop();
	double sf_time = TW.usertime();
	
	for (size_t i = 0; i < result.size(); i++) {
		R.write(std::cout, result[i]) << std::endl;
	}
	
	std::cout << bm_time << " " << cv_time <<" " << sf_time << std::endl;

	return 0;
}


