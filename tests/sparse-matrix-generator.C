#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "givaro/givtimer.h"

#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"

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

typedef NTL_zz_pX PolynomialRing;
typedef typename PolynomialRing::Element Polynomial;
typedef typename PolynomialRing::Coeff Coeff;

Givaro::Timer TW;

int main(int argc, char** argv) {
	size_t p = 7;
	size_t n = 10;
	int seed = time(NULL);
	
	double sparsity = 0.05;
	
	std::string bumpFile;
	std::string matrixFile;
	int fsnum = 1;
	std::string outFile;
	std::string mofile;

	static Argument args[] = {
		{ 'm', "-m M", "Name of file for bumps", TYPE_STR, &bumpFile},
		{ 'F', "-F F", "Name of output file for matrix", TYPE_STR, &mofile},
		{ 'l', "-l L", "aka fsnum, Choose L-th divisor sequence (you also must set n)", TYPE_INT, &fsnum},
		{ 'n', "-n N", "Dimension of matrix", TYPE_INT, &n},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 's', "-s S", "Target sparsity of matrix", TYPE_DOUBLE, &sparsity},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	srand(seed);

	Field F(p);
	PolynomialRing R(p);
	SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
	TestPolySmithFormUtil<Field> util(F);
	
	/*
	SparseMat M(F);
	if (bumpFile == "") {
		std::vector<Polynomial> fs;
		Polynomial det;
	
        Gen.invariants(fs, n, fsnum);
        Gen.generate(M, det, fs, sparsity);
	} else {
		Polynomial det;
		
		// create sparse matrix from bumps and compute determinant
		M.resize(n, n);
		Gen.generate(M, det, bumpFile, sparsity);
	}
	//*/
	
	std::vector<integer> v;
	Polynomial x, xm1;
	R.init(x, v = {0, 1});
	R.assign(xm1, x);
	R.subin(xm1, R.one);

	Polynomial det;
	std::vector<Polynomial> fs;
	Gen.addTriangle(fs, n/2, 6.0 / n, x);
	Gen.addTriangle(fs, n/2, 1, xm1);
	
	SparseMat M(F);
	Gen.generate(M, det, fs, sparsity);
	
	if (mofile != "") {
		std::ofstream out(mofile);
		M.write(out);
		out.close();
		return 0;
	}
	
	return 0;
}
