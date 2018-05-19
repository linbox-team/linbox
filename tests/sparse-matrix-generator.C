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
	size_t p = 3;
	size_t n = 10;
	size_t t = 0;
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
		{ 't', "-t T", "Trial matrix number", TYPE_INT, &t},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	srand(seed);

	Field F(p);
	PolynomialRing R(p);
	SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
	TestPolySmithFormUtil<Field> util(F);
	
	Polynomial det;
	std::vector<Polynomial> fs;
	
	SparseMat M(F);
	if (bumpFile == "" && t == 0) {
		std::vector<Polynomial> fs;
		Polynomial det;
	
        Gen.invariants(fs, n, fsnum);
        Gen.generate(M, det, fs, sparsity);
	} else if (bumpFile == "" && t != 0) {
		std::vector<integer> v;
		Polynomial x, xm1, xp1, x2, x3;
		R.init(x, v = {0, 1});
		R.assign(xm1, x);
		R.subin(xm1, R.one);
		R.assign(xp1, x);
		R.addin(xp1, R.one);
		R.mul(x2, x, x);
		R.mul(x3, x2, x);
		
		if (t == 1) {
			// trial 1 - identity
			Gen.addTriangle(fs, n, 1.0 / n, xm1);
		} else if (t == 2) {
			// trial 2 - 
			Gen.addTriangle(fs, n, 3.0 / n, xm1);
		} else if (t == 3) {
			// trial 3
			Gen.addTriangle(fs, n, 1, xm1);
		} else if (t == 4) {
			// trial 4
			Gen.addTriangle(fs, n, n / 3.0, xm1);
		} else if (t == 5) {
			// trial 5
			Gen.addTriangle(fs, n, n, xm1);
		} else if (t == 6) {
			// trial 6
			Gen.addTriangle(fs, n/2, 6.0 / n, xm1);
			Gen.addTriangle(fs, n/2, 1, xp1);
		} else if (t == 7) {
			// trial 7
			Gen.addTriangle(fs, n/2, 6.0 / n, x);
			Gen.addTriangle(fs, n/2, 1, xm1);
		} else if (t == 8) {
			// trial 8
			Gen.addTriangle(fs, n, 1.0 / n, x2);
		} else if (t == 9) {
			// trial 9
			Gen.addTriangle(fs, n, 1.0 / n, x3);
		} else if (t == 10) {
			// trial 10
			Gen.addTriangle(fs, n/2, 1.0 / n, x2);
			Gen.addTriangle(fs, n/2, 1, xm1);
		} else if (t == 11) {
			Gen.addTriangle(fs, n/2, n / 6.0, xm1);
			Gen.addTriangle(fs, n/2, 6.0 / n, x);
		} else if (t == 12) {
			Gen.addTriangle(fs, n/2, n / 6.0, xm1);
			Gen.addTriangle(fs, n/2, 1, xp1);
		} else if (t == 13) {
			Gen.addTriangle(fs, n/2, n / 6.0, xm1);
			Gen.addTriangle(fs, n/2, 0.5, xp1);
		}
		Gen.generate(M, det, fs, sparsity);
	} else if (bumpFile != "") {
		Polynomial det;
		
		// create sparse matrix from bumps and compute determinant
		M.resize(n, n);
		Gen.generate(M, det, bumpFile, sparsity);
	}
	
	if (mofile != "") {
		std::ofstream out(mofile);
		M.write(out);
		out.close();
		return 0;
	}
	
	return 0;
}
