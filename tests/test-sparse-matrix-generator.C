#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "sparse-matrix-generator.h"
#include "test-poly-smith-form.h"

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

typedef NTL_zz_pX PolynomialRing;
typedef typename PolynomialRing::Element Polynomial;

int main(int argc, char** argv)
{
	int p = 7;
	int n = 10;
	double sparsity = 0.5;
	std::string bumpFile;

	static Argument args[] = {
		{ 'n', "-n N", "Dimension of matrix", TYPE_INT, &n},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'b', "-b B", "Name of file for bumps", TYPE_STR, &bumpFile},
		{ 's', "-s S", "Target sparsity", TYPE_DOUBLE, &sparsity},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	srand(time(NULL));

	Field F(p);
	PolynomialRing R(p);
	SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
	
	std::vector<Polynomial> fs;
	Gen.readFile(fs, bumpFile);
	
	std::cout << "Finished reading" << std::endl;

	for (size_t i = 0; i < fs.size(); i++) {
		R.write(std::cout, fs[i]) << std::endl;
	}
	
	SparseMat M(F, n, n);
	
	Gen.build(M, fs);
	
	TestPolySmithFormUtil<Field> util(F);
	util.printMatrix(M);
	
	Gen.fillIn(M, sparsity);
	util.printMatrix(M);
	
	std::cout << "sparsity: " << Gen.sparsity(M) << std::endl;

	return 0;
}


