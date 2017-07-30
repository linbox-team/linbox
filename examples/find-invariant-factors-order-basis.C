#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

// Polynomial Matrix / Order Basis
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"

//#define LINBOX_USES_OMP 1
#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

//#include "linbox/algorithms/coppersmith-invariant-factors.h"

// Computes the invariant factors of a sparse matrix (given in Matrix Market Format)
// Effectively times: TPL_omp, BlockCoppersmithDomain and KannanBachem

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::TPL> SparseMat;
//typedef SparseMatrix<Field, SparseMatrixFormat::TPL_omp> SparseMat;

typedef typename Field::RandIter RandIter;
typedef RandomDenseMatrix<RandIter, Field> RandomMatrix;

// Polynomial with matrix coefficients
typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> MatrixP;

// Matrix with polynomial coefficients
typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> PMatrix;

typedef MatrixP::Matrix Matrix;
typedef MatrixDomain<Field> MatDomain;

int main(int argc, char** argv)
{
	int earlyTerm = 10;
	int p = 97, b = 3;
	std::string mFname,oFname;

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 't', "-t T", "Early term threshold", TYPE_INT, &earlyTerm},
		{ 'b', "-b B", "Blocking factor", TYPE_INT, &b},
		{ 'm', "-m M", "Name of file for matrix M", TYPE_STR, &mFname},
		{ 'o', "-o O", "Name of file for output", TYPE_STR, &oFname},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	Field F(p);
	SparseMat M(F);

	{
		std::ifstream iF(mFname);
		M.read(iF);
		M.finalize();
		iF.close();
	}

	std::cout << "Finished reading" << std::endl;

	int n = M.coldim();
	MatrixP Series(F, b, b, 2*n);
	
	Matrix U(F, b, n);
	Matrix V(F, n, b);
	
	RandIter RI(F);
	RandomMatrix RDM(F, RI);
	RDM.random(U);
	RDM.random(V);
	
	MatDomain MD(F);
	for (int i = 0; i < 2*n; i++) {
		Matrix tmp(F, n, b);
		MD.copy(tmp, V);
		for (int j = 0; j < i; j++) {
			Matrix tmp2(F, n, b);
			MD.copy(tmp2, tmp);
			MD.mul(tmp, M, tmp2);
		}
		
		MD.mul(Series[i], U, tmp);
		MD.write(std::cout, Series[i]) << std::endl;
	}
	Series.write(std::cout) << std::endl;
	
	OrderBasis<Field> SB(F);
	
	MatrixP Sigma(F, n, n, (2*n)+1);
	std::vector<size_t> shift(n,0);
	
	SB.M_Basis(Sigma, Series, 2*n, shift);
	
	Sigma.write(std::cout) << std::endl;
	
	return 0;
}


