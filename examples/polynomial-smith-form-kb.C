#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

// Polynomial Matrix / Order Basis
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/smith-form-kannan-bachem.h"

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
	size_t p = 7;
	Field F(p);
	RandIter RI(F);
	PolynomialMatrixMulDomain<Field> PMD(F);
	PolynomialSmithFormKannanBachemDomain<Field> PKB(F);
	
	PMatrix A(F, 3, 3, 3);
	Field::Element tmp;
	
	F.assign(A.ref(0, 0, 0), F.one);
	F.assign(A.ref(1, 1, 1), F.one);
	F.assign(A.ref(2, 2, 0), F.zero);
	F.assign(A.ref(2, 2, 1), F.assign(tmp, 3));
	F.assign(A.ref(2, 2, 2), F.one);
	
	PMatrix L(F, 3, 3, 3);
	F.assign(L.ref(0, 0, 0), F.one);
	F.assign(L.ref(1, 1, 0), F.one);
	F.assign(L.ref(2, 2, 0), F.one);
	for (int i = 0; i < 3; i++) {
		for (int j = i+1; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				F.assign(L.ref(j, i, k), RI.random(tmp));
			}
		}
	}
	
	PMatrix R(F, 3, 3, 3);
	
	F.assign(R.ref(0, 0, 0), F.one);
	F.assign(R.ref(1, 1, 0), F.one);
	F.assign(R.ref(2, 2, 0), F.one);
	for (int i = 0; i < 3; i++) {
		for (int j = i+1; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				F.assign(R.ref(i, j, k), RI.random(tmp));
			}
		}
	}
	
	/*
	A.write(std::cout) << std::endl;
	L.write(std::cout) << std::endl;
	R.write(std::cout) << std::endl;
	std::cout << "--------------------------------" << std::endl;
	//*/
	
	PMatrix M(F, 3, 3, 1);
	PMatrix AR(F, 3, 3, 1);
	
	PMD.mul(AR, A, R);
	PMD.mul(M, L, AR);
	
	M.write(std::cout) << std::endl;
	
	PKB.solveTextbook(M);
	
	M.write(std::cout) << std::endl;
	
	return 0;
}
