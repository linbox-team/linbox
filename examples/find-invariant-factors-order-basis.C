#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

// Polynomial Matrix / Order Basis
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/order-basis.h"

// Smith Form Stuff
#include "linbox/ring/givaro-poly.h"
#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"

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

typedef GivaroPoly<Field> PolyDom;
typedef typename PolyDom::Element Polynomial;
typedef MatrixDomain<PolyDom> PolyMatDom;
typedef typename PolyMatDom::OwnMatrix GPMatrix;
typedef SmithFormKannanBachemDomain<PolyDom> SmithDom;

void printMatrix(Field F, const Matrix &A) {
	std::cout << "[" << std::endl;
	for (size_t i = 0; i < A.rowdim(); i++) {
		std::cout << "\t[";
		for (size_t j = 0; j < A.coldim(); j++) {
			F.write(std::cout, A.getEntry(i, j));
			if (j < A.coldim() - 1) {
				std::cout << ", ";
			}
		}
		std::cout << "]" << std::endl;
	}
	std::cout << "]" << std::endl;
}

void printGPMatrix(PolyDom PD, const GPMatrix &A) {
	std::cout << "[" << std::endl;
	for (size_t i = 0; i < A.rowdim(); i++) {
		std::cout << "\t[";
		for (size_t j = 0; j < A.coldim(); j++) {
			PD.write(std::cout, A.getEntry(i, j));
			if (j < A.coldim() - 1) {
				std::cout << ", ";
			}
		}
		std::cout << "]" << std::endl;
	}
	std::cout << "]" << std::endl;
}

template<class PMatrix>
void convertToBlasMatrix(const PolyDom &PD, GPMatrix &G, const PMatrix &M) {
	M.write(std::cout) << std::endl;
	
	for (size_t i = 0; i < M.rowdim(); i++) {
		for (size_t j = 0; j < M.coldim(); j++) {
			Polynomial tmp;
			PD.init(tmp, M(i, j));
			G.setEntry(i, j, tmp);
		}
	}
}

template<class PMatrix>
void writeToPolyMat(const Field F, PMatrix &PM, const Matrix &M, size_t k) {
	for (size_t i = 0; i < PM.rowdim(); i++) {
		for (size_t j = 0; j < PM.coldim(); j++) {
			F.assign(PM.ref(i, j, k), M.getEntry(i, j));
		}
	}
}

int main(int argc, char** argv)
{
	size_t p = 97, b = 3;
	std::string mFname,oFname;

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'b', "-b B", "Blocking factor", TYPE_INT, &b},
		{ 'm', "-m M", "Name of file for matrix M", TYPE_STR, &mFname},
		{ 'o', "-o O", "Name of file for output", TYPE_STR, &oFname},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	Field F(p);
	PolyDom PD(F, "x");
	PolyMatDom PMD(PD);
	
	SparseMat M(F);

	{
		std::ifstream iF(mFname);
		M.read(iF);
		M.finalize();
		iF.close();
	}

	std::cout << "Finished reading" << std::endl;

	size_t n = M.coldim();
	size_t d = 2*n;
	
	MatrixP Series(F, b, b, d);
	
	Matrix U(F, b, n);
	Matrix V(F, n, b);
	
	RandIter RI(F);
	RandomMatrix RDM(F, RI);
	RDM.random(U);
	RDM.random(V);
	
	MatDomain MD(F);
	Matrix tmp(F, n, b);
	MD.copy(tmp, V);
	
	Matrix UMiV(F, b, b);
	MD.mul(UMiV, U, tmp);
	writeToPolyMat(F, Series, UMiV, 0);
	
	for (size_t i = 1; i < d; i++) {
		Matrix tmp2(tmp);
		MD.mul(tmp, M, tmp2);
		
		MD.mul(UMiV, U, tmp);
		printMatrix(F, UMiV);
		
		writeToPolyMat(F, Series, UMiV, i);
		//MD.write(std::cout, Series[i]) << std::endl;
		Series.write(std::cout) << std::endl;
	}
	// Series.write(std::cout) << std::endl;
	
	GPMatrix gpSeries(PD, b, b);
	convertToBlasMatrix(PD, gpSeries, Series);
	std::cout << "Series Matrix: " << std::endl;
	printGPMatrix(PD, gpSeries);
	
	std::cout << std::endl << "d: " << d << std::endl;
	
	OrderBasis<Field> SB(F);
	
	MatrixP Sigma(F, b, b, d+1);
	std::vector<size_t> shift(b,0);
	
	SB.PM_Basis(Sigma, Series, d, shift);
	
	// Sigma.write(std::cout) << std::endl;
	
	GPMatrix A(PD, b, b);
	convertToBlasMatrix(PD, A, Sigma);
	printGPMatrix(PD, A);
	
	SmithDom SFD(PD);
	
	std::vector<Polynomial> lst;
	SFD.solveAdaptive(lst, A);
	
	for (size_t i = 0; i < b; i++) {
		std::cout << lst[i] << std::endl;
	}
	
	return 0;
}


