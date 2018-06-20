#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

// Polynomial Matrix / Order Basis
#include "linbox/ring/modular.h"
#include "linbox/ring/givaro-poly.h"
#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"

//#define LINBOX_USES_OMP 1
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef GivaroPoly<Field> PolyRing;
typedef typename Field::Element Element;
typedef typename PolyRing::Element Polynomial;
typedef typename PolyRing::RandIter RandIter;

typedef MatrixDomain<PolyRing> MatDom;

typedef BlasMatrix<PolyRing> Matrix;

typedef SmithFormKannanBachemDomain<PolyRing> SmithDom;

void solveTextBook(PolyRing PD, SmithDom &SKB, Matrix &A) {
	typename MatDom::OwnMatrix ACopy(A);
	
	std::vector<Polynomial> lst;
	SKB.solveTextBook(lst, ACopy);
	
	std::cout << "Text Book Solve:" << std::endl;
	for (size_t i = 0; i < lst.size(); i++) {
		PD.write(std::cout << i << ": ", lst[i]) << std::endl;
	}
}

void solve(PolyRing PD, SmithDom &SKB, Matrix &A) {
	typename MatDom::OwnMatrix ACopy(A);
	
	std::vector<Polynomial> lst;
	SKB.solve(lst, ACopy);
	
	std::cout << "Solve:" << std::endl;
	for (size_t i = 0; i < lst.size(); i++) {
		PD.write(std::cout << i << ": ", lst[i]) << std::endl;
	}
}

void solveAdaptive(PolyRing PD, SmithDom &SKB, Matrix &A) {
	typename MatDom::OwnMatrix ACopy(A);
	
	std::vector<Polynomial> lst;
	SKB.solveAdaptive(lst, ACopy);
	
	std::cout << "Solve Adaptive:" << std::endl;
	for (size_t i = 0; i < lst.size(); i++) {
		PD.write(std::cout << i << ": ", lst[i]) << std::endl;
	}
}

void halfSolve(PolyRing PD, SmithDom &SKB, Matrix &A) {
	typename MatDom::OwnMatrix ACopy(A);
	
	std::vector<Polynomial> lst;
	SKB.halfSolve(lst, ACopy);
	
	std::cout << "Half Solve:" << std::endl;
	for (size_t i = 0; i < lst.size(); i++) {
		PD.write(std::cout << i << ": ", lst[i]) << std::endl;
	}
}

void makeTriangularMatrix(PolyRing PD, Matrix &M, bool upper, size_t p) {
	for (size_t i = 0; i < 3; i++) {
		M.setEntry(i, i, PD.one);
		
		for (size_t j = i+1; j < 3; j++) {
			std::vector<integer> pol; 
			for (size_t k = 0; k < 5; k++) {
				pol.push_back(rand() % p);
			}
			Polynomial tmp;
			PD.init(tmp, pol);
			
			if (upper) {
				M.setEntry(i, j, tmp);
			} else {
				M.setEntry(j, i, tmp);
			}
		}
	}
}

int main(int argc, char** argv) {
	size_t p = 7;
	Field F(p);
	GivaroPoly<Field> PD(F, "x");
	RandIter RI(PD);
	MatDom MD(PD);
	SmithDom SKB(PD);
	
	std::vector<int> v;
	Polynomial f, f2;
	PD.init(f, v = {1, 1, 0, 1});
	PD.mul(f2, f, f);
	
	Matrix A(PD, 3, 3);
	A.setEntry(0, 0, f2);
	A.setEntry(1, 1, PD.one);
	A.setEntry(2, 2, f);
	
	srand(time(NULL));
	
	Matrix L(PD, 3, 3);
	makeTriangularMatrix(PD, L, true, p);
	
	Matrix R1(PD, 3, 3);
	makeTriangularMatrix(PD, R1, false, p);
	
	Matrix R2(PD, 3, 3);
	makeTriangularMatrix(PD, R2, false, p);
	
	// A = L * R1 * A * R2
	MD.leftMulin(R1, A);
	MD.leftMulin(L, A);
	MD.mulin(A, R2);
	
	solveTextBook(PD, SKB, A);
	solve(PD, SKB, A);
	halfSolve(PD, SKB, A);
	solveAdaptive(PD, SKB, A);
	
	return 0;
}
