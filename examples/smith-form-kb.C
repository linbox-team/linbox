#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

// Polynomial Matrix / Order Basis
#include "linbox/ring/modular.h"
#include "linbox/ring/givaro-poly.h"
#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem2.h"

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

typedef SmithFormKannanBachem2Domain<PolyRing> SmithDom;

void solveTextBook(PolyRing PD, SmithDom &SKB, Matrix &A) {
	typename MatDom::OwnMatrix ACopy(A);
	
	std::vector<Polynomial> lst;
	SKB.solveTextBook(lst, ACopy);
	
	for (size_t i = 0; i < lst.size(); i++) {
		PD.write(std::cout << i << ": ", lst[i]) << std::endl;
	}
}

void solve(PolyRing PD, SmithDom &SKB, Matrix &A) {
	typename MatDom::OwnMatrix ACopy(A);
	
	std::vector<Polynomial> lst;
	SKB.solve(lst, ACopy);
	
	for (size_t i = 0; i < lst.size(); i++) {
		PD.write(std::cout << i << ": ", lst[i]) << std::endl;
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
	
	Matrix M(PD, 3, 3);
	M.setEntry(0, 0, f2);
	M.setEntry(1, 1, PD.one);
	M.setEntry(2, 2, f);
	
	srand (time(NULL));
	
	Matrix L(PD, 3, 3);
	for (size_t i = 0; i < 3; i++) {
		L.setEntry(i, i, PD.one);
		
		for (size_t j = i+1; j < 3; j++) {
			std::vector<integer> pol; 
			for (size_t k = 0; k < 5; k++) {
				pol.push_back(rand() % p);
			}
			Polynomial tmp;
			PD.init(tmp, pol);
			L.setEntry(i, j, tmp);
		}
	}
	
	Matrix R(PD, 3, 3);
	for (size_t i = 0; i < 3; i++) {
		R.setEntry(i, i, PD.one);
		
		for (size_t j = i+1; j < 3; j++) {
			std::vector<integer> pol; 
			for (size_t k = 0; k < 5; k++) {
				pol.push_back(rand() % p);
			}
			Polynomial tmp;
			PD.init(tmp, pol);
			R.setEntry(j, i, tmp);
		}
	}
	
	for (size_t i = 0; i < 3; i++) {
		for (size_t j = 0; j < 3; j++) {
			Polynomial tmp;
			PD.write(std::cout, L.getEntry(tmp, i, j)) << std::endl;
		}
	}
	
	Matrix A(PD, 3, 3), LA(PD, 3, 3);
	MD.mul(LA, L, M);
	MD.mul(A, LA, R);
	
	solveTextBook(PD, SKB, A);
	solve(PD, SKB, A);
	
	return 0;
}
