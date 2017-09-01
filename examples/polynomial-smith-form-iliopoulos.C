#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

// Polynomial Matrix / Order Basis
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"
#include "linbox/algorithms/polynomial-matrix/smith-form-iliopoulos.h"

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
typedef PolynomialMatrix<PMType::matfirst,PMStorage::plain,Field> PMatrix;

// Matrix with polynomial coefficients
//typedef PolynomialMatrix<PMType::polfirst,PMStorage::plain,Field> PMatrix;

typedef GivaroPoly<Field> PolyDomain;
typedef typename GivaroPoly<Field>::Element Polynomial;

void writePolynomialToMatrix(PolyDomain &PD, PMatrix &M, size_t r, size_t c, const Polynomial &p) {
	Givaro::Degree d = PD.degree(p).value();
	
	for (size_t k = 0; k <= M.degree(); k++) {
		if (d >= k) {
			Element e;
			PD.getEntry(e, Givaro::Degree(k), p);
			PD.subdomain().assign(M.ref(r, c, k), e);
		} else {
			PD.subdomain().assign(M.ref(r, c, k), PD.subdomain().zero);
		}
	}
}

int main(int argc, char** argv) {
	size_t p = 7;
	Field F(p);
	PolyDomain PD(F, "x");
	RandIter RI(F);
	PolynomialMatrixMulDomain<Field> PMD(F);
	PolynomialSmithFormIliopoulosDomain<Field> PKB(F);
	
	std::vector<int> v;
	Polynomial f, f2;
	PD.init(f, v = {1, 1, 0, 1});
	PD.mul(f2, f, f);
	
	PD.write(std::cout << "f: ", f) << std::endl;
	PD.write(std::cout << "f^2: ", f2) << std::endl;
	
	PMatrix A(F, 3, 3, PD.degree(f2).value() + 1);
	writePolynomialToMatrix(PD, A, 1, 1, PD.one);
	writePolynomialToMatrix(PD, A, 2, 2, f);
	writePolynomialToMatrix(PD, A, 0, 0, f2);
	
	A.write(std::cout << "A: ") << std::endl;
	
	PMatrix L(F, 3, 3, 3);
	F.assign(L.ref(0, 0, 0), F.one);
	F.assign(L.ref(1, 1, 0), F.one);
	F.assign(L.ref(2, 2, 0), F.one);
	for (int i = 0; i < 3; i++) {
		for (int j = i+1; j < 3; j++) {
			for (int k = 0; k < 3; k++) {
				Field::Element tmp;
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
				Field::Element tmp;
				F.assign(R.ref(i, j, k), RI.random(tmp));
			}
		}
	}
	
	PMatrix M(F, 3, 3, 1);
	PMatrix AR(F, 3, 3, 1);
	
	PMD.mul(AR, A, R);
	PMD.mul(M, L, AR);
	
	M.setsize(M.real_degree() + 1);
	
	M.write(std::cout << "before:" << std::endl) << std::endl;
	
	Polynomial d;
	PD.mul(d, f, f);
	PD.mulin(d, f);
	PD.write(std::cout << "d: ", d) << std::endl;
		
	std::vector<Polynomial> result;
	PKB.solve(result, M, d);
	
	std::cout << "result: " << std::endl;
	for (size_t i = 0; i < result.size(); i++) {
		PD.write(std::cout << i << ": ", result[i]) << std::endl;
	}
	
	return 0;
}
