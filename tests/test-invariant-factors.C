#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "givaro/givtimer.h"

#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/wiedemann.h"

#include "linbox/algorithms/poly-smith-form.h"
#include "linbox/algorithms/invariant-factors.h"

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

typedef MatrixDomain<PolynomialRing> PolyMatrixDom;
typedef typename PolyMatrixDom::OwnMatrix PolyMatrix;

Givaro::Timer TW;

class TestInvariantFactorsHelper {
public:
	const size_t _p;
	const Field F;
	const PolynomialRing R;
	const MatrixDomain<Field> MD;
	
	TestInvariantFactorsHelper(size_t p) : _p(p), F(p), R(p), MD(F) {};
	
	template<class Matrix>
	double nnz(const Matrix &M) const {
		double nnz = 0;
		
		for (size_t i = 0; i < M.rowdim(); i++) {
			for (size_t j = 0; j < M.coldim(); j++) {
				if (F.isZero(M.getEntry(i,j))) {
					continue;
				}
				
				nnz++;
			}
		}
		
		return nnz;
	}

	void writeInvariantFactors(
		std::ostream &os,
		const std::vector<Polynomial> &factors) const {
	
		for (size_t i = 0; i < factors.size(); i++) {
			Polynomial f;
			R.monic(f, factors[i]);
			R.write(os, f) << std::endl;
		}
	}
	
	double extWiedemann(SparseMat &M, size_t extend) const {
		typedef Givaro::GFqDom<int64_t> ExtField;
		typedef typename ExtField::RandIter ExtRandIter;
		
		// uint64_t extend1 = (uint64_t)Givaro::FF_EXPONENT_MAX((uint64_t)p, (uint64_t)LINBOX_EXTENSION_DEGREE_MAX);
		uint64_t extend1 = extend;
		//std::cout << extend1 << "\t" << std::flush;
		
		ExtField EF((uint64_t) _p, extend1);
		ExtRandIter RI(EF);
		
		typedef typename SparseMat::template rebind<ExtField>::other FBlackbox;
		FBlackbox EM(M, EF);
		
		typedef BlackboxContainer<ExtField, FBlackbox> BBContainer;
		BBContainer seq(&EM, EF, RI);
		MasseyDomain<ExtField, BBContainer> WD(&seq, 10);
		
		TW.clear();
		TW.start();
		
		BlasVector<ExtField> phi(EF);
		unsigned long deg;
		WD.minpoly(phi, deg);
		
		TW.stop();
		double bm_time = TW.usertime();
		std::cout << bm_time << "\t" << (phi.size() - 1) << std::endl;
		
		return bm_time;
	}
}; // End of TestInvariantFactorsHelper

void time1(auto fun) {
	TW.clear();
	TW.start();
	
	bool success = fun();
	
	TW.stop();
	std::cout << TW.usertime();
	if (!success) {
		std::cout << "(failed)";
	}
	std::cout << "\t" << std::flush;
}

int main(int argc, char** argv) {
	size_t p = 7;
	size_t extend = 1;
	size_t b = 5;
	
	std::string matrixFile;
	std::string outFile;
	
	int seed = time(NULL);

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'e', "-e E", "Extension field exponent (p^e)", TYPE_INT, &extend},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
		
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	srand(seed);

	Field F(p);
	PolynomialRing R(p);
	
	SparseMat M(F);
	{
		std::ifstream iF(matrixFile);
		M.read(iF);
		M.finalize();
		iF.close();
	}
		
	assert(M.rowdim() == M.coldim());
	size_t n = M.rowdim();
	if (n <= 20) M.write(std::cout) << std::endl;
	
	TestInvariantFactorsHelper helper(p);
	std::cout << n << "\t" << helper.nnz(M) << "\t" << b << "\t" << std::flush;
	
	if (b == 1) {
		if (extend > 1) {
			helper.extWiedemann(M, extend);
		} else {
			// TODO: implement
		}
		
		return 0;
	}
	
	PolySmithFormDomain<PolynomialRing> PSFD(R);
	InvariantFactors<Field, PolynomialRing> IFD(F, R);
		
	// Generate random left and right projectors
	std::vector<size_t> degree;
	std::vector<Matrix> minpoly;
	time1([&](){IFD.computeGenerator(minpoly, M, b); return true;});
	
	// Convert to matrix with polynomial entries
	PolyMatrix G(R, b, b);
	IFD.convert(G, minpoly);
	
	// Compute smith form of generator
	size_t exponent_limit = PSFD.detLimit(G);
	std::cout << exponent_limit << "\t" << std::flush;
	
	Polynomial det, mp;
	std::vector<Polynomial> result;
	time1([&](){return PSFD.dixon(mp, G);});
	time1([&](){PSFD.detPopov(det, G); return true;});
	time1([&](){PSFD.detLocalX(det, G); return true;});
	time1([&](){PSFD.solve(result, G, det); return true;});
	std::cout << std::endl;
	
	if (outFile != "") {
		std::ofstream out(outFile + ".txt");
		helper.writeInvariantFactors(out, result);
		out.close();
	}
	
	return 0;
}
