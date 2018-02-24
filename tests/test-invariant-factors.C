#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "givaro/givtimer.h"

#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"
#include "linbox/ring/polynomial-local-x.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/block-coppersmith-domain.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/blackbox-block-container-spmv.h"
#include "linbox/algorithms/blackbox-block-container-smmx.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/algorithms/poly-smith-form-local-x.h"
#include "linbox/algorithms/weak-popov-form.h"
#include "linbox/algorithms/invariant-factors.h"
#include "linbox/algorithms/wiedemann.h"

#include "linbox/algorithms/poly-smith-form.h"

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

typedef MatrixDomain<PolynomialRing> PolyMatrixDom;
typedef typename PolyMatrixDom::OwnMatrix PolyMatrix;
typedef SmithFormKannanBachemDomain<PolynomialRing> SmithFormDom;

typedef NTL_zz_pE QuotientRing;
typedef typename QuotientRing::Element QPolynomial;
typedef MatrixDomain<QuotientRing> QuotMatrixDom;
typedef typename QuotMatrixDom::OwnMatrix QuotMatrix;
typedef SmithFormLocal<QuotientRing> SmithFormLocalDom;

typedef PolynomialLocalX<PolynomialRing> LocalRing;
typedef typename LocalRing::Element LocalPolynomial;
typedef MatrixDomain<LocalRing> LocalMatrixDom;
typedef typename LocalMatrixDom::OwnMatrix LocalMatrix;
typedef PolySmithFormLocalXDomain<PolynomialRing> LocalSmithFormDom;
typedef WeakPopovFormDomain<PolynomialRing> WeakPopovFormDom;

Givaro::Timer TW;

class TestInvariantFactorsHelper {
public:
	const size_t _p;
	const Field F;
	const PolynomialRing R;
	const MatrixDomain<Field> MD;
	
	TestInvariantFactorsHelper(size_t p) : _p(p), F(p), R(p), MD(F) {};
	
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
	
	double computeMinpoly(
		std::vector<size_t> &degree,
		std::vector<Matrix> &minpoly,
		const SparseMat &M, 
		size_t b) const {
	
		size_t n = M.rowdim();
		
		RandIter RI(F);
		RandomMatrix RM(F, RI);
		
		Matrix U(F, b, n);
		Matrix V(F, n, b);
		
		RM.random(U);
		RM.random(V);
		
		// Construct block sequence to input to BM
		typedef BlackboxBlockContainerSmmx<Field, SparseMat> Sequence;
		Sequence seq(&M, F, U, V);
		BlockCoppersmithDomain<MatrixDom, Sequence> BCD(MD, &seq, 10);
				
		TW.clear();
		TW.start();
		
		degree = BCD.right_minpoly(minpoly);
		
		TW.stop();
		double bm_time = TW.usertime();
		std::cout << bm_time << "\t" << std::flush;
		
		return bm_time;
	}
	
	void convertMinPolyToPolyMatrix(
		PolyMatrix &G, 
		const std::vector<Matrix> &minpoly) const {
	
		size_t b = G.rowdim();
		for (size_t i = 0; i < b; i++) {
			for (size_t j = 0; j < b; j++) {
				std::vector<long> coeffs;
				for (size_t k = 0; k < minpoly.size(); k++) {
					long coeff;
					F.convert(coeff, minpoly[k].getEntry(i, j));
					coeffs.push_back(coeff);
				}
				Polynomial tmp;
				R.init(tmp, coeffs);
				
				G.setEntry(i, j, tmp);
			}
		}
	}
	
	double timeIliopoulos(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) const {
		
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
	
		TW.clear();
		TW.start();
		
		SFD.solveIliopoulos(result, G, det);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << sf_time << "\t" << std::flush;
		
		return sf_time;
	}
	
	void local(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &f,
		long multiplicity) const {
	
		SmithFormLocalDom SFD;
		
		Polynomial modulus;
		R.pow(modulus, f, multiplicity);
		
		QuotientRing QR(_p, modulus);
		
		QuotMatrix QM(M, QR);
		
		std::list<QPolynomial> L;
		SFD(L, QM, QR);
		
		Hom<PolynomialRing, QuotientRing> hom(R, QR);
		
		size_t j = 0;
		std::list<QPolynomial>::const_iterator it;
		for (it = L.begin(); it != L.end(); it++) {
			Polynomial tmp;
			hom.preimage(tmp, *it);
								
			if (R.isOne(tmp)) {
				// noop
			} else if (R.isZero(tmp)) {
				R.mulin(result[j], modulus);
			} else {
				R.mulin(result[j], tmp);
			}
			j++;
		}
	}
	
	void localRank(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &f) const {
	
		SmithFormLocalDom SFD;
		
		QuotientRing QR(_p, f);
		
		QuotMatrix QM(M, QR);
		
		std::list<QPolynomial> L;
		SFD(L, QM, QR);
		
		std::list<QPolynomial>::reverse_iterator it = L.rbegin();
		it++;
		if (QR.isZero(*it)) {
			R.mulin(result[result.size() - 2], f);
			R.mulin(result[result.size() - 1], f);
		} else {
			Polynomial f2;
			R.pow(f2, f, 2);
			R.mulin(result[result.size() - 1], f2);
		}
	}
	
	void localFactored(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &sf_factor,
		long multiplicity) const {
	
		std::vector<std::pair<Polynomial, long>> factors;
		R.factor(factors, sf_factor);
		
		for (size_t i = 0; i < factors.size(); i++) {
			local(result, M, factors[i].first, factors[i].second * multiplicity);
		}
	}
	
	void localFactoredRank(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &sf_factor) const {
	
		std::vector<std::pair<Polynomial, long>> factors;
		R.factor(factors, sf_factor);
		
		for (size_t i = 0; i < factors.size(); i++) {
			localRank(result, M, factors[i].first);
		}
	}
	
	double timeFactoredLocal(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) const {
		
		std::vector<std::pair<Polynomial, long>> factors;
		
		TW.clear();
		TW.start();
		
		R.squareFree(factors, det);
		
		result.clear();
		for (size_t i = 0; i < M.rowdim(); i++) {
			result.push_back(R.one);
		}
		
		for (size_t i = 0; i < factors.size(); i++) {
			if (factors[i].second == 1) {
				R.mulin(result[result.size() - 1], factors[i].first);
			} else if (factors[i].second == 2) {
				localFactoredRank(result, M, factors[i].first);
			} else {
				localFactored(result, M, factors[i].first, factors[i].second);
			}
		}
		
		TW.stop();
		double fp_time = TW.usertime();
		std::cout << fp_time << "\t" << std::flush;
		
		return fp_time;
	}
	
	template<class Iterator>
	size_t limit(const Iterator &begin, const Iterator &end) const {
		auto comp = [&](auto &a, auto &b){return R.deg(a) < R.deg(b);};
		auto acc = [&](size_t v, auto it) {
			return v + R.deg(*std::max_element(it.begin(), it.end(), comp));
		};
		return std::accumulate(begin, end, 0, acc);
	}
	
	size_t detLimit(const PolyMatrix &M, size_t dim) const {
		size_t limit1 = limit(M.rowBegin(), M.rowEnd());
		size_t limit2 = limit(M.colBegin(), M.colEnd());
		
		return std::min(std::min(limit1, limit2), dim) + 1;
	}
	
	double timeLocalX(
		Polynomial &det, 
		const PolyMatrix &M, 
		size_t exponent) const {
	
		LocalSmithFormDom SFD(R, exponent);
		LocalRing L(R, exponent);
				
		LocalMatrix G(M, L);
		
		TW.clear();
		TW.start();
		
		SFD.solveDet(det, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << sf_time << "\t" << std::flush;
		
		NTL::MakeMonic(det);
		
		return sf_time;
	}
	
	double timePopov(Polynomial &det, const PolyMatrix &M) const {
		WeakPopovFormDom PFD(R);
		
		TW.clear();
		TW.start();
		
		PFD.solveDet(det, M);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << sf_time << "\t" << std::flush;
		
		NTL::MakeMonic(det);
		
		return sf_time;
	}
	
	double timeDixon(Polynomial &mp, const PolyMatrix &M, size_t m) const {
		PolySmithFormDomain<PolynmialRing> PSFD(R);
		
		TW.clear();
		TW.start();
		
		PSFD.dixon(mp, M, m);
		
		TW.stop();
		double ds_time = TW.usertime();
		std::cout << ds_time << "\t" << std::flush;
		
		NTL::MakeMonic(mp);
		
		return ds_time;
	}
}; // End of TestInvariantFactorsHelper

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
	SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
	TestPolySmithFormUtil<Field> util(F);
	
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
	std::cout << n << "\t" << Gen.nnz(M) << "\t" << b << "\t" << std::flush;
	
	if (b == 1) {
		if (extend > 1) {
			helper.extWiedemann(M, extend);
		} else {
			// TODO: implement
		}
		
		return 0;
	}
	
	PolySmithFormDomain<PolynomialRing> PSFD(R);
		
	// Generate random left and right projectors
	std::vector<size_t> degree;
	std::vector<Matrix> minpoly;
	helper.computeMinpoly(degree, minpoly, M, b);
	
	// Convert to matrix with polynomial entries
	PolyMatrix G(R, b, b);
	helper.convertMinPolyToPolyMatrix(G, minpoly);
	
	// Compute smith form of generator
	size_t exponent_limit = helper.detLimit(G, n);
	std::cout << exponent_limit << "\t" << std::flush;
	
	Polynomial det;
	std::vector<Polynomial> result;
	//helper.timeLocalX(det2, G, exponent_limit);
	//helper.timeFactoredLocal(result3, G, det2);
	std::cout << std::endl;
		
	if (outFile != "") {
		std::ofstream out(outFile + ".txt");
		helper.writeInvariantFactors(out, result);
		out.close();
	}
	
	return 0;
}
