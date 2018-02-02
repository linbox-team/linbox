#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

#include "givaro/givtimer.h"

#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"
#include "linbox/ring/polynomial-local-x.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/blackbox/scalar-matrix.h"

#include "linbox/algorithms/block-coppersmith-domain.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/blackbox-block-container-spmv.h"
#include "linbox/algorithms/blackbox-block-container-smmx.h"
#include "linbox/algorithms/block-massey-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/algorithms/poly-smith-form-local-x.h"
#include "linbox/algorithms/weak-popov-form.h"
#include "linbox/algorithms/poly-dixon.h"
#include "linbox/algorithms/invariant-factors.h"
#include "linbox/algorithms/wiedemann.h"

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

typedef BlackboxBlockContainer<Field, SparseMat> Sequence;
typedef BlockMasseyDomain<Field, Sequence> MasseyDom;
typedef BlockCoppersmithDomain<MatrixDom, Sequence> CoppersmithDom;

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
typedef PolyDixonDomain<PolynomialRing, QuotientRing> DixonDom;

Givaro::Timer TW;

class TestInvariantFactorsHelper {
public:
	const size_t _p;
	const Field F;
	const PolynomialRing R;
	const MatrixDomain<Field> MD;
	
	TestInvariantFactorsHelper(size_t p) : _p(p), F(p), R(p), MD(F) {};
	
	void writeInvariantFactors(std::ostream &os, std::vector<Polynomial> &factors) const {
		for (size_t i = 0; i < factors.size(); i++) {
			Polynomial f;
			R.monic(f, factors[i]);
			R.write(os, f) << std::endl;
		}
	}
	
	void computeMinpoly(
		std::vector<size_t> &degree,
		std::vector<Matrix> &minpoly,
		std::vector<size_t> &degree2,
		std::vector<Matrix> &minpoly2,
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
		Sequence seq(&M, F, U, V);
		
		// Compute minimal generating polynomial matrix
		// MasseyDom BMD(&seq); // pascal
		CoppersmithDom BCD(MD, &seq, 10); // george
		
		TW.clear();
		TW.start();
		
		//BMD.left_minpoly_rec(minpoly, degree);
		degree = BCD.right_minpoly(minpoly);
		
		TW.stop();
		double bm_time = TW.usertime();
		std::cout << bm_time << " " << std::flush;
		
		// Construct block sequence to input to BM
		typedef BlackboxBlockContainerSmmx<Field, SparseMat> FflasSequence;
		FflasSequence fseq(&M, F, U, V);
		
		// Compute minimal generating polynomial matrix
		// BlockMasseyDomain<Field, FflasSequence> BMD(&seq); // pascal
		BlockCoppersmithDomain<MatrixDom, FflasSequence> FBCD(MD, &fseq, 10); // george
		
		TW.clear();
		TW.start();
		
		//BMD.left_minpoly_rec(minpoly, degree);
		degree2 = FBCD.right_minpoly(minpoly2);
		
		TW.stop();
		double bm2_time = TW.usertime();
		std::cout << bm2_time << " " << std::flush;
	}
	
	double computeMinpolyFflas(std::vector<size_t> &degree, std::vector<Matrix> &minpoly, const SparseMat &M, size_t b) const {
		size_t n = M.rowdim();
		
		RandIter RI(F);
		RandomMatrix RM(F, RI);
		
		Matrix U(F, b, n);
		Matrix V(F, n, b);
		
		RM.random(U);
		RM.random(V);
		
		// Construct block sequence to input to BM
		typedef BlackboxBlockContainerSmmx<Field, SparseMat> FflasSequence;
		FflasSequence seq(&M, F, U, V);
		
		// Compute minimal generating polynomial matrix
		// BlockMasseyDomain<Field, FflasSequence> BMD(&seq); // pascal
		BlockCoppersmithDomain<MatrixDom, FflasSequence> BCD(MD, &seq, 10); // george
				
		TW.clear();
		TW.start();
		
		//BMD.left_minpoly_rec(minpoly, degree);
		degree = BCD.right_minpoly(minpoly);
		
		TW.stop();
		double bm_time = TW.usertime();
		std::cout << bm_time << " " << std::flush;
		
		return bm_time;
	}
	
	double convertMinPolyToPolyMatrix(PolyMatrix &G, const std::vector<Matrix> &minpoly) {
		size_t b = G.rowdim();
		
		TW.clear();
		TW.start();
		
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
		
		TW.stop();
		double cv_time = TW.usertime();
		// std::cout << cv_time << " " << std::flush;
		
		return cv_time;
	}
	
	void computeDet(Polynomial &det, std::vector<Polynomial> &factors) {
		R.assign(det, R.one);
		for (size_t i = 0; i < factors.size(); i++) {
			R.mulin(det, factors[i]);
		}
		
		Polynomial monic_det;
		R.monic(monic_det, det);
		R.assign(det, monic_det);
	}
	
	double timeKannanBachem(std::vector<Polynomial> &result, const PolyMatrix &M) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		bool mem_error = false;
		try {
			SFD.solve(result, G);
		} catch (const std::bad_alloc &e) {
			mem_error = true;
		}
		
		TW.stop();
		double sf_time = TW.usertime();
		
		if (!mem_error) {
			std::cout << sf_time << " " << std::flush;
		} else {
			std::cout << "MEM " << std::flush;
		}
		
		return mem_error ? -1 : sf_time;
	}
	
	void iliopoulos(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) {
		
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		SFD.solveIliopoulos(result, G, det);
	}
	
	double timeIliopoulos(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) {
	
		TW.clear();
		TW.start();
		
		iliopoulos(result, M, det);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << sf_time << " " << std::flush;
		
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
	
	double timeFactoredLocal(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) {	
	
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
			} else {
				localFactored(result, M, factors[i].first, factors[i].second);
			}
		}
		
		TW.stop();
		double fp_time = TW.usertime();
		std::cout << fp_time << " " << std::flush;
		
		return fp_time;
	}
	
	double timeFullyFactoredLocal(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) {	
	
		std::vector<std::pair<Polynomial, long>> factors;
		
		TW.clear();
		TW.start();
		
		R.factor(factors, det);
		
		result.clear();
		for (size_t i = 0; i < M.rowdim(); i++) {
			result.push_back(R.one);
		}
		
		for (size_t i = 0; i < factors.size(); i++) {
			if (factors[i].second == 1) {
				R.mulin(result[result.size() - 1], factors[i].first);
			} else {
				local(result, M, factors[i].first, factors[i].second);
			}
		}
		
		TW.stop();
		double fp_time = TW.usertime();
		std::cout << fp_time << " " << std::flush;
		
		return fp_time;
	}
	
	double timeFactoredIlio(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) {	
	
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
			} else if (R.isIrreducible(factors[i].first)) {
				local(result, M, factors[i].first, factors[i].second);
			} else {
				Polynomial modulus;
				R.pow(modulus, factors[i].first, factors[i].second);
				
				std::vector<Polynomial> part;
				iliopoulos(part, M, modulus);
				
				for (size_t i = 0; i < part.size(); i++) {
					if (R.isZero(part[i])) {
						R.mulin(result[i], modulus);
					} else {
						Polynomial tmp;
						R.gcd(tmp, part[i], modulus);
						
						R.mulin(result[i], tmp);
					}
				}
			}
		}
		
		TW.stop();
		double fp_time = TW.usertime();
		std::cout << fp_time << " " << std::flush;
		
		return fp_time;
	}
	
	size_t detLimit(const PolyMatrix &M, size_t dim) {
		size_t limit1 = 0;
		for (size_t i = 0; i < M.rowdim(); i++) {
			size_t max_degree = 0;
			for (size_t j = 0; j < M.coldim(); j++) {
				size_t deg = R.deg(M.getEntry(i, j));
				if (deg > max_degree) {
					max_degree = deg;
				}
			}
			
			limit1 += max_degree;
		}
		
		size_t limit2 = 0;
		for (size_t i = 0; i < M.coldim(); i++) {
			size_t max_degree = 0;
			for (size_t j = 0; j < M.rowdim(); j++) {
				size_t deg = R.deg(M.getEntry(j, i));
				if (deg > max_degree) {
					max_degree = deg;
				}
			}
			
			limit2 += max_degree;
		}
				
		return std::min(std::min(limit1, limit2), dim) + 1;
	}
	
	double timeLocalX(Polynomial &det, const PolyMatrix &M, size_t exponent) {
		LocalSmithFormDom SFD(R, exponent);
		LocalRing L(R, exponent);
				
		LocalMatrix G(M, L);
		
		TW.clear();
		TW.start();
		
		SFD.solveDet(det, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << sf_time << " " << std::flush;
		
		NTL::MakeMonic(det);
		
		return sf_time;
	}
	
	bool dixon(Polynomial &minpoly, const PolyMatrix &M, const Polynomial &f, size_t max_deg) const {
		SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
		QuotientRing QR(_p, f);
		
		DixonDom DD(R, QR);
		PolyMatrix Mx(M);
				
		PolyMatrix b(R, Mx.rowdim(), 1);
		for (size_t i = 0; i < b.rowdim(); i++) {
			Polynomial e;
			Gen.randomPolynomial(e, max_deg);
			
			b.setEntry(i, 0, e);
		}
		
		PolyMatrix x(R, Mx.rowdim(), 1);
		size_t m = 4 * ((max_deg / R.deg(f)) + 1);
		
		bool success = DD.solve(x, Mx, b, f, m);
		
		if (!success) {
			return false;
		}
		
		Polynomial fm;
		R.pow(fm, f, m);
		
		R.assign(minpoly, R.one);
		for (size_t i = 0; i < x.rowdim(); i++) {
			Polynomial numer, denom, tmp;
			DD.rat_recon(numer, denom, x.getEntry(i, 0), fm);
			
			R.lcm(tmp, minpoly, denom);
			R.assign(minpoly, tmp);
		}
		
		return true;
	}
	
	double timeDixon(Polynomial &minpoly, const PolyMatrix &M, size_t m) const {
		SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
		
		size_t max_deg = 0;
		for (size_t i = 0; i < M.rowdim(); i++) {
			for (size_t j = 0; j < M.coldim(); j++) {
				Polynomial tmp;
				M.getEntry(tmp, i, j);
				
				if (R.isZero(tmp)) {
					continue;
				}
				
				size_t deg = R.deg(tmp);
				max_deg = max_deg < deg ? deg : max_deg;
			}
		}
		
		TW.clear();
		TW.start();
		
		size_t d = 1;
		bool success = false;
		Polynomial f;
		for (size_t i = 0; i < 10 && !success; i++) {
			Gen.randomIrreducible(f, d++);
			
			success = dixon(minpoly, M, f, m);
		}
		
		TW.stop();
		double d_time = TW.usertime();
		std::cout << d_time << " " << std::flush;
		
		if (!success) {
			std::cout << "(failed) " << std::flush;
		} else {
			// R.write(std::cout, f) << " " << std::flush;
		}
		
		return d_time;
	}
	
	double timePopov(Polynomial &det, const PolyMatrix &M) {
		WeakPopovFormDom PFD(R);
		
		TW.clear();
		TW.start();
		
		PFD.solveDet(det, M);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << sf_time << " " << std::flush;
		
		NTL::MakeMonic(det);
		
		return sf_time;
	}
}; // End of TestInvariantFactorsHelper

int main(int argc, char** argv) {
	size_t p = 7;
	size_t n = 10;
	size_t b = 5;
	int seed = time(NULL);
	size_t t = 2;
	double probability = 0.5;
	size_t rank = 0;
	
	double sparsity = 0.05;
	size_t extend = 1;
	
	std::string bumpFile;
	std::string matrixFile;
	int fsnum = 1;
	std::string outFile;
	std::string mofile;

	static Argument args[] = {
		{ 'm', "-m M", "Name of file for bumps", TYPE_STR, &bumpFile},
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'F', "-F F", "Name of output file for matrix", TYPE_STR, &mofile},
		{ 'l', "-l L", "aka fsnum, Choose L-th divisor sequence (you also must set n)", TYPE_INT, &fsnum},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		{ 'n', "-n N", "Dimension of matrix", TYPE_INT, &n},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 's', "-s S", "Target sparsity of matrix", TYPE_DOUBLE, &sparsity},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
		{ 't', "-t T", "Target t-th largest invariant factor", TYPE_INT, &t},
		{ 'c', "-c C", "Choose b such that prob of t-th being correct > c", TYPE_DOUBLE, &probability},
		{ 'R', "-R R", "Random matrix with rank R", TYPE_INT, &rank},
		{ 'e', "-e E", "Extension field exponent (p^e)", TYPE_INT, &extend},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	srand(seed);

	Field F(p);
	PolynomialRing R(p);
	SparseMatrixGenerator<Field, PolynomialRing> Gen(F, R);
	TestPolySmithFormUtil<Field> util(F);
	
	SparseMat M(F);
	Polynomial det;
	
	TW.clear();
	TW.start();
	
	if (rank != 0) {
		M.resize(n, n);
		Gen.randomMatrix(M, n, rank, sparsity);
	} else if (matrixFile == "" && bumpFile == "") {
		std::vector<Polynomial> fs;
        Gen.invariants(fs, n, fsnum);
        Gen.generate(M, det, fs, sparsity);
//		std::cout << "Must provide either matrix or bumps input" << std::endl;
//		return -1;
	} else if (matrixFile == "") {
		// create sparse matrix from bumps and compute determinant
		M.resize(n, n);
		Gen.generate(M, det, bumpFile, sparsity);
	} else {
		std::ifstream iF(matrixFile);
		M.read(iF);
		M.finalize();
		iF.close();
	}
	
	TW.stop();
	//double mg_time = TW.usertime();
	// std::cout << mg_time << " " << std::flush;
		
	assert(M.rowdim() == M.coldim());
	n = M.rowdim();
	if (n <= 20) M.write(std::cout) << std::endl;
	
	if (mofile != "") {
		std::ofstream out(mofile);
		M.write(out);
		out.close();
		return 0;
	}
	
#if 0
	std::vector<Polynomial> lifs;
	InvariantFactors<Field, PolynomialRing, QuotientRing> IFD(F, R);
	IFD.largestInvariantFactors(lifs, M, t, probability);
	
	for (size_t i = 0; i < lifs.size(); i++) {
		R.write(std::cout, lifs[i]) << std::endl;
	}
	
	Element detM;
	IFD.det(detM, M, t, probability);
	F.write(std::cout, detM) << std::endl;
#else 
	TestInvariantFactorsHelper helper(p);
	std::cout << n << " " << Gen.nnz(M) << " " << b << " " << std::flush;
	
	if (b == 1) {
		if (extend > 1) {
			typedef Givaro::GFqDom<int64_t> ExtField;
			typedef typename ExtField::RandIter ExtRandIter;
			
			// uint64_t extend1 = (uint64_t)Givaro::FF_EXPONENT_MAX((uint64_t)p, (uint64_t)LINBOX_EXTENSION_DEGREE_MAX);
			uint64_t extend1 = extend;
			//std::cout << extend1 << " " << std::flush;
			
			ExtField EF((uint64_t) p, extend1);
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
			std::cout << bm_time << " " << (phi.size() - 1) << std::endl;
		} else {
			std::vector<size_t> degree;
			std::vector<Matrix> minpoly;
			std::vector<size_t> degree2;
			std::vector<Matrix> minpoly2;
			helper.computeMinpoly(degree, minpoly, degree2, minpoly2, M, b);
			
			std::cout << extend << " " << std::flush;
			
			RandIter RI(F);
			
			typedef BlackboxContainer<Field, SparseMat> BBContainer;
			BBContainer seq(&M, F, RI);
			MasseyDomain<Field, BBContainer> WD(&seq, 10);
						
			TW.clear();
			TW.start();
			
			BlasVector<Field> phi(F);
			unsigned long deg;
			WD.minpoly(phi, deg);
			
			TW.stop();
			double bm_time = TW.usertime();
			std::cout << bm_time << " " << (phi.size() - 1) << std::endl;
		}
		
		return 0;
	}
	
	Polynomial mp;
	std::vector<Polynomial> result;
	Polynomial det2;
	std::vector<Polynomial> result2;
	std::vector<Polynomial> result3;
	std::vector<Polynomial> result4;
	std::vector<Polynomial> result5;
	
	Polynomial detFflas;
	std::vector<Polynomial> resultFflas;
	
	// Generate random left and right projectors
	std::vector<size_t> degree;
	std::vector<Matrix> minpoly;
	std::vector<size_t> degree2;
	std::vector<Matrix> minpoly2;
	//helper.computeMinpoly(degree, minpoly, degree2, minpoly2, M, b);
	helper.computeMinpolyFflas(degree, minpoly, M, b);
	
	// Convert to matrix with polynomial entries
	PolyMatrix G(R, b, b);
	helper.convertMinPolyToPolyMatrix(G, minpoly);
	
	TestPolySmithFormUtil<PolynomialRing> putil(R);
	//putil.printMatrix(G);
	//std::cout << std::endl;
	
	// Compute smith form of generator
	size_t exponent_limit = helper.detLimit(G, n);
	std::cout << exponent_limit << " " << std::flush;
	
	helper.timeDixon(mp, G, exponent_limit);
	helper.timePopov(det, G);
	helper.timeLocalX(det2, G, exponent_limit);
	helper.timeFactoredLocal(result3, G, det2);
	helper.timeFactoredIlio(result4, G, det2);
	helper.timeFullyFactoredLocal(result5, G, det2);
	helper.timeIliopoulos(result2, G, det2);
	
	//Polynomial t1, t2;
	//R.monic(t1, mp);
	//R.monic(t2, result3[result3.size() - 1]);
	//std::cout << std::endl;
	//R.write(std::cout << "dixon = ", t1) << std::endl;
	//R.write(std::cout << "mp = ", t2) << std::endl;
	//std::string mpPass = (R.areEqual(t1, t2) ? "Pass" : "Fail");
	//std::cout << mpPass << " " << std::flush;
	
	//double kb_time = helper.timeKannanBachem(result, G);
	//timeHybrid(R, result, G);
	//helper.computeDet(det, result);
	
	//std::cout << (kb_time / total_time) << " ";
	//std::cout << (kb_time / total2_time) << " " << std::flush;
	
	//R.write(std::cout << "det1: ", det) << std::endl;
	//R.write(std::cout << "det2: ", det2) << std::endl;
	//std::cout << (R.areEqual(det, det2) ? "Pass" : "Fail");
	
	std::cout << std::endl;
		
	if (outFile == "") {
		//helper.writeInvariantFactors(std::cout, result3);
		//helper.writeInvariantFactors(std::cout, resultFflas);
	} else {
		std::ofstream out1(outFile + "1.txt");
		std::ofstream out2(outFile + "2.txt");
		std::ofstream out3(outFile + "3.txt");
		std::ofstream out4(outFile + "4.txt");
		
		helper.writeInvariantFactors(out1, result);
		helper.writeInvariantFactors(out2, result2);
		helper.writeInvariantFactors(out3, result3);
		helper.writeInvariantFactors(out4, result4);
		
		out1.close();
		out2.close();
		out3.close();
		out4.close();
	}
#endif
	
	return 0;
}


