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
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/permutation.h"
#include "linbox/blackbox/scalar-matrix.h"

#include "linbox/algorithms/block-coppersmith-domain.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-massey-domain.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"
#include "linbox/algorithms/smith-form-local.h"
#include "linbox/algorithms/poly-smith-form-local-x.h"
#include "linbox/algorithms/invariant-factors.h"
#include "linbox/solutions/det.h"

#include "sparse-matrix-generator.h"
#include "test-poly-smith-form.h"
#include "linbox/solutions/det.h"

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

typedef x ExtensionField;
typedef SparseMatrix<ExtensionField, SparseMatrixFormat::CSR> SparseMatExtension;

typedef Field::RandIter RandIter;
typedef MatrixDomain<Field> MatrixDom;
typedef typename MatrixDom::OwnMatrix Matrix;
typedef RandomDenseMatrix<RandIter, Field> RandomMatrix;
typedef DenseMatrix<Field> DM;
//#define PRECONDITION

typedef SparseMat Preconditioner;
//typedef Permutation<Field> Preconditioner; // ..there will be others

typedef Compose<SparseMat,Preconditioner> ProductMat; // yes preconditioner

typedef BlackboxBlockContainer<Field, ProductMat> Sequence;
typedef BlockMasseyDomain<Field, Sequence> MasseyDom;
typedef BlockCoppersmithDomain<MatrixDom, Sequence> CoppersmithDom;

typedef NTL_zz_pX PolynomialRing;
typedef typename PolynomialRing::Element Polynomial;

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

Givaro::Timer TW;

class TestInvariantFactorsHelper {
public:
	const size_t _p;
	const Field F;
	const PolynomialRing R;
	const MatrixDomain<Field> MD;
	
	TestInvariantFactorsHelper(size_t p) : _p(p), F(p), R(p), MD(F) {};
	
	void writeInvariantFactors(std::ostream &os, std::vector<Polynomial> &factors, bool full=true) const {
		for (size_t i = 0; i < factors.size(); i++) {
			Polynomial f;
			R.monic(f, factors[i]);
			if (full)
				R.write(os, f) << std::endl;
			else
				os << R.deg(f) << std::endl;
		}
	}

	const PolynomialRing& ring() const { return R; }

	// count nontrivial (not 1 or x) invariant factors.
	size_t countInvariantFactors(std::vector<Polynomial> &factors) const {
		Polynomial x; 
		std::vector<long> coeffs; coeffs.push_back(0); coeffs.push_back(1);
		R.init(x, coeffs);

		size_t count = 0;
		for (size_t i = 0; i < factors.size(); i++) {
			if (not R.areEqual(factors[i], R.one) and not R.areEqual(factors[i], x))
				count++;
		}
		return count;
	}
	
	// total degree of invariant factors.
	size_t degInvariantFactors(std::vector<Polynomial> &factors) const {
		size_t d = 0;
		for (size_t i = 0; i < factors.size(); i++) 
			d += R.deg(factors[i]);
		return d;
	}
	
	double computeMinpoly(std::vector<size_t> &degree, std::vector<Matrix> &minpoly, const ProductMat &M, size_t b) const {
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
		MasseyDom BMD(&seq); // pascal
		CoppersmithDom BCD(MD, &seq, 10); // george
				
		TW.clear();
		TW.start();
		
		//BMD.left_minpoly_rec(minpoly, degree);
		degree = BCD.right_minpoly(minpoly);
		
		TW.stop();
		double bm_time = TW.usertime();
		std::cout << "seq->gen " << bm_time << " " << std::flush;
		
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
	
	void timeTextbook(std::vector<Polynomial> &result, const PolyMatrix &M) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		SFD.solveTextBook(result, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << ", text_time " << sf_time << " " << std::flush;
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
			std::cout << ", kb_time " << sf_time << " " << std::flush;
		} else {
			std::cout << "MEM " << std::flush;
		}
		
		return mem_error ? -1 : sf_time;
	}
	
	void timeHybrid(std::vector<Polynomial> &result, const PolyMatrix &M) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		SFD.solveAdaptive(result, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << ", hybrid " << sf_time << " " << std::flush;
	}
	
	void timeHybrid2(std::vector<Polynomial> &result, const PolyMatrix &M) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		SFD.solveAdaptive2(result, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		std::cout << ", hybrid2 " << sf_time << " " << std::flush;
	}
	
	double timeIliopoulos(std::vector<Polynomial> &result, const PolyMatrix &M, const Polynomial &det) {
		SmithFormDom SFD(R);
		result.clear();
		PolyMatrix G(M);
		
		TW.clear();
		TW.start();
		
		SFD.solveIliopoulos(result, G, det);
		
		TW.stop();
		double sf_time = TW.usertime();
		//std::cout << sf_time << " " << std::flush;
		
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
		
		std::vector<Polynomial> result;
		result.clear();
		
		LocalMatrix G(M, L);
		
		TW.clear();
		TW.start();
		
		SFD.solveDet(det, G);
		
		TW.stop();
		double sf_time = TW.usertime();
		//std::cout << ", local " << sf_time << " " << std::flush;
		
		NTL::MakeMonic(det);
		
		return sf_time;
	}
}; // End of TestInvariantFactorsHelper

// preconditioners

void identity(SparseMat & P)
{
	for (size_t i=0; i < P.rowdim(); ++i)
		P.setEntry(i,i,P.field().one);
	P.finalize();
}

void randomVec(std::vector<size_t> & V, size_t n) 
{
	V.resize(n);
	for (size_t i = 0; i < n; ++i) V[i] = i;
	// Knuth construction
	for (size_t i = 0; i < n-1; ++i) {
		size_t j = i + rand()%(n-i);
		std::swap(V[i], V[j]);
	}
}

void randomCycle(SparseMat & P)
{
	size_t n = P.rowdim();
	std::vector<size_t> Perm;
	randomVec(Perm, n);
//	for (size_t i = 0; i < Perm.size(); ++i) std::cout << Perm[i] << " ";
//	std::cout << std::endl;
	size_t first = Perm.back(); Perm.pop_back();
	size_t i = first;
	while (Perm.size() > 0) {
		size_t j = Perm.back(); Perm.pop_back();
		P.setEntry(i,j, P.field().one);
		i = j;
	}
	P.setEntry(i,first, P.field().one);
	P.finalize();
} 

void randomPerm(SparseMat & P)
{
	size_t n = P.rowdim();
	std::vector<size_t> Perm;
	randomVec(Perm, n);
	for (size_t i = 0; i < n; ++i) P.setEntry(i, Perm[i], P.field().one);
	P.finalize();
}

void fillBlock(SparseMat & P, size_t s, size_t b)
{
	const Field& F = P.field();
	BlasMatrix<Field> B(F,b,b);
	Field::Element x; F.init(x);
//	SparseMat B(P.field(),b,b);

	Field::Element d; F.init(d);
	do { // hack because BlasMatrix::random() is broken 
		for (size_t i = 0; i < b; ++i)
			for (size_t j = 0; j < b; ++j) {
				F.init(x,rand());
				B.setEntry(i,j,x);
			}
	} while (F.isZero(det(d,B))); 

	for (size_t i = 0; i < b; ++i)
		for (size_t j = 0; j < b; ++j)
			P.setEntry(i+s,j+s,B.getEntry(x,i,j));
}

void blockDiag(SparseMat & P, size_t b=5)
{
	BlasMatrix<Field> B(P.field(), b, b);
	size_t n = P.rowdim(); size_t q = n/b; size_t r = n - q*b;
	for (size_t k = 0; k < q-1; ++k)
		fillBlock(P,b*k,b);
	fillBlock(P,b*(q-1),b+r);
	P.finalize();
}

int main(int argc, char** argv) {
	size_t p = 3;
	size_t n = 10;
	size_t t = 2;
	double probability = .99;
	double sparsity = 0.05;
	int seed = time(NULL);
	size_t times = 1;
	size_t exponent = 0;
	int pre = 0;

	size_t b = 0;
	
	std::string bumpFile;
	std::string matrixFile;
	int fsnum = 1;
	std::string outFile;

	static Argument args[] = {
		{ 'v', "-v V", "Version of preconditioner", TYPE_INT, &pre},
		{ 'm', "-m M", "Name of file for bumps", TYPE_STR, &bumpFile},
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'l', "-l L", "aka fsnum, Choose L-th divisor sequence (you also must set n)", TYPE_INT, &fsnum},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		{ 'n', "-n N", "Dimension of matrix", TYPE_INT, &n},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 's', "-s S", "Target sparsity of matrix", TYPE_DOUBLE, &sparsity},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
		{ 't', "-t T", "Run iliopoulos with t-th largest invariant factor", TYPE_INT, &t},
		{ 'k', "-k K", "Repeat computation K times", TYPE_INT, &times},
		{ 'e', "-e E", "Compute local at x^e", TYPE_INT, &exponent},
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
	
	if (matrixFile == "" && bumpFile == "") {
		std::vector<Polynomial> fs;
		Gen.invariants(fs, n, fsnum);
		Gen.generate(M, det, fs, sparsity);
		//std::cout << "Must provide either matrix or bumps input" << std::endl;
		//return -1;
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
		
	assert(M.rowdim() == M.coldim());
	n = M.rowdim();

	// Permutation preconditioner //
	Preconditioner P(F,n,n); 
	switch (pre) {
		case 0: identity(P); break;
		case 1: randomPerm(P); break;
		case 2: randomCycle(P); break;
		case 3: blockDiag(P,5); break;
	}

	//P.write(std::cout) << std::endl;

	ProductMat MP(M,P);
	/*
	if (n <= 20) {
		DM U(F,n,n), V(F,n,n);
		for (size_t i = 0; i < n; ++i) V.setEntry(i,i,F.one);
		MP.applyRight(U,V);
		U.write(std::cout) << std::endl;
	}
	*/

	std::vector<Polynomial> result;
	
	InvariantFactors<Field, PolynomialRing, QuotientRing> IFD(F,R);
	b = IFD.min_block_size(t,probability);
	TestInvariantFactorsHelper helper(p);

	for (size_t i = 0; i < times; i++) {
		std::cout << "dim " << n << ", blk " << b << ", nnz " << Gen.nnz(M) << std::endl;

		// det by b3if
		TW.clear(); 
		TW.start();
		IFD.largestInvariantFactors(result, MP, t, probability);
		TW.stop();
		double iftime = TW.usertime();

		typename PolynomialRing::Coeff d1;
		Polynomial pol;
		helper.computeDet(pol, result);
		d1 = pol[0];
		if (R.deg(pol) == n or R.getCoeffField().isZero(d1)) std::cout << "b3 det ";
		else std::cout << "bogus det ";
		if (pre == 2 and n & 1) R.getCoeffField().negin(d1);
		R.getCoeffField().write(std::cout, d1) << " in time " << iftime << std::endl;

		// det by other means
		Element d2; F.init(d2);
		TW.clear(); 
		TW.start();
		Method::Wiedemann Meth;
		LinBox::det(d2,MP, Meth);
		TW.stop();
		double dtime = TW.usertime();
		F.write(std::cout << "wied det ", d2)
			<< " in time " << dtime << std::endl;
		
	}

	
	std::cout << "nontriv: " << helper.countInvariantFactors(result)
		<< ", minpoldeg: " << helper.ring().deg(result.back()) << std::endl;

	if (outFile == "") {
		helper.writeInvariantFactors(std::cout, result, 0);
	} else {
		
	  //if (helper.degInvariantFactors(result) > M.rowdim()) {
		std::ofstream out(outFile);
		helper.writeInvariantFactors(out, result);
		out.close();
	  //}
	}
	
	return 0;
}


