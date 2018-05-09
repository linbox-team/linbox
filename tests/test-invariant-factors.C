#include "linbox/linbox-config.h"

#include <algorithm>
#include <iostream>
#include <vector>

#include "givaro/givtimer.h"

#include "givaro/gfqext.h"
#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"
#include "linbox/ring/polynomial-ring.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/wiedemann.h"

#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sum.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/toeplitz.h"
#include "linbox/blackbox/diagonal.h"

#include "linbox/algorithms/poly-smith-form.h"
#include "linbox/algorithms/invariant-factors.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"

#include "linbox/algorithms/frobenius-small.h"

#include "sparse-matrix-generator.h"

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;

typedef NTL_zz_pX Ring;
typedef typename Ring::Element Polynomial;
typedef BlasMatrix<Ring> PolyMatrix;

Givaro::Timer TW;

void time2(std::function<bool()> fun) {
	TW.clear();
	TW.start();
	
	bool s = fun();
	
	TW.stop();
	std::cout << TW.usertime() << (s ? "" : "(failed)") << "\t" << std::flush;
}

void time1(std::function<void()> fun) {
	time2([fun](){fun(); return true;});
}

template<class Field1>
void writeInvariantFactor(
	const std::string outFile,
	const BlasVector<Field1> &result) {

	if (outFile == "") {
		return;
	}
	
	Field1 F(result.field());
	std::ofstream out(outFile);
	
	for (size_t i = 0; i < result.size(); i++) {
		F.write(out, result[i]) << std::endl;
	}
	
	out.close();
}

template<class Blackbox1, class Blackbox2>
void copy(Blackbox1 &A, const Blackbox2 &B) {
	typedef typename Blackbox1::Field Field1;
	typedef typename Blackbox2::Field Field2;
	
	for (size_t i = 0; i < B.rowdim(); i++) {
		for (size_t j = 0; j < B.coldim(); j++) {
			typename Field2::Element elm;
			integer tmp;
			typename Field1::Element oelm;
			
			B.getEntry(elm, i, j);
			B.field().convert(tmp, elm);
			A.field().init(oelm, tmp);
			A.setEntry(i, j, oelm);
		}
	}
	A.finalize();
}

template<class Ring, class Blackbox>
void wiedemann1(BlasVector<Ring> &result, const Blackbox &M) {
	typedef typename Blackbox::Field Field1;
	typedef typename Field1::RandIter RandIter;
	typedef BlackboxContainer<Field1, Blackbox> Sequence;
	
	BlasVector<Field1> phi(M.field());
	unsigned long deg;
	
	time1([&](){
		RandIter RI(M.field());
		Sequence seq(&M, M.field(), RI);
		MasseyDomain<Field1, Sequence> WD(&seq, 20);
		
		WD.minpoly(phi, deg);
	});
	std::cout << (phi.size() - 1) << std::endl;
	result.field().init(result.refEntry(0), phi);
}

template<class Blackbox>
void computeGenerator(
	std::vector<BlasMatrix<typename Blackbox::Field>> &gen,
	const Blackbox &M,
	size_t b,
	int earlyTerm = 10) {
	typedef typename Blackbox::Field Field1;
	typedef typename Field1::RandIter RandIter1;
	typedef MatrixDomain<Field1> MatrixDom1;
	typedef typename MatrixDom1::OwnMatrix Matrix1;
	
	Field1 EF(M.field());
	RandIter1 RI(EF);
	RandomDenseMatrix<RandIter1, Field1> RDM(EF, RI);
	MatrixDom1 MD(EF);
	
	size_t n = M.rowdim();
	Matrix1 U(EF, b, n);
	Matrix1 V(EF, n, b);
	
	RDM.random(U);
	RDM.random(V);
	
	typedef BlackboxBlockContainer<Field1, Blackbox> Sequence;
	Sequence blockSeq(&M, EF, U, V);
	BlockCoppersmithDomain<MatrixDom1, Sequence> coppersmith(MD, &blockSeq, earlyTerm);
	
	coppersmith.right_minpoly(gen);
}

template<class Blackbox>
void wiedemannb(
	BlasVector<NTL_zz_pEX> &result,
	const Blackbox &M,
	const typename NTL_zz_pEX::Element &modulus) {

	typedef NTL_zz_pE ExtField;
	typedef NTL_zz_pEX ExtRing;
	
	ExtField EF(result.field().getCoeffField());
	ExtRing ER(result.field());
	
	size_t b = result.size();
		
	std::vector<BlasMatrix<ExtField>> minpoly;
	time1([&](){computeGenerator(minpoly, M, b);});
	
	BlasMatrix<ExtRing> G(ER, b, b);
	for (size_t i = 0; i < b; i++) {
		for (size_t j = 0; j < b; j++) {
			typename ExtRing::Element g;
			ER.init(g);
			
			for (size_t k = 0; k < minpoly.size(); k++) {
				ER.setCoeff(g, k, minpoly[k].getEntry(i, j));
			}
			
			G.setEntry(i, j, g);
		}
	}
	
	std::vector<typename ExtRing::Element> tmpResult;
	SmithFormKannanBachemDomain<ExtRing> SFD(ER);
	time1([&](){SFD.solveIliopoulos(tmpResult, G, modulus);});	
	
	for (size_t i = 0; i < b; i++) {
		result.setEntry(i, tmpResult[i]);
	}
}

class TestInvariantFactorsHelper {
public:
	const size_t _p;
	
	TestInvariantFactorsHelper(size_t p) : _p(p) {};
	
	template<class Field1>
	void writeInvariantFactor(std::string outFile, BlasVector<Field1> result) const {
		if (outFile == "") {
			return;
		}
		
		Field1 F(result.field());
		std::ofstream out(outFile);
		
		for (size_t i = 0; i < result.size(); i++) {
			if (F.isZero(result[i])) {
				continue;
			}
			
			if (i == 0 || !F.isOne(result[i])) {
				F.write(out << "(", result[i]) << ")";
			}
			if (i > 0) {
				out << (F.isOne(result[i]) ? "" : "*") << "x^" << i;
			}
			if (i < result.size() - 1) {
				out << "+";
			}
		}
		out << std::endl;
		
		out.close();
	}
	
	template<class Blackbox>
	void wiedemann(std::string outFile, Blackbox &M) const {
		typedef typename Blackbox::Field Field1;
		typedef typename Field1::RandIter RandIter;
		typedef BlackboxContainer<Field1, Blackbox> Sequence;
		
		BlasVector<Field1> phi(M.field());
		unsigned long deg;
		
		time1([&](){
			RandIter RI(M.field());
			Sequence seq(&M, M.field(), RI);
			MasseyDomain<Field1, Sequence> WD(&seq, 20);
			
			WD.minpoly(phi, deg);
		});
		std::cout << (phi.size() - 1) << std::endl;
		writeInvariantFactor(outFile, phi);
	}
	
	template<class Field1, class Rep>
	void wiedemann(std::string outFile, SparseMatrix<Field1, Rep> &A, size_t precond) const {
		if (precond == 0) {
			wiedemann(outFile, A);
			return;
		}
		
		typedef SparseMatrix<Field1, Rep> Sp;
		typedef Diagonal<Field1> Diag;
		typedef Compose<Diagonal<Field1>, Sp> Scaled;
		Diag D(A.field(), A.rowdim());
		Scaled M(D, A);
		
		if (precond == 1) { // D A -- Determinant preconditioner
			wiedemann(outFile, M);
		} else if (precond == 2) { // D A^T D A -- Rank preconditioner
			Transpose<Sp> T(A);
			Diag E(A.field(), A.rowdim());
			
			typedef Compose<Diagonal<Field1>, Transpose<Sp>> ScaledTranspose;
			ScaledTranspose W(E, T);
			
			Compose<ScaledTranspose, Scaled>  C(W, M);
			
			wiedemann(outFile, C);
		}
	}
	
	void extWiedemann(std::string outFile, SparseMat &M, size_t extend, size_t precond) const {
		if (extend <= 1) {
			wiedemann(outFile, M);
			return;
		}
		
		uint64_t maxZechExt = Givaro::FF_EXPONENT_MAX((uint64_t)_p, (uint64_t)LINBOX_EXTENSION_DEGREE_MAX);
		if (extend <= maxZechExt) {
			typedef Givaro::GFqDom<int64_t> ExtField;
			typedef typename SparseMat::template rebind<ExtField>::other FBlackbox;
			
			ExtField EF((uint64_t) _p, extend);
			FBlackbox EM(M, EF);
			EM.finalize();
			
			wiedemann(outFile, EM, precond);
		} else {
			typedef NTL_zz_pE ExtField;
			typedef SparseMatrix<ExtField, SparseMatrixFormat::CSR> FBlackbox;
			
			ExtField EF(_p, extend);
			FBlackbox EM(EF, M.rowdim(), M.coldim());
			copy(EM, M);
			
			wiedemann(outFile, EM, precond);
		}
	}
}; // End of TestInvariantFactorsHelper

void randomVec(std::vector<size_t> & V, size_t n)  {
	V.resize(n);
	for (size_t i = 0; i < n; ++i) V[i] = i;
	// Knuth construction
	for (size_t i = 0; i < n-1; ++i) {
		size_t j = i + rand()%(n-i);
		std::swap(V[i], V[j]);
	}
}

template<class Sp>
void permuteRows(Sp & MP, const std::vector<size_t> & P, const Sp& M) {
	Element x;
	M.field().init(x);
	for (size_t i = 0; i < M.rowdim(); ++i) {
		for (size_t j = 0; j < M.coldim(); ++j) {
			if (not M.field().isZero(M.getEntry(x,i,j))) {
				MP.setEntry(P[i],j,x);
			}
		}
	}
	
	MP.finalize();
}

void readMatrix(SparseMat &M, const std::string matrixFile, bool perm) {
	if (!perm) {
		std::ifstream iF(matrixFile);
		M.read(iF);
		M.finalize();
		iF.close();
	} else {
		SparseMat OldM(M.field());
		std::ifstream iF(matrixFile);
		OldM.read(iF);
		OldM.finalize();
		iF.close();
		
		assert(M.rowdim() == M.coldim());
		size_t n = OldM.rowdim();
		M.resize(n, n);
		
		std::vector<size_t> v;
		randomVec(v, n);
		permuteRows(M, v, OldM);
	}
}

void randomNonzero(const Field &F, Element &elm) {
	typename Field::RandIter RI(F);
	do {
		RI.random(elm);
	} while (F.isZero(elm));
}

void randomTriangular(SparseMat &T, size_t s, bool upper, bool randomDiag = false) {
	Field F(T.field());
	typename Field::RandIter RI(F);
	
	for (size_t i = 0; i < T.rowdim(); i++) {
		if (randomDiag) {
			Element elm;
			randomNonzero(F, elm);
			T.setEntry(i, i, elm);
		} else {
			T.setEntry(i, i, F.one);
		}
	}
	
	for (size_t r = 0; r < T.rowdim() - 1; r++) {
		for (size_t k = 0; k < s; k++) {
			size_t c = (rand() % (T.coldim() - r - 1)) + r + 1;
			
			Element elm;
			randomNonzero(F, elm);
			
			if (upper) {
				T.setEntry(r, c, elm);
			} else {
				T.setEntry(c, r, elm);
			}
		}
	}
	
	T.finalize();
}

int main(int argc, char** argv) {
	size_t p = 3;
	size_t extend = 1;
	size_t b = 4;
	size_t modIndex = 0;
	
	int precond = 0;
	size_t s = 0;
	size_t perm = 0;
	size_t k = 0;
	
	std::string matrixFile;
	std::string outFile;
	
	int seed = time(NULL);
	
	size_t alg = 0;

	static Argument args[] = {
		{ 'a', "-a A", "Algs to run after BM", TYPE_INT, &alg},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'e', "-e E", "Extension field exponent (p^e)", TYPE_INT, &extend},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
		{ 't', "-t T", "Use t-th LIF as modulus", TYPE_INT, &modIndex},
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		{ 'k', "-k K", "Compute the (k+1)-th invariant factor", TYPE_INT, &k},
		
		{ 's', "-s S", "Number of nonzeros in random triangular preconditioner", TYPE_INT, &s},
		{ 'c', "-c C", "Choose what preconditioner to apply", TYPE_INT, &precond},
		{ 'z', "-z Z", "Permute rows of input", TYPE_INT, &perm},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	if (matrixFile == "") {
		std::cout << "an input matrix must be provided" << std::endl;
		return -1;
	}
	
	srand(seed);

	Field F(p);
	Ring R(p);
	MatrixDomain<Field> MD(F);
	InvariantFactors<Field, Ring> IFD(F, R);
	PolySmithFormDomain<Ring> PSFD(R);
	
	SparseMat M(F);                                           
	readMatrix(M, matrixFile, perm != 0);
	
	assert(M.rowdim() == M.coldim());
	size_t n = M.rowdim();
	
	SparseMat PreR(F, n, n);
	SparseMat PreL(F, n, n);
	
	bool precondL = false;
	bool precondR = false;
	
	if (precond == 1) { // determinant
		randomTriangular(PreR, s, false, false);
		randomTriangular(PreL, s, true, false);
		precondL = precondR = true;
	} else if (precond == 2) { // rank
		randomTriangular(PreR, s, false, true);
		randomTriangular(PreL, s, true, true);
		precondL = precondR = true;
	} else if (precond == 3) { // k-th invariant
		typedef NTL_zz_pE ExtField;
		typedef NTL_zz_pEX ExtRing;
		typedef SparseMatrix<ExtField, SparseMatrixFormat::CSR> ExtBlackbox;
		
		ExtField EF(p, extend);
		ExtRing ER(EF);
		
		ExtBlackbox EM(EF, n, n);
		copy(EM, M);
		
		typedef typename ExtField::Element Coeff;
		typedef typename ExtField::RandIter RandIter;
		typedef typename ExtRing::Element Poly;
		
		RandIter RI(EF);
		
		Poly u, v;
		ER.init(u);
		ER.init(v);
		for (size_t i = 0; i < n+k-1; i++) {
			Coeff c1, c2;
			do {
				RI.random(c1);
			} while (EF.isZero(c1));
			ER.setCoeff(u, i, c1);
			
			do {
				RI.random(c2);
			} while (EF.isZero(c2));
			ER.setCoeff(v, i, c2);
		}
		
		typedef Toeplitz<ExtField, ExtRing> Toep;
		Toep U(ER, u, n, k);
		Toep V(ER, v, k, n);
		Compose<Toep, Toep> B(U, V);
		Sum<ExtBlackbox, Compose<Toep, Toep>> Mk(EM, B);
		
		BlasVector<ExtRing> result(ER, b);
		wiedemann1(result, EM);
		
		Poly minpoly;
		result.getEntry(minpoly, 0);
		ER.write(std::cout << "minpoly: ", minpoly) << std::endl;
		
		//*
		if (b == 1) {
			wiedemann1(result, Mk);
		} else {
			wiedemannb(result, Mk, minpoly);
			std::cout << std::endl;
		}
		//*/
		
		//*
		for (size_t i = 0; i < b; i++) {
			Poly fi;
			ER.gcd(fi, result[i], minpoly);
			result.setEntry(i, fi);
		}
		//*/
		
		ER.write(std::cout << "(k+1)-th LIF: ", result[result.size() - 1]) << std::endl;
		writeInvariantFactor(outFile, result);
		
		return 0;
	}
	
	TestInvariantFactorsHelper helper(p);
	std::cout << n << "\t" << M.size() << "\t" << b << "\t" << std::flush;
	
	if (b == 1) {
		helper.extWiedemann(outFile, M, extend, precond);
		return 0;
	}
		
	// Generate random left and right projectors
	std::vector<BlasMatrix<Field>> minpoly;
	
	if (precondL && precondR) {
		time1([&](){IFD.computeGenerator(minpoly, PreL, M, PreR, b);});
	} else if (precondR) {
		time1([&](){IFD.computeGenerator(minpoly, M, PreR, b);});
	} else if (precondL) {
		time1([&](){IFD.computeGenerator(minpoly, PreL, M, b);});
	} else {
		time1([&](){IFD.computeGenerator(minpoly, M, b);});
	}
	
	// Convert to matrix with polynomial entries
	PolyMatrix G(R, b, b);
	IFD.convert(G, minpoly);
	
	// Compute smith form of generator
	std::cout << PSFD.detLimit(G) << "\t" << std::flush;
	
	Polynomial det, mp;
	std::vector<Polynomial> result;
		
	if (alg & 8) {
		time2([&](){return PSFD.dixon(mp, G);});
	}
	if (alg & 4) {
		time1([&](){PSFD.detPopov(det, G);});
	}
	if (alg & 2) {
		time1([&](){PSFD.detLocalX(det, G);});
	}
	if ((alg & 1) && !R.isZero(det)) {
		time1([&](){PSFD.solve(result, G, det);});
	}
	
	size_t total = 0;
	for (auto it = result.begin(); it != result.end(); it++) {
		if (R.deg(*it) > 1) {
			total++;
		}
	}
	
	std::cout << " \t" << total << std::endl;
	
	if (outFile != "") {
		std::ofstream out(outFile);
		std::for_each(result.begin(), result.end(), [&](const Polynomial &v) {
			Polynomial f;
			R.write(out, R.monic(f, v)) << std::endl;
		});
		out.close();
	}
	
	return 0;
}
