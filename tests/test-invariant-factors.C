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
#include "linbox/blackbox/transpose.h"

#include "linbox/algorithms/poly-smith-form.h"
#include "linbox/algorithms/invariant-factors.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"

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

class TestInvariantFactorsHelper {
public:
	const size_t _p;
	const Field F;
	const Ring R;
	const MatrixDomain<Field> MD;
	
	TestInvariantFactorsHelper(size_t p) : _p(p), F(p), R(p), MD(F) {};
	
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
		
		// M = DA
		Field1 F(A.field());
		typename Field1::RandIter RI(F);
		SparseMatrix<Field1, Rep> M(F, A.rowdim(), A.coldim());
		std::cout << "scaling:" << std::flush;
		time1([&](){
			for (size_t i = 0; i < M.coldim(); i++) {
				typename Field1::Element factor;
				do {
					RI.random(factor);
				} while (F.isZero(factor));
				
				for (size_t j = 0; j < M.rowdim(); j++) {
					typename Field1::Element elm, scaledElm;
					A.getEntry(elm, i, j);
					if (!F.isZero(elm)) {
						F.mul(scaledElm, elm, factor);
						M.setEntry(i, j, scaledElm);
					}
				}
			}
			M.finalize();
		});
		
		if (precond == 1) { // D A -- Determinant preconditioner
			wiedemann(outFile, M);
		} else if (precond == 2) { // A^T D A -- Rank preconditioner
			Transpose<SparseMatrix<Field1, Rep>> T(A);
			Compose<Transpose<SparseMatrix<Field1, Rep>>, SparseMatrix<Field1, Rep>>  C(T, M);
			
			wiedemann(outFile, C);
		}
	}
	
	template<class Blackbox1, class Blackbox2>
	void copy(Blackbox1 &A, const Blackbox2 &B) const {
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

	template<class Blackbox>
	void computeGenerator(
		std::vector<BlasMatrix<typename Blackbox::Field>> &gen,
		const Blackbox &M,
		size_t b,
		int earlyTerm = 10) const
	{
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
	
	void extCoppersmith(std::string outFile, SparseMat &M, size_t b, size_t extend, Polynomial modulus1) const {
		typedef NTL_zz_pE ExtField;
		typedef NTL_zz_pEX ExtRing;
		typedef SparseMatrix<ExtField, SparseMatrixFormat::CSR> FBlackbox;
		
		ExtField EF((uint64_t) _p, extend);
		ExtRing ER(EF);
		FBlackbox EM(EF, M.rowdim(), M.coldim());
		copy(EM, M);
		
		typename ExtRing::Element modulus = NTL::conv<NTL::zz_pEX>(modulus1);
		
		std::vector<BlasMatrix<ExtField>> minpoly;
		time1([&](){computeGenerator(minpoly, EM, b);});
		
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
		
		SmithFormKannanBachemDomain<ExtRing> SFD(ER);
		
		std::vector<typename ExtRing::Element> result;
		time1([&](){SFD.solveIliopoulos(result, G, modulus);});
		
		size_t total = 0;
		for (auto it = result.begin(); it != result.end(); it++) {
			if (ER.deg(*it) > 1) {
				total++;
			}
		}
		std::cout << " \t" << total << std::endl;
		
		if (outFile != "") {
			std::ofstream out(outFile);
			std::for_each(result.begin(), result.end(), [&](const typename ExtRing::Element &v) {
				typename ExtRing::Element f;
				ER.write(out, ER.monic(f, v)) << std::endl;
			});
			out.close();
		}
		
		/*
		PolySmithFormDomain<ExtRing> PSFD(ER);
		std::cout << PSFD.detLimit(G) << "\t" << std::flush;
		typename ExtRing::Element det;
		time1([&](){PSFD.detLocalX(det, G);});
		
		std::vector<typename ExtRing::Element> result;
		time1([&](){PSFD.solve(result, G, det);});
		
		size_t total = 0;
		for (auto it = result.begin(); it != result.end(); it++) {
			if (ER.deg(*it) > 1) {
				total++;
			}
		}
		std::cout << " \t" << total << std::endl;
		
		if (outFile != "") {
			std::ofstream out(outFile);
			std::for_each(result.begin(), result.end(), [&](const typename ExtRing::Element &v) {
				typename ExtRing::Element f;
				ER.write(out, ER.monic(f, v)) << std::endl;
			});
			out.close();
		}
		*/
	}
	//*/
}; // End of TestInvariantFactorsHelper

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

void randomPerm(SparseMat & P)
{
	size_t n = P.rowdim();
	std::vector<size_t> Perm;
	randomVec(Perm, n);
	for (size_t i = 0; i < n; ++i) P.setEntry(i, Perm[i], P.field().one);
	P.finalize();
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

template<class Sp>
void permuteCols(Sp & MP, const std::vector<size_t> & P, const Sp& M) {
	Element x;
	M.field().init(x);
	for (size_t i = 0; i < M.rowdim(); ++i) {
		for (size_t j = 0; j < M.coldim(); ++j) {
			if (not M.field().isZero(M.getEntry(x,i,j))) {
				MP.setEntry(i,P[j],x);
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

void randomSparse(SparseMat &M, size_t s) {
	Field F(M.field());
	
	for (size_t i = 0; i < M.rowdim(); i++) {
		for (size_t k = 0; k < s; k++) {
			size_t j = rand() % M.coldim();
			
			Element elm;
			randomNonzero(F, elm);
			
			M.setEntry(i, j, elm);
		}
	}
	M.finalize();
}

void randomToeplitz(SparseMat &T, bool upper, bool scaled = false) {
	Field F(T.field());
	typename Field::RandIter RI(F);
	
	SparseMatrix<Field, SparseMatrixFormat::SMM> S(F, T.rowdim(), T.rowdim());
	
	Element elm;
	do {
		RI.random(elm);
	} while (F.isZero(elm));
		
	for (size_t i = 0; i < T.rowdim(); i++) {
		S.setEntry(i, i, elm);
	}
	
	
	for (size_t i = 1; i < T.rowdim(); i++) {
		std::cout << "i: " << i << std::endl;
		
		RI.random(elm);
		if (F.isZero(elm)) {
			continue;
		}
		
		for (size_t j = 0; i+j < T.rowdim(); j++) {
			if (upper) {
				S.setEntry(j, i+j, elm);
			} else {
				S.setEntry(i+j, j, elm);
			}
		}
	}
	S.finalize();
	
	for (size_t i = 0; i < T.rowdim(); i++) {
		do {
			RI.random(elm);
		} while (F.isZero(elm));
		
		for (size_t j = 0; j < T.rowdim(); j++) {
			Element tmp;
			S.getEntry(tmp, i, j);
			F.mulin(tmp, elm);
			S.setEntry(i, j, tmp);
		}
	}
	S.finalize();
	
	for (size_t i = 0; i < T.rowdim(); i++) {
		for (size_t j = 0; j < T.rowdim(); j++) {
			Element tmp;
			S.getEntry(tmp, i, j);
			T.setEntry(i, j, tmp);
		}
	}
	T.finalize();
}

void randomFatDiag(SparseMat &T, size_t s, bool upper, bool randomDiag = false) {
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
	
	for (size_t i = 1; i <= s; i++) {
		for (size_t k = 0; k < T.rowdim(); k++) {
			if (i + k >= T.rowdim()) {
				continue;
			}
			
			Element elm;
			randomNonzero(F, elm);
			if (upper) {
				T.setEntry(i + k, k, elm);
			} else {
				T.setEntry(k, i + k, elm);
			}
		}
	}
	
	T.finalize();
}

double randf() {
	return (double)(rand() * 1.0 / RAND_MAX);
}

void randomWiedemann(SparseMat &T, double k, bool right) {
	Field F(T.field());
	typename Field::RandIter RI(F);
	size_t p = F.cardinality();
	size_t n = T.rowdim();
		double pr1 = (p-1.0)/p;
	
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			double pr2 = k * log(1.0 * n) / ((right ? j : i) + 1.0);
			
			Element elm;
			F.assign(elm, F.zero);
			if (pr1 < pr2) {
				RI.random(elm);
			} else if (randf() <= pr2) {
				randomNonzero(F, elm);
			}
			
			if (F.isZero(elm)) {
				continue;
			}
			
			T.setEntry(i, j, elm);
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
	
	size_t n = M.rowdim();
	if (n <= 20) M.write(std::cout) << std::endl;
	
	SparseMat PreR(F, n, n);
	SparseMat PreL(F, n, n);
	
	bool precondL = false;
	bool precondR = false;
	
	if (abs(precond) == 1) {
		randomFatDiag(PreR, s, true);
		precondR = true;
	} else if (abs(precond) == 2) {
		randomFatDiag(PreR, s, true, precond < 0);
		randomFatDiag(PreL, s, false, precond < 0);
		precondL = precondR = true;
	} else if (abs(precond) == 3) {
		randomFatDiag(PreR, s, false, precond < 0);
		randomFatDiag(PreL, s, true, precond < 0);
		precondL = precondR = true;
	} else if (abs(precond) == 4) {
		randomTriangular(PreR, s, true, precond < 0);
		randomTriangular(PreL, s, false, precond < 0);
		precondL = precondR = true;
	} else if (abs(precond) == 5) {
		randomTriangular(PreR, s, false, precond < 0);
		randomTriangular(PreL, s, true, precond < 0);
		precondL = precondR = true;
	} else if (abs(precond) == 6) {
		randomToeplitz(PreR, false);
		randomToeplitz(PreL, true, precond < 0);
		precondL = precondR = true;
	} else if (precond == 7) {
		randomSparse(PreR, s);
		precondR = true;
	} else if (precond == 8) {
		double k = 0.5;
		randomWiedemann(PreR, k, true);
		randomWiedemann(PreL, k, false);
		
		for (size_t i = 0; i < 20; i++) {
			for (size_t j = 0; j < 20; j++) {
				F.write(std::cout, PreL.getEntry(i,j)) << " ";
			}
			std::cout << std::endl;
		}
		precondL = precondR = true;
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
	
	if (extend > 1) {
		helper.extCoppersmith(outFile, M, b, extend, result[result.size() - modIndex - 1]);
		return 0;
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
