#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "givaro/givtimer.h"

#include "givaro/extension.h"
#include "givaro/gfqext.h"
#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"
#include "linbox/ring/polynomial-ring.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/matrix/random-matrix.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/wiedemann.h"

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
		
		Field1 F = result.field();
		std::ofstream out(outFile);
		
		for (size_t i = 0; i < result.size(); i++) {
			if (F.isZero(result[i])) {
				continue;
			}
			
			if (i == 0 || !F.isOne(result[i])) {
				F.write(out, result[i]);
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
	
	template<class Field1, class Rep>
	void wiedemann(std::string outFile, SparseMatrix<Field1, Rep> &M) const {
		typedef typename Field1::RandIter RandIter;
		typedef BlackboxContainer<Field1, SparseMatrix<Field1, Rep>> Sequence;
		
		BlasVector<Field1> phi(M.field());
		unsigned long deg;
		
		time1([&](){
			RandIter RI(M.field());
			Sequence seq(&M, M.field(), RI);
			MasseyDomain<Field1, Sequence> WD(&seq, 10);
			
			WD.minpoly(phi, deg);
		});
		writeInvariantFactor(outFile, phi);
		
		std::cout << (phi.size() - 1) << std::endl;
	}
	
	void extWiedemann(std::string outFile, SparseMat &M, size_t extend) const {
		if (extend == 1) {
			wiedemann(outFile, M);
			return;
		}
		
		typedef Givaro::GFqDom<int64_t> ExtField;
		typedef typename SparseMat::template rebind<ExtField>::other FBlackbox;
		
		// uint64_t extend1 = (uint64_t)Givaro::FF_EXPONENT_MAX((uint64_t)p, (uint64_t)LINBOX_EXTENSION_DEGREE_MAX);
		uint64_t extend1 = extend;
		//std::cout << extend1 << "\t" << std::flush;
		
		ExtField EF((uint64_t) _p, extend1);
		FBlackbox EM(M, EF);
		
		wiedemann(outFile, EM);
	}
	
	//*
	template<class Field1, class Rep>
	void coppersmith(SparseMatrix<Field1, Rep> &M, size_t b) const {
		InvariantFactors<Field1, Ring> IFD(M.field(), R);
			
		// Generate random left and right projectors
		std::vector<BlasMatrix<Field1>> minpoly;
		time1([&](){IFD.computeGenerator(minpoly, M, b);});
		std::cout << std::endl;
		
		/*
		typedef PolynomialRing<Field1> Ring1;
		typedef typename Ring1::Element Poly;
		
		Ring1 ER(M.field());
		MatrixDomain<Ring1> ERMD(ER);
		BlasMatrix<Ring1> G(ER, b, b);
		for (size_t i = 0; i < b; i++) {
			for (size_t j = 0; j < b; j++) {
				std::vector<typename Field1::Element> lst;				
				for (size_t k = 0; k < minpoly.size(); k++) {
					lst.push_back(minpoly[k].getEntry(i, j));
				}
				
				Poly tmp;
				ER.init(tmp, lst);
				
				G.setEntry(i, j, tmp);
				ER.write(std::cout, G.getEntry(i, j)) << std::endl;
			}
			std::cout << std::endl;
		}
		
		SmithFormKannanBachemDomain<Ring1> SFD(ER);
		std::vector<typename Ring1::Element> result;
		SFD.solveTextBook(result, G);
		
		for (size_t i = 0; i < result.size(); i++) {
			ER.write(std::cout, result[i]) << std::endl;
		}
		//*/
	}
	
	void extCoppersmith(SparseMat &M, size_t b, size_t extend) const {
		//*
		typedef Givaro::GFqDom<int64_t> ExtField;
		typedef typename SparseMat::template rebind<ExtField>::other FBlackbox;
		
		ExtField EF((uint64_t) _p, extend);
		FBlackbox EM(M, EF);
		
		coppersmith(EM, b);
		//*/
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
void permuteRows(Sp & MP, std::vector<size_t> & P, Sp& M) {
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
void permuteCols(Sp & MP, std::vector<size_t> & P, Sp& M) {
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

void randomNonzero(const Field &F, Element &elm) {
	typename Field::RandIter RI(F);
	do {
		RI.random(elm);
	} while (F.isZero(elm));
}

void randomTriangular(SparseMat &T, size_t s) {
	typename Field::RandIter RI(T.field());
	
	for (size_t i = 0; i < T.rowdim(); i++) {
		T.setEntry(i, i, T.field().one);
	}
	
	for (size_t r = 0; r < T.rowdim() - 1; r++) {
		for (size_t k = 0; k < s; k++) {
			size_t c = (rand() % (T.coldim() - r - 1)) + r + 1;
			
			Element elm;
			randomNonzero(T.field(), elm);
			
			T.setEntry(r, c, elm);
		}
	}
	
	T.finalize();
}

void randomLowerTriangular(SparseMat &T, size_t s) {
	typename Field::RandIter RI(T.field());
	
	for (size_t i = 0; i < T.rowdim(); i++) {
		T.setEntry(i, i, T.field().one);
	}
	
	for (size_t r = 0; r < T.rowdim() - 1; r++) {
		for (size_t k = 0; k < s; k++) {
			size_t c = (rand() % (T.coldim() - r - 1)) + r + 1;
			
			Element elm;
			randomNonzero(T.field(), elm);
			
			T.setEntry(c, r, elm);
		}
	}
	
	T.finalize();
}

int main(int argc, char** argv) {
	size_t p = 3;
	size_t extend = 1;
	size_t b = 4;
	
	size_t n = 1000;
	size_t s = 10;
	double fillIn = 0;
	
	std::string matrixFile;
	std::string outFile;
	
	int seed = time(NULL);
	
	size_t alg = 0;

	static Argument args[] = {
		{ 'a', "-a A", "Algs to run after BM", TYPE_INT, &alg},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'e', "-e E", "Extension field exponent (p^e)", TYPE_INT, &extend},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		
		{ 'n', "-n N", "Dimension of input", TYPE_INT, &n},
		{ 's', "-s S", "Number of nonzeros in random triangular preconditioner", TYPE_INT, &s},
		{ 'i', "-i I", "Sparsity of input matrix", TYPE_DOUBLE, &fillIn},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	
	//*
	if (matrixFile == "") {
		std::cout << "an input matrix must be provided" << std::endl;
		return -1;
	}
	//*/
	
	srand(seed);

	Field F(p);
	Ring R(p);
	
	//*
	SparseMat M(F);                                           
	{
		std::ifstream iF(matrixFile);
		M.read(iF);
		M.finalize();
		iF.close();
	}
	//*/
	
	/*
	SparseMat OldM(F, n, n);
	{
		SparseMatrixGenerator<Field, Ring> Gen(F, R);
		
		std::vector<integer> v;
		Polynomial x, xm1;
		R.init(x, v = {0, 1});
		R.assign(xm1, x);
		R.subin(xm1, R.one);
		
		Polynomial det;
		std::vector<Polynomial> fs;
		Gen.addTriangle(fs, n/2, 6.0 / n, x);
		Gen.addTriangle(fs, n, 1, xm1);
		Gen.generate(OldM, det, fs, fillIn);
	}
	SparseMat M(OldM);
	//*/
	
	assert(M.rowdim() == M.coldim());
	n = M.rowdim();
	if (n <= 20) M.write(std::cout) << std::endl;
	
	SparseMat PreR(F, n, n);
	randomTriangular(PreR, s);
	
	SparseMat PreL(F, n, n);
	randomLowerTriangular(PreL, s);
	
	TestInvariantFactorsHelper helper(p);
	std::cout << n << "\t" << M.size() << "\t" << b << "\t" << std::flush;
	
	if (b == 1) {
		helper.extWiedemann(outFile, M, extend);
		return 0;
	}
	
	if (extend > 1) {
		helper.extCoppersmith(M, b, extend);
		return 0;
	}
	
	InvariantFactors<Field, Ring> IFD(F, R);
	PolySmithFormDomain<Ring> PSFD(R);
		
	// Generate random left and right projectors
	std::vector<BlasMatrix<Field>> minpoly;
	
	if (s == 0) {
		time1([&](){IFD.computeGenerator(minpoly, M, b);});
	} else {
		time1([&](){IFD.computeGenerator(minpoly, PreL, M, PreR, b);});
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
