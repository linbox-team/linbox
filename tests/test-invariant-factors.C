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


using namespace LinBox;

typedef Givaro::Modular<double> Field;
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
		*/
	}
	
	void extCoppersmith(SparseMat &M, size_t b, size_t extend) const {
		/*
		typedef Givaro::GFqDom<int64_t> ExtField;
		typedef typename SparseMat::template rebind<ExtField>::other FBlackbox;
		
		ExtField EF((uint64_t) _p, extend);
		FBlackbox EM(M, EF);
		
		coppersmith(EM, b);
		//*/
	}
	//*/
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
	
	if (matrixFile == "") {
		std::cout << "an input matrix must be provided" << std::endl;
		return -1;
	}
	
	srand(seed);

	Field F(p);
	Ring R(p);
	
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
	std::cout << n << "\t" << M.size() << "\t" << b << "\t" << std::flush;
	
	if (b == 1) {
		if (extend > 1) {
			helper.extWiedemann(outFile, M, extend);
		} else {
			helper.wiedemann(outFile, M);
		}
		
		return 0;
	}
	
	PolySmithFormDomain<Ring> PSFD(R);
	
	if (extend > 1) {
		helper.extCoppersmith(M, b, extend);
		
		return 0;
	}
	
	InvariantFactors<Field, Ring> IFD(F, R);
		
	// Generate random left and right projectors
	std::vector<BlasMatrix<Field>> minpoly;
	time1([&](){IFD.computeGenerator(minpoly, M, b);});
	
	// Convert to matrix with polynomial entries
	PolyMatrix G(R, b, b);
	IFD.convert(G, minpoly);
	
	// Compute smith form of generator
	std::cout << PSFD.detLimit(G) << "\t" << std::flush;
	
	Polynomial det, mp;
	std::vector<Polynomial> result;
	time2([&](){return PSFD.dixon(mp, G);});
	time1([&](){PSFD.detPopov(det, G);});
	time1([&](){PSFD.detLocalX(det, G);});
	time1([&](){PSFD.solve(result, G, det);});
	std::cout << std::endl;
	
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
