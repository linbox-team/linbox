#include "linbox/linbox-config.h"

#include <algorithm>
#include <iostream>
#include <vector>

#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/frobenius-large.h"

#include "givaro/givtimer.h"

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

void readMatrix(SparseMat &M, const std::string matrixFile) {
	std::ifstream iF(matrixFile);
	M.read(iF);
	M.finalize();
	iF.close();
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
void test(const std::string &outFile, const Ring &R, const Blackbox &A, size_t k) {
	typedef typename Ring::Element Poly;
	
	FrobeniusLarge<Ring> FLD(R);
	
	std::vector<Poly> fs;
	std::vector<size_t> ms;
	time1([&](){FLD.solve(fs, ms, A, k);});
	std::cout << std::endl;
	
	if (outFile != "") {
		std::ofstream out(outFile);
		for (size_t i = 0; i < fs.size(); i++) {
			Poly f;
			R.write(out << ms[i] << " ", R.monic(f, fs[i])) << std::endl;
		}
		out.close();
	}
}

int main(int argc, char** argv) {
	size_t p = 3;
	size_t e = 1;
	size_t k = 0;
	std::string matrixFile, outFile;
	int seed = time(NULL);

	static Argument args[] = {
		{ 'k', "-k K", "Number of invariant factors to compute", TYPE_INT, &k},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'e', "-e E", "Size of extension field", TYPE_INT, &e},
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'o', "-o O", "Name of file for output", TYPE_STR, &outFile},
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
	readMatrix(M, matrixFile);
	assert(M.rowdim() == M.coldim());
	size_t n = M.rowdim();
	assert(e >= 1);
	
	if (e == 1) {
		std::cout << "Not implemented" << std::endl;
	} else {
		typedef NTL_zz_pE ExtField;
		typedef NTL_zz_pEX ExtRing;
		typedef SparseMatrix<ExtField, SparseMatrixFormat::CSR> ExtBlackbox;
		
		ExtField EF(p, e);
		ExtRing ER(EF);
		
		ExtBlackbox EM(EF, n, n);
		copy(EM, M);
		
		test(outFile, ER, EM, k);
	}
	
	return 0;
}
