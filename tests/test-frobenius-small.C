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

void readMatrix(SparseMat &M, const std::string matrixFile) {
	std::ifstream iF(matrixFile);
	M.read(iF);
	M.finalize();
	iF.close();
}

int main(int argc, char** argv) {
	size_t p = 3;
	size_t b = 8;
	size_t k = 4;
	std::string matrixFile;
	int seed = time(NULL);

	static Argument args[] = {
		{ 'b', "-b B", "Block size for minpolyvec", TYPE_INT, &b},
		{ 'k', "-k K", "Number of invariant factors to compute", TYPE_INT, &k},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
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
	
	SparseMat MT(F, n, n);
	for (size_t i = 0; i < n; i++) {
		for (size_t j = 0; j < n; j++) {
			MT.setEntry(j, i, M.getEntry(i, j));
		}
	}
	MT.finalize();
	
	FrobeniusSmall<Field, Ring> FSD(F, R);
	
	std::vector<Polynomial> fs;
	FSD.solve(fs, M, k);
	
	for (size_t i = 0; i < fs.size(); i++) {
		R.write(std::cout, fs[i]) << std::endl;
	}
	
	return 0;
}
