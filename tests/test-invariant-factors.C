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

#include "linbox/blackbox/blockbb.h"
#include "linbox/blackbox/block-compose.h"
#include "linbox/blackbox/fflas-csr.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sum.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/toeplitz.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/blackbox/polynomial.h"

#include "linbox/algorithms/poly-smith-form.h"
#include "linbox/algorithms/invariant-factors.h"
#include "linbox/algorithms/smith-form-kannan-bachem.h"

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

template<class Ring1>
void writeLifs(Ring1 R, std::string filename, std::vector<typename Ring1::Element> &result) {
	std::ofstream out(filename);
	std::for_each(result.begin(), result.end(), [&](const typename Ring1::Element &v) {
		typename Ring1::Element f;
		R.write(out, R.monic(f, v)) << std::endl;
	});
	out.close();
}

void readMatrix(SparseMat &M, const std::string matrixFile) {
	std::ifstream iF(matrixFile);
	M.read(iF);
	M.finalize();
	iF.close();
}

int main(int argc, char** argv) {
	size_t p = 3;
	size_t extend = 1;
	size_t b = 4, m = 0, l = 1;
	size_t modIndex = 0;
	size_t exponent = 0;
	
	int precond = 0;
	size_t s = 0;
	size_t perm = 0;
	size_t k = 0;
	size_t spmv = 0;
	
	std::string matrixFile;
	std::string outFile;
	
	int seed = time(NULL);
	
	size_t alg = 0;

	static Argument args[] = {
		{ 'a', "-a A", "Algs to run after BM", TYPE_INT, &alg},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'e', "-e E", "Extension field exponent (p^e)", TYPE_INT, &extend},
		{ 'v', "-v V", "Enable spmv instead of spmm", TYPE_INT, &spmv},
		
		{ 'b', "-b B", "Block size", TYPE_INT, &b},
		{ 'm', "-m M", "Step size for iterative lifs", TYPE_INT, &m},
		{ 'l', "-l L", "Limit for iterative lifs", TYPE_INT, &l},
		{ 't', "-t T", "Use t-th LIF as modulus", TYPE_INT, &modIndex},
		{ 'd', "-d D", "Compute rank of ((x-1)^d)(A)", TYPE_INT, &exponent},
		
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		{ 'k', "-k K", "Compute the (k+1)-th invariant factor", TYPE_INT, &k},
		{ 's', "-s S", "Number of nonzeros in random triangular preconditioner", TYPE_INT, &s},
		{ 'c', "-c C", "Choose what preconditioner to apply", TYPE_INT, &precond},
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
	
	typename Field::RandIter RI(F);
	RandomDenseMatrix<typename Field::RandIter, Field> RDM(F, RI);
	
	SparseMat M(F);                                           
	readMatrix(M, matrixFile);
	
	assert(M.rowdim() == M.coldim());
	
	FflasCsr<Field> FM(&M);
	std::vector<Polynomial> result;
	if (precond == 1) { // determinant
		Element det;
		if (IFD.det(det, FM, b)) {
			F.write(std::cout << "det: ", det) << std::endl;
		} else {
			std::cout << "det failed" << std::endl;
		}
	} else if (precond == 2) { // rank
		size_t rank;
		if (IFD.rank(rank, FM, b)) {
			std::cout << "rank: " << rank << std::endl;
		} else {
			std::cout << "rank failed" << std::endl;
		}
	} else {
		time1([&](){IFD.largestInvariantFactors(result, FM, b);});
	}
	
	if (outFile != "") {
		writeLifs(R, outFile, result);
	}
	
	return 0;
}
