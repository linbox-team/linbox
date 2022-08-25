#include "linbox/linbox-config.h"

#include <algorithm>
#include <iostream>
#include <vector>

//#include "givaro/givtimer.h"

#include "linbox/ring/modular.h"
#include "linbox/util/commentator.h"
#include "linbox/ring/ntl.h"

#include "linbox/algorithms/invariant-factors.h"
#include "test-frobenius.h"

using namespace LinBox;

#if 0
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
#endif

int main(int argc, char** argv) {
	size_t p = 3;
	size_t k = 0;
	int seed = time(NULL);

	static Argument args[] = {
		{ 'k', "-k K", "Number of invariant factors to compute (default k=0 for all)", TYPE_INT, &k},
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	srand(seed);

   typedef Givaro::Modular<double> Field;
	Field F(p);

   typedef NTL_zz_pX PolyRing;
	PolyRing R(p);

	InvariantFactors<Field, PolyRing> IFD(F, R);

   bool pass = testFrobenius(IFD, F, R, k);
   return pass ? 0 : -1;

#if 0
	size_t b = 4;
	int precond = 0;
	size_t s = 0;
	std::string matrixFile;
	std::string outFile;
	

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 'b', "-b B", "Block size", TYPE_INT, &b},		
		{ 'f', "-f F", "Name of file for matrix", TYPE_STR, &matrixFile},
		{ 'o', "-o O", "Name of output file for invariant factors", TYPE_STR, &outFile},
		{ 'r', "-r R", "Random seed", TYPE_INT, &seed},
		{ 's', "-s S", "Number of nonzeros in random triangular preconditioner", TYPE_INT, &s},
		{ 'c', "-c C", "Choose what preconditioner to apply", TYPE_INT, &precond},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);
	if (matrixFile == "") {
		std::cout << "an input matrix must be provided" << std::endl;
		return -1;
	}
	
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
		time2([&](){return IFD.det(det, FM, b);});
		F.write(std::cout << "det: ", det) << std::endl;
	} else if (precond == 2) { // rank
		size_t rank;
		time2([&](){return IFD.rank(rank, FM, b);});
		std::cout << "rank: " << rank << std::endl;
	} else if (precond == 0) {
		time1([&](){IFD.largestInvariantFactors(result, FM, b);});
	} else if (precond == 4) {
		time1([&](){IFD.largestInvariantFactors2(result, FM, b);});
	}
	
	if (outFile != "") {
		writeLifs(R, outFile, result);
	}
	
	return 0;
#endif
}
