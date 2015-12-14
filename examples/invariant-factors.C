#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
//#include <omp.h>

//#define LINBOX_USES_OMP 1
#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/coppersmith-invariant-factors.h"

// Computes the invariant factors of a sparse matrix (given in Matrix Market Format)
// Effectively times: TPL_omp, BlockCoppersmithDomain and KannanBachem

using namespace LinBox;

typedef Givaro::Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::TPL> SparseMat;
//typedef SparseMatrix<Field, SparseMatrixFormat::TPL_omp> SparseMat;

typedef CoppersmithInvariantFactors<Field,SparseMat> FactorDomain;
typedef typename FactorDomain::PolyDom PolyDom;
typedef typename FactorDomain::PolyRing PolyRing;
typedef DenseVector<PolyRing> FactorVector;

int main(int argc, char** argv)
{
	int earlyTerm;
	int p = 3, b = 3;
	std::string mFname,oFname;

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 't', "-t T", "Early term threshold", TYPE_INT, &earlyTerm},
		{ 'b', "-b B", "Blocking factor", TYPE_INT, &b},
		{ 'm', "-m M", "Name of file for matrix M", TYPE_STR, &mFname},
		{ 'o', "-o O", "Name of file for output", TYPE_STR, &oFname},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	Field F(p);
	SparseMat M(F);

	{
		std::ifstream iF(mFname);
		M.read(iF);
		M.finalize();
		iF.close();
	}

	std::cout << "Finished reading" << std::endl;

	PolyDom PD(F,"x");
	PolyRing R(PD);
	FactorVector factorList(R);
	FactorDomain CIF(F,M,b);

	size_t numFactors=CIF.computeFactors(factorList,earlyTerm);
	std::cout << "Finished computing factors" << std::endl;

	{
		std::ofstream out(oFname);
		for (size_t i=0;i<numFactors;++i) {
			R.write(out,factorList[i]);
			out << std::endl;
		}
		out.close();
	}

	return 0;
}


