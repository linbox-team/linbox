#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
#include <omp.h>

#define LINBOX_USES_OMP 1
#include "linbox/field/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-coppersmith-domain.h"

#include "linbox/solutions/det.h"

// Computes the minimal polynomial of a sparse matrix (given in Matrix Market Format)
// Times BlockCoppersmithDomain using TPL_omp

using namespace LinBox;

typedef Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field, SparseMatrixFormat::TPL_omp> SparseMat;
typedef MatrixDomain<Field> Domain;
typedef typename Domain::OwnMatrix Block;

void benchmarkBCD(Field& F,
                  Domain& MD,
                  SparseMat& M,
                  Block& U,
                  Block& V,
                  std::vector<Block>& gen,
                  std::vector<size_t>& deg,
                  int t)
{
	BlackboxBlockContainer<Field,SparseMat> blockseq(&M,F,U,V);
	BlockCoppersmithDomain<Domain,BlackboxBlockContainer<Field,SparseMat> >
		BCD(MD,&blockseq,t);

	double start=omp_get_wtime();
	deg=BCD.right_minpoly(gen);
	double time=omp_get_wtime()-start;
	std::cout << time << std::endl;
}

int main(int argc, char** argv)
{
	int earlyTerm;
	int p;
	std::string uFname,vFname,mFname;

	static Argument args[] = {
		{ 'p', "-p P", "Set the field GF(p)", TYPE_INT, &p},
		{ 't', "-t T", "Early term threshold", TYPE_INT, &earlyTerm},
		{ 'm', "-m M", "Name of file for matrix M", TYPE_STR, &mFname},
		{ 'u', "-u U", "Name of file for matrix U", TYPE_STR, &uFname},
		{ 'v', "-v V", "Name of file for matrix V", TYPE_STR, &vFname},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	Field F(p);
	Domain MD(F);
	SparseMat M(F);
	Block U(F),V(F);

	{
		ifstream iF(mFname);
		M.read(iF);
		M.finalize();
		iF.close();
	}
	{
		ifstream iF(uFname);
		U.read(iF);
		iF.close();
	}
	{
		ifstream iF(vFname);
		V.read(iF);
		iF.close();
	}

	std::vector<Block> gen;
	std::vector<size_t> deg;

	benchmarkBCD(F,MD,M,U,V,gen,deg,earlyTerm);

	return 0;
}
