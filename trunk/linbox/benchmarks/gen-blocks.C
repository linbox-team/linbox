#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <omp.h>
#include <time.h>
#include <set>
#include <utility>
#include <sstream>

#include "linbox/field/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-coppersmith-domain.h"

#include "linbox/solutions/det.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/methods.h"

#include "linbox/algorithms/wiedemann.h"

#include "examples/map-sparse.h"

// Generates random dense matrices U and V and sparse matrix M for use
// with block-coppersmith-benchmark and invariant-factors-benchmark.
// Also computes the min-poly and saves that as well

using namespace LinBox;

typedef Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field> SparseMat;
typedef MatrixDomain<Field> Domain;
typedef typename Domain::OwnMatrix Block;

int randRange(int start, int end)
{
        double rval = rand();
        static const double NORMALIZING_CONSTANT = 1.0/(1.0+RAND_MAX);
        double normedRVal = rval*NORMALIZING_CONSTANT;
        double rangeSize = end-start;
        int offset = rangeSize*normedRVal;
        return start+offset;
}

void randomBlock(Block& block, Field& field, int q, int m, int n)
{
  Element d;
  for (int i=0;i<m;++i) {
    for (int j=0;j<n;++j) {
      field.init(d,randRange(0,q));
      block.setEntry(i,j,d);
    }
  }
}

void randomNonSingular(Block& block, Field& field, int q, int m, int n)
{
  randomBlock(block,field,q,m,n);
  long unsigned r;
  long unsigned b=(long unsigned)(m<n?m:n);
  std::cerr << "starting" << std::endl;
  while (LinBox::rank(r,block,Method::BlasElimination()) < b) {
    randomBlock(block,field,q,m,n);
    std::cerr << "loop" << std::endl;
  }
}

std::string fileDesc(int n, int b, double s, int q)
{
  std::stringstream ss;
  ss << "n" << n;
  ss << "b" << b;
  int sInt=10000*(s-floor(s));
  ss << "s" << setfill('0') << setw(4) << sInt << setw(0);
  ss << "q" << q;
  ss << ".txt";

  return ss.str();
}

void computeMinPoly(Domain& MD,
		    SparseMat& M,
		    Field& F,
		    Block& U,
		    Block& V,
		    std::vector<Block>& gen,
		    std::vector<size_t>& deg,
		    int t)
{
  BlackboxBlockContainer<Field,SparseMat> blockseq(&M,F,U,V);
  BlockCoppersmithDomain<Domain,BlackboxBlockContainer<Field,SparseMat> > BCD(MD,&blockseq,t);

  deg=BCD.right_minpoly(gen);
}

int main(int argc, char** argv) {
	int seed=501;
	int n;
	int b;
	int nnz;
	double sparsity;
	int t;
	int q=3;

	srand(seed);

	static Argument args[] = {
		{ 'n', "-n N", "Set row/col dimension of test matrix n.", TYPE_INT, &n },
		{ 'b', "-b B", "Set block size b.", TYPE_INT, &b },
		{ 's', "-s S", "Set the sparsity [0-1).", TYPE_DOUBLE, &sparsity },
		{ 'q', "-q Q", "Set the field GF(p)", TYPE_INT, &q},
		{ 't', "-t T", "Early term threshold", TYPE_INT, &t},
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	nnz=(int)((double)(n*n)*sparsity);

	Field F(q);
	Domain MD(F);
	SparseMat M(F,n,n);
	Block U(F,b,n),V(F,n,b);

	std::vector<Block> gen;
	std::vector<size_t> deg;

	randomNonSingular(U,F,q,b,n);
	randomNonSingular(V,F,q,n,b);

	std::cerr << "done1" << std::endl;
	MapSparse<Field> sparse(F,n,n);
	MapSparse<Field>::generateSparseNonSingular(sparse,nnz,seed);
	sparse.copy(M);
	std::cerr << "done2" << std::endl;

	{
	  ofstream oF("U"+fileDesc(n,b,sparsity,q));
	  U.write(oF);
	  oF.close();
	}

	{
	  ofstream oF("V"+fileDesc(n,b,sparsity,q));
	  V.write(oF);
	  oF.close();
	}

	{
	  ofstream oF("M"+fileDesc(n,b,sparsity,q));
	  M.write(oF);
	  oF.close();
	}

	std::cerr << "done3" << std::endl;
	computeMinPoly(MD,M,F,U,V,gen,deg,t);

	{
	  ofstream oF("MP"+fileDesc(n,b,sparsity,q));
	  for (int i=0;i<deg.size()-1;++i) {
	    oF << deg[i] << " ";
	  }
	  if (deg.size() > 0) {
	    oF << deg[deg.size()-1];
	  }
	  oF << std::endl;

	  for (int i=0;i<gen.size();++i) {
	    oF << gen[i] << std::endl;
	  }
	}

	std::cout << fileDesc(n,b,sparsity,q) << std::endl;
	return 0;
}
