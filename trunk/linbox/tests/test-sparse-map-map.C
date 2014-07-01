
#include "linbox/linbox-config.h"

#include "linbox/field/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/SparseMatrix/sparse-map-map-matrix.h"
#include "linbox/solutions/charpoly.h"
#include "linbox/algorithms/coppersmith-invariant-factors.h"
#include "linbox/vector/blas-vector.h"

#include "tests/test-blackbox.h"

#include <fflas-ffpack/ffpack/ffpack.h>

using namespace LinBox;

typedef Modular<double> Field;
typedef typename Field::Element Element;
typedef SparseMatrix<Field,SparseMatrixFormat::SMM> SparseMat;

#define PRIME 67108859

void makeMat(SparseMat& M)
{
	Field F=M.field();
	int m=M.rowdim(),n=M.coldim();
	Element d;
	F.init(d,0);
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			F.addin(d,F.one);
			M.setEntry(i,j,d);
			if (j==i-1) {
				// test overwriting
				M.setEntry(i,j,F.zero);
			}
		}
	}
}

bool testNNZ(Field& F,int n, int m)
{
	SparseMat M(F,n,m);
	size_t nnz=0;
	if (M.nnz() != 0) {
		return false;
	}
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			M.setEntry(i,j,F.one);
			++nnz;
			if (M.nnz() != nnz) {
				return false;
			}
		}
	}
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			M.setEntry(i,j,F.zero);
			--nnz;
			if (M.nnz() != nnz) {
				return false;
			}
		}
	}
	if (!M.verify()) {
		return false;
	}
	return true;
}

bool testPattern(Field& F, int n, int m)
{
	SparseMat M(F,n,m);
	Element d,e;
	F.init(d,0);
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			F.addin(d,F.one);
			M.setEntry(i,j,d);
			if (j==i-1) {
				M.setEntry(i,j,F.zero);
			}
		}
	}
	F.init(e,0);
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			F.addin(e,F.one);
			M.getEntry(d,i,j);
			if (j==i-1) {
				if (!F.isZero(d)) {
					return false;
				}
			} else if (!F.areEqual(e,d)) {
				return false;
			}
		}
	}

	//Write over old matrix with values offset by 1
	F.init(d,1);
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			F.addin(d,F.one);
			M.setEntry(i,j,d);
			if (j==i) {
				M.setEntry(i,j,F.zero);
			}
		}
	}
	F.init(e,1);
	for (int i=0;i<m;++i) {
		for (int j=0;j<n;++j) {
			F.addin(e,F.one);
			M.getEntry(d,i,j);
			if (j==i) {
				if (!F.isZero(d)) {
					return false;
				}
			} else if (!F.areEqual(e,d)) {
				return false;
			}
		}
	}

	if (!M.verify()) {
		return false;
	}
	return true;
}

bool testSMMBlackbox(Field& F, int n, int m)
{
	SparseMat M(F,m,n);
	makeMat(M);
	return testBlackbox(M,true);
}


bool testSim(Field& F, int n, int b)
{
	SparseMat M(F,n,n);
	Element d;
	F.init(d,0);
	for (int i=0;i<n;++i) {
		F.addin(d,F.one);
		M.setEntry(i,i,d);
	}
	
	typedef CoppersmithInvariantFactors<Field,SparseMat> FactorDomain;
	typedef typename FactorDomain::PolyDom PolyDom;
	typedef typename FactorDomain::PolyRing PolyRing;
	typedef BlasVector<PolyRing> FactorVector;
	PolyDom PD(F,"x");
	PolyRing R(PD);

	FactorVector MFactors(R);
	FactorDomain MCFD(F,M,b);
	int numMFactors=MCFD.computeFactors(MFactors);

	linbox_check(n>=2);
	M.randomSim(2*n);

	FactorVector QFactors(R);
	FactorDomain QCFD(F,M,b);
	int numQFactors=QCFD.computeFactors(QFactors);


	if (numMFactors != numQFactors) {
		return false;
	}

	for (int i=0;i<numQFactors;++i) {
		if (!R.areEqual(QFactors[i],MFactors[i])) {
			return false;
		}
	}

	return true;
}



int main(int argc, char** argv) 
{
	Field F(PRIME);

	static Argument args[] = {
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	bool pass=true;

	pass=pass&&testNNZ(F,10,10);
	pass=pass&&testNNZ(F,1,10);
	pass=pass&&testNNZ(F,2,10);
	pass=pass&&testNNZ(F,10,2);
	pass=pass&&testNNZ(F,10,1);

	pass=pass&&testPattern(F,10,10);
	pass=pass&&testPattern(F,1,10);
	pass=pass&&testPattern(F,2,10);
	pass=pass&&testPattern(F,10,2);
	pass=pass&&testPattern(F,10,1);

	pass=pass&&testSMMBlackbox(F,10,10);
	pass=pass&&testSMMBlackbox(F,1,10);
	pass=pass&&testSMMBlackbox(F,2,10);
	pass=pass&&testSMMBlackbox(F,10,2);
	pass=pass&&testSMMBlackbox(F,10,1);

	pass=pass&&testSim(F,10,5);

	return pass?0:-1;
}

