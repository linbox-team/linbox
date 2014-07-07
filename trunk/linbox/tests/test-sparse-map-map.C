
#include "linbox/linbox-config.h"
#include "linbox/util/commentator.h"

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

template <class Field>
class SMMTests {
public:
	typedef typename Field::Element Element;
	typedef SparseMatrix<Field,SparseMatrixFormat::SMM> SparseMat;

	static bool runBasicTests(Field& F, int n, int m)
	{
		bool pass=true;
		stringstream ss;
		ss << "Running basic tests with n=" << n << " m=" << m;
		commentator().start(ss.str().c_str(),0,3);
		pass=pass&&testNNZ(F,n,m);
		pass=pass&&testPattern(F,n,m);
		pass=pass&&testSMMBlackbox(F,n,m);
		commentator().stop(MSG_STATUS(pass));
		return pass;
	}

	static void setMatPattern(SparseMat& M, const Element& startVal)
	{
		Field F=M.field();
		int m=M.rowdim(),n=M.coldim();
		Element d;
		F.assign(d,startVal);
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

	static bool testMatPattern(SparseMat& M, const Element& startVal)
	{
		Field F=M.field();
		int m=M.rowdim(),n=M.coldim();
		Element d,e;
		F.assign(d,startVal);
		for (int i=0;i<m;++i) {
			for (int j=0;j<n;++j) {
				F.addin(d,F.one);
				M.getEntry(e,i,j);
				if (j==i-1) {
					if (!F.isZero(e)) {
						return false;
					}
				} else {
					if (!F.areEqual(d,e)) {
						return false;
					}
				}
			}
		}
		return true;
	}

	static bool testNNZ(Field& F,int n, int m)
	{
		commentator().start("Testing nnz()");
		SparseMat M(F,n,m);
		size_t nnz=0;
		if (M.nnz() != 0) {
			commentator().stop(MSG_FAILED);
			return false;
		}
		for (int i=0;i<m;++i) {
			for (int j=0;j<n;++j) {
				M.setEntry(i,j,F.one);
				++nnz;
				if (M.nnz() != nnz) {
					commentator().stop(MSG_FAILED);
					return false;
				}
			}
		}
		for (int i=0;i<m;++i) {
			for (int j=0;j<n;++j) {
				M.setEntry(i,j,F.zero);
				--nnz;
				if (M.nnz() != nnz) {
					commentator().stop(MSG_FAILED);
					return false;
				}
			}
		}
		if (!M.verify()) {
			commentator().stop(MSG_FAILED);
			return false;
		}
		commentator().stop(MSG_PASSED);
		return true;
	}

	static bool testPattern(Field& F, int n, int m)
	{
		commentator().start("Testing matrix overwrites");
		SparseMat M(F,n,m);
		bool pass=true;

		setMatPattern(M,F.zero);
		pass=pass&&testMatPattern(M,F.zero);
		pass=pass&&M.verify();

		//Write over old matrix with values offset by 1
		setMatPattern(M,F.one);
		pass=pass&&testMatPattern(M,F.one);
		pass=pass&&M.verify();

		commentator().stop(MSG_STATUS(pass));
		return pass;
	}

	static bool testSMMBlackbox(Field& F, int n, int m)
	{
		commentator().start("Running testBlackbox()");
		SparseMat M(F,m,n);
		setMatPattern(M,F.zero);
		bool pass=testBlackbox(M,true);
		commentator().stop(MSG_STATUS(pass));
		return pass;
	}

	static bool testSim(Field& F, int n, int b)
	{
		stringstream ss;
		ss << "Testing randomSim() with n=" << n << " b=" << b;
		commentator().start(ss.str().c_str());

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
			commentator().stop(MSG_FAILED);
			return false;
		}

		for (int i=0;i<numQFactors;++i) {
			if (!R.areEqual(QFactors[i],MFactors[i])) {
				commentator().stop(MSG_FAILED);
				return false;
			}
		}
		commentator().stop(MSG_PASSED);
		return true;
	}
};

template <class Field>
bool basicTestSuite(Field& F)
{
	bool pass=true;
	commentator().start("Basic tests",0,5);
	pass=pass&&SMMTests<Field>::runBasicTests(F,10,10);
	pass=pass&&SMMTests<Field>::runBasicTests(F,1,10);
	pass=pass&&SMMTests<Field>::runBasicTests(F,2,10);
	pass=pass&&SMMTests<Field>::runBasicTests(F,10,2);
	pass=pass&&SMMTests<Field>::runBasicTests(F,10,1);
	commentator().stop(MSG_STATUS(pass));
	return pass;
}

template <class Field>
bool matGenTestSuite(Field& F)
{
	bool pass=true;
	commentator().start("Testing randomSim()");
	pass=pass&&SMMTests<Field>::testSim(F,10,10);
	commentator().stop(MSG_STATUS(pass));
	return pass;
}

template <class Field>
bool fieldTestSuite(Field& F)
{
	bool pass=true;
	stringstream ss;
	ss << "Testing SMM with Field ";
	F.write(ss);
	commentator().start(ss.str().c_str());
	pass=pass&&basicTestSuite(F);
	pass=pass&&matGenTestSuite(F);
	commentator().stop(MSG_STATUS(pass));
	return pass;
}

bool testSuite()
{
	bool pass=true;
	Modular<double> DoubleSmallField(3),DoubleLargeField(67108859);
	pass=pass&&fieldTestSuite(DoubleSmallField);
	pass=pass&&fieldTestSuite(DoubleLargeField);

	GivaroExtension<Modular<double> > ExtField(DoubleSmallField,3);
	pass=pass&&fieldTestSuite(ExtField);
	return pass;
}

int main(int argc, char** argv) 
{
	static Argument args[] = {
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	commentator().start("Testing sparse-map-map-matrix");

	bool pass=testSuite();

	commentator().stop(MSG_STATUS(pass));

	return pass?0:-1;
}



