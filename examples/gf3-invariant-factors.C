#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
#include <omp.h>

#define LINBOX_USES_OMP 1


#include "linbox/blackbox/pascal.h"

#include "linbox/ring/modular.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/sliced3.h"
#include "linbox/algorithms/coppersmith-invariant-factors.h"

using namespace LinBox;

template<class Field>
std::istream& readVector(std::istream& is,
                         const Field& F,
                         std::vector<typename Field::Element>& out)
{
	char c;
	is >> std::ws;
	is >> c;
	linbox_check(c=='[');
	is >> std::ws;
	out.clear();
	c=is.peek();
	while (is.good() && (c != ']')) {
		typename Field::Element d;
		F.read(is,d);
		out.push_back(d);
		is >> std::ws;
		c=is.peek();
	}
	return is;
}

typedef SlicedField<Givaro::Modular<int64_t>,uint64_t> Field;
typedef typename Field::Element Element;
typedef PascalBlackbox<Field> SparseMat;

typedef CoppersmithInvariantFactors<Field,SparseMat,Givaro::Modular<int64_t> > FactorDomain;
typedef typename FactorDomain::PolyDom PolyDom;
typedef typename FactorDomain::PolyRing PolyRing;
typedef DenseVector<PolyRing> FactorVector;

int main(int argc, char** argv)
{
	int earlyTerm;
	int b;
	std::string mFname,oFname;

	static Argument args[] = {
		{ 't', "-t T", "Early term threshold", TYPE_INT, &earlyTerm},
		{ 'b', "-b B", "Blocking factor", TYPE_INT, &b},
		{ 'm', "-m M", "Name of file for coefficients", TYPE_STR, &mFname},
		{ 'o', "-o O", "Name of file for output", TYPE_STR, &oFname},
		END_OF_ARGUMENTS
	};

	parseArguments(argc,argv,args);

	Field F(3);
	std::vector<Element> coeffs;
	{
		std::ifstream iF(mFname);
		readVector(iF,F,coeffs);
		iF.close();
	}
	int n=coeffs.size();
	coeffs.resize(2*n);
	SparseMat M(n,n,coeffs,F);

#if 0
	{
		MatrixDomain<Field> MD;
		typename MatrixDomain<Field>::OwnMatrix I(F,n,n),O(F,n,n);
		for (int i=0;i<n;++i)
			for (int j=0;j<n;++j)
				I.setEntry(i,j,i==j?1:0);
		M.applyLeft(O,I);
		for (int i=0;i<n;++i) {
			for (int j=0;j<n;++j) {
				Element d;
				O.getEntry(d,i,j);
				std::cout << (int)d;
			}
			std::cout << std::endl;
		}
	}
#else
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
#endif
	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
