
#include "linbox/linbox-config.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/solutions/rank.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;


bool testPrimeField(int p, int n, int m)
{
	typedef Modular<double> Field;
	typedef MatrixDomain<Field> Domain;
	typedef typename Domain::OwnMatrix Matrix;
	typedef typename Field::RandIter RandIter;

	Field F(p);
	Matrix M(F,n,m);
	RandIter RI(F);
	RandomDenseMatrix<RandIter,Field> RDM(F,RI);
	long unsigned r;

	RDM.randomFullRank(M);
	rank(r,M,Method::BlasElimination());

	return r==(long unsigned)(n<m?n:m);
}

int main(int argc, char** argv)
{
	bool pass=true;

	static Argument args[] = {
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	pass=pass&& testPrimeField(101,10,10);
	pass=pass&& testPrimeField(101,2,10);
	pass=pass&& testPrimeField(101,10,2);
	pass=pass&& testPrimeField(101,1,10);

	return pass?0:-1;
}
