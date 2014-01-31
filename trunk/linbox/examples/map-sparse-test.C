#include "examples/map-sparse.h"
#include <tests/test-common.h> //bb: is this supposed to be installed ?

using namespace LinBox;

// why this test here ?
int main(int argc, char* argv[])
{
	typedef Modular<double> Field;
	typedef Field::Element Element;
	int q=65537;
	int n=10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		//{ 's', "-s S", "Sparse matrices with density S.", TYPE_DOUBLE,     &sparsity },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	srand ((unsigned)time (NULL));
	Field F (q);
	MapSparse<Field> ms(F,n,n);
	Element d;
	for (int i=0;i<n;++i) {
		F.init(d,1);
		ms.setEntry(i,i,d);
		F.init(d,2);
		ms.setEntry(n-1-i,i,d);
	}
	F.init(d,2);
	ms.addCol(d,1,3);
	F.init(d,1);
	ms.addCol(d,3,1);
	F.init(d,2);
	ms.timesRow(d,3);
	F.init(d,3);
	ms.timesCol(d,1);
	ms.swapRows(4,5);
	ms.swapCols(1,8);
	if (n <= 30) ms.print(std::cerr);
	std::cerr << ms.nnz() << std::endl;

	ms.randomSim(5*n);
	if (n <= 30) ms.print(std::cerr);
	std::cerr << ms.nnz() << std::endl;
	ms.randomEquiv(n*n/2);
	if (n <= 30) ms.print(std::cerr);
	std::cerr << ms.nnz() << std::endl;

	ms.write(std::cerr);
	return 0;
}
