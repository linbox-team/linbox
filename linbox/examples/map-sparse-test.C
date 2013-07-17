#include "examples/map-sparse.h"

using namespace LinBox;

int main()
{
	typedef Modular<double> Field;
	typedef Field::Element Element;
	int q=65537;
	Field F (q);
	MapSparse<Field> ms(F,10,10);
	Element d;
	for (int i=0;i<10;++i) {
		F.init(d,1);
		ms.setEntry(i,i,d);
		F.init(d,2);
		ms.setEntry(9-i,i,d);
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
	ms.print(std::cerr);
	std::cerr << ms.nnz() << std::endl;
	return 0;
}
