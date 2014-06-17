#include <iostream>

//#include "linbox/field/modular-int.h"
#include "linbox/field/Givaro/givaro-gfq.h"
#include "linbox/algorithms/gauss.h"

using namespace LinBox;

int main (int argc, char **argv)
{

	if ( argc <  2 || argc > 4) {
		std::cerr << "Usage to get a random null space basis over GF(p,k):  <matrix-file-in-SMS-format> p [k]" << std::endl;
		return -1;
	}

	std::ifstream input (argv[1]);
	if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }

	//typedef Modular<int> Field;
	typedef GivaroGfq Field;
	Field F(atoi(argv[2]),argc>3?atoi(argv[3]):1);
	SparseMatrix<Field, SparseMatrixFormat::SparseSeq > B (F);
	B.read (input);
	std::cout << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;

    BlasMatrix<Field> NullSpace(F,B.coldim(),B.coldim());
    GaussDomain<Field> GD(F);
    
    GD.nullspacebasisin(NullSpace, B);
    
    NullSpace.write( std::cerr << "X:=", Tag::FileFormat::Maple ) << ';' << std::endl;

    std::cerr << "NullsSpace dimensions:" << NullSpace.rowdim() << 'x' << NullSpace.coldim() << std::endl;
    
    return 0;
    

}
