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
	SparseMatrix<Field, SparseMatrixFormat::SparseSeq > A (F);
	A.read (input);
	std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

    GaussDomain<Field> GD(F);
    
    typename Field::Element Det;
    unsigned long Rank;
    size_t Ni(A.rowdim()),Nj(A.coldim());

    Permutation<Field> P((int)Nj,F);

    GD.InPlaceLinearPivoting(Rank, Det, A, P, Ni, Nj );

    for(size_t i=0; i< Ni; ++i) {
        if (A[i].size() == 0) {
            size_t j(i);
            if (nextnonzero(j,Ni,A)) {
                A[i] = A[j];
                A[j].resize(0);
            }
            else {
                break;
            }
        }
    }
    size_t nullity = A.coldim()-Rank;
    BlasMatrix<Field> NullSpace(F,A.coldim(),nullity);
	GD.nullspacebasis(NullSpace, Rank, A, P);

    NullSpace.write( std::cerr << "X:=", Tag::FileFormat::Maple ) << ';' << std::endl;
    
    std::cerr << "NullsSpace dimensions:" << NullSpace.rowdim() << 'x' << NullSpace.coldim() << std::endl;
    
    return 0;
    

}
