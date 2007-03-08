/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
/**\file examples/valence.C
\brief Valence of sparse matrix over Z or Zp.
\ingroup examples
*/
//#include "linbox-config.h"
#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/field/gmp-integers.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/valence.h"
#include "linbox/util/matrix-stream.h"

using namespace LinBox;

int main (int argc, char **argv)
{

	if (argc < 2 || argc > 3) {
		std::cerr << "Usage: valence <matrix-file-in-supported-format> [-ata|-aat]" << std::endl;
		return -1;
	}

	std::ifstream input (argv[1]);
	if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }

        PID_integer ZZ;
        MatrixStream< PID_integer > ms( ZZ, input );
        typedef SparseMatrix<PID_integer>  Blackbox;
        Blackbox A (ms);
        
        std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;
        
        PID_integer::Element val_A;


        if (argc == 3) {
            Transpose<Blackbox> T(&A);
            if (strcmp(argv[2],"-aat")) {
                Compose< Transpose<Blackbox>, Blackbox > C (&T, &A);
                std::cout << "A^T A is " << C.rowdim() << " by " << C.coldim() << std::endl;
                valence(val_A, C);
            } else {
                Compose< Blackbox, Transpose<Blackbox> > C (&A, &T);
                std::cout << "A A^T is " << C.rowdim() << " by " << C.coldim() << std::endl;
                valence(val_A, C);                
            }
        } else {
            if (A.rowdim() != A.coldim()) {
                std::cerr << "Valence works only on square matrices, try either to change the dimension in the matrix file, or to compute the valence of A A^T or A^T A, via the -aat or -ata options."  << std::endl;
                exit(0);
            } else
                valence (val_A, A);
        }
        
        std::cout << "Valence is " << val_A << std::endl;
        
	return 0;
}
