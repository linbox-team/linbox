/** 
 * examples/valence.C
 *
 * Copyright (C) 2005, 2010  J-G Dumas
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */

/**\file examples/valence.C
\brief Valence of sparse matrix over Z or Zp.
\ingroup examples
*/
//#include "linbox-config.h"
#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/valence.h"
#include "linbox/util/matrix-stream.h"

using namespace LinBox;

int main (int argc, char **argv)
{
    commentator.setMaxDetailLevel (-1);
    commentator.setMaxDepth (-1);
    commentator.setReportStream (std::cerr);


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
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
