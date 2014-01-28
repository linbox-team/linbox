/*
 * examples/valence.C
 *
 * Copyright (C) 2005, 2010  J-G Dumas
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/**\file examples/valence.C
 * \brief Valence of sparse matrix over Z or Zp.
 * \ingroup examples
 * @example  examples/valence.C
 */
#include <iostream>

#include <linbox/field/modular.h>
#include <linbox/blackbox/transpose.h>
#include <linbox/blackbox/compose.h>
#include <linbox/matrix/sparse.h>
#include <linbox/solutions/valence.h>
#include <linbox/util/matrix-stream.h>

using namespace LinBox;

int main (int argc, char **argv)
{
	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::cerr);


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
		}
		else {
			Compose< Blackbox, Transpose<Blackbox> > C (&A, &T);
			std::cout << "A A^T is " << C.rowdim() << " by " << C.coldim() << std::endl;
			valence(val_A, C);
		}
	}
	else {
		if (A.rowdim() != A.coldim()) {
			std::cerr << "Valence works only on square matrices, try either to change the dimension in the matrix file, or to compute the valence of A A^T or A^T A, via the -aat or -ata options."  << std::endl;
			exit(0);
		}
		else
			valence (val_A, A);
	}

	std::cout << "Valence is " << val_A << std::endl;

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
