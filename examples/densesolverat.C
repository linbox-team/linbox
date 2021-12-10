/*
 * examples/densesolverat.C
 *
 * Copyright (C) 2012 J-G Dumas
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

/**\file examples/densesolverat.C
 @example examples/densesolverat.C
 @author Jean-Guillaume.Dumas@univ-grenoble-alpes.fr
 * \brief Direct dense solver over the rationals 
 * \ingroup examples
 */
#include <iostream>
#include "givaro/modular.h"

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/methods.h"

using namespace LinBox;
using namespace std;


int main (int argc, char **argv)
{

	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::cerr);


	if (argc < 2 || argc > 4) {
		cerr << "Usage: solve <matrix-file-in-supported-format> [<dense-vector-file>]" << endl;
		return 0;
	}
	std::ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file " << argv[1] << endl; return -1; }
	std::ifstream invect;

	bool createB = false;
	if (argc == 2) {
		createB = true;
	}

	if (argc == 3) {
		invect.open (argv[2], std::ifstream::in);
		if (!invect) {
			createB = true;
		}
		else {
			createB = false;
		}
	}

    {
        typedef Givaro::QField<Givaro::Rational> Rats;
        Rats QQ;
        typedef DenseVector<Rats> RVector;

        MatrixStream<Rats> ms( QQ, input );
        DenseMatrix<Rats> A ( ms );
		std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

		RVector X(QQ, A.coldim()),B(QQ, A.rowdim());

		if (createB) {
			cerr << "Creating a random {-1,1} vector " << endl;
			for(auto it=B.begin(); it != B.end(); ++it)
				if (drand48() <0.5)
					*it = -1;
				else
					*it = 1;
		} else {
			for(auto&& it:B) invect >> it;
            invect.close();
		}

		std::cout << "B is [";
		for(auto it:B) QQ.write(cout, it) << ' ';
		std::cout << ']' << std::endl;

		Timer chrono;

            // DenseElimination dense
		std::cout << "DenseElimination" << std::endl;
		chrono.start();
		solve (X, A, B, Method::DenseElimination());
		chrono.stop();

		std::cout << "(DenseElimination) Solution is [";
		for(auto it:X) QQ.write(cout, it) << " ";
		std::cout << ']' << std::endl;

		std::cout << "CPU time (seconds): " << chrono.usertime() << std::endl;


	}

	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
