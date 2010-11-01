/** 
 * examples/det.C
 *
 * Copyright (C) 2005, 2010 D. Saunders,  J-G. Dumas 
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

/**\file examples/det.C examples/det.C
\brief Determinant of sparse matrix over Z or Zp.
\ingroup examples
*/

//#include "linbox-config.h"
#include <iostream>

#include "linbox/field/modular-double.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/blas-blackbox.h"
#include "linbox/solutions/det.h"
#include "linbox/util/matrix-stream.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
    commentator.setMaxDetailLevel (-1);
    commentator.setMaxDepth (-1);
    commentator.setReportStream (std::cerr);


	if (argc <= 1 || argc > 3) {
		cerr << "Usage: det <matrix-file-in-supported-format> [<p>]" << endl;
		return -1;
	}

	if (argc == 2 ) {
	
		// For a small integer matrix test, do "det data/mat2S". 
		// It is a 2 by 2 matrix with determinant = -2.

		typedef PID_integer Integers;		
		Integers ZZ;

		ifstream input (argv[1]);
		if (!input) 
		{ cerr << "Error opening matrix file " << argv[1] << endl; 
			return -1; 
		}
		MatrixStream< Integers> ms ( ZZ, input );
		BlasBlackbox<Integers> A(ms);
		cout << "Matrix is " << A.rowdim() << " by " << A.coldim() << endl;

		Integers::Element det_A;
		det (det_A, A);

		cout << "Determinant is ";
		ZZ.write(cout, det_A) << endl;
	}
	if (argc == 3) { 

		typedef Modular<double> Field;
		double q = atof(argv[2]);
		Field F(q);

		ifstream input (argv[1]);
		if (!input) 
		{ cerr << "Error opening matrix file " << argv[1] << endl; 
		  return -1; 
		}
		MatrixStream< Field > ms ( F, input );
		SparseMatrix<Field> B (ms);
		cout << "Matrix is " << B.rowdim() << " by " << B.coldim() << endl;

		Field::Element det_B;
		det (det_B, B);

		cout << "Determinant is ";
		F.write(cout, det_B) << " mod " << q << endl;
	}

	return 0;
} 
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
