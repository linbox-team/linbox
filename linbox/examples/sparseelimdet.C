/*
 * examples/sparseelimdet.C
 *
 * Copyright (C) 2006, 2010  J-G Dumas
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

/** \file examples/sparseelimdet.C
 * @example  examples/sparseelimdet.C
 \brief Gaussian elimination determinant of sparse matrix over Z or Zp.
 \ingroup examples
 */
//#include "linbox-config.h"

#include <iostream>
#include <vector>
#include <utility>

#include "linbox/field/modular.h"
#include "linbox/field/gf2.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/zero-one.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/det.h"
#include "linbox/util/matrix-stream.h"

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
	commentator.setMaxDetailLevel (-1);
	commentator.setMaxDepth (-1);
	commentator.setReportStream (std::cerr);

	if (argc < 2 || argc > 3)
	{	cerr << "Usage: sparseelimdet <matrix-file-in-supported-format> [<p>]" << endl; return -1; }

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

	Method::SparseElimination SE;

	if (argc == 2) { // determinant over the integers.

		PID_integer ZZ;
		PID_integer::Element d;
		MatrixStream<PID_integer> ms( ZZ, input );
		SparseMatrix<PID_integer, Vector<PID_integer>::SparseSeq > A ( ms );
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		SE.strategy(Specifier::PIVOT_LINEAR);
		det (d, A, SE);

		ZZ.write(cout << "Determinant is ", d) << endl;
	}
	if (argc == 3) {
		double q = atof(argv[2]);
		typedef Modular<double> Field;
		Field::Element d;
		Field F(q);
		MatrixStream<Field> ms( F, input );
		SparseMatrix<Field, Vector<Field>::SparseSeq > B (ms);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;


		SE.strategy(Specifier::PIVOT_NONE);
		// using Sparse Elimination
		det (d, B, SE);
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;
		F.write(cout << "Determinant is ", d) << endl;

		SE.strategy(Specifier::PIVOT_LINEAR);
		// using Sparse Elimination with reordering
		detin (d, B, SE);
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;
		F.write(cout << "Determinant is ", d) << endl;


	}

	return 0;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

