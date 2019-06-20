/*
 * examples/sparseelimrank.C
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

/** \file examples/sparseelimrank.C
 * @example  examples/sparseelimrank.C
 \brief Gaussian elimination Rank of sparse matrix over Z or Zp.
 \ingroup examples
 */

#include <linbox/linbox-config.h>

#include <iostream>
#include <vector>
#include <utility>
#include <givaro/zring.h>
#include <givaro/givrational.h>

#include <givaro/modular.h>
#include <linbox/field/gf2.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/blackbox/zero-one.h>
#include <linbox/solutions/rank.h>
#include <linbox/util/matrix-stream.h>

using namespace LinBox;
using namespace std;

int main (int argc, char **argv)
{
	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::cerr);
    commentator().setBriefReportStream (std::cout);

	if (argc < 2 || argc > 3)
	{	cerr << "Usage: rank <matrix-file-in-supported-format> [<p>]" << endl; return -1; }

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

	size_t r;

	if (argc == 2) { // rank over the integers.

		/* We could pick a random prime and work mod that prime, But the point here
		   is that the rank function in solutions/ handles that issue.  Our matrix here
		   is an integer matrix and our concept is that we are getting the rank of that
		   matrix by some blackbox magic inside linbox.
		   */
		Givaro::QField<Givaro::Rational> QQ;
		MatrixStream<Givaro::QField<Givaro::Rational>> ms( QQ, input );
		SparseMatrix<Givaro::QField<Givaro::Rational>, SparseMatrixFormat::SparseSeq > A ( ms );
        if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::cerr << "A:=",Tag::FileFormat::Maple) << ';' << std::endl;


		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		LinBox::rank (r, A, Method::SparseElimination() );

		cout << "Z Rank is " << r << endl;
	}
	if (argc == 3) {
		double q = atof(argv[2]);
		typedef Givaro::Modular<double> Field;
		Field F(q);
		MatrixStream<Field> ms( F, input );
		SparseMatrix<Field, SparseMatrixFormat::SparseSeq > B (ms);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;


		Method::SparseElimination SE;
// 		SE.pivotStrategy = PivotStrategy::None;
// 		// using Sparse Elimination
// 		LinBox::rank (r, B, SE);
// 		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;
// 		cout << "Rank is " << r << endl;

		SE.pivotStrategy = PivotStrategy::Linear;
		// using Sparse Elimination
        Givaro::Timer chrono; chrono.start();
		LinBox::rankInPlace (r, B, SE);
        chrono.stop();
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;
		F.write(cout << "Rank is " << r << " over ") << endl;
        std::cerr << chrono << std::endl;

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
