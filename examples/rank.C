/*
 * examples/rank.C
 *
 * Copyright (C) 2005, 2010 D. Saunders, J-G Dumas
 *
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

/**\file examples/rank.C
 * @example  examples/rank.C
  \brief Rank of sparse matrix over Z or Zp.
  \ingroup examples
  */

#include <linbox/linbox-config.h>
#include <iostream>
#include <sstream>
#include <givaro/givrational.h>

#include <linbox/ring/modular.h>
#include <linbox/field/gf2.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/blackbox/zero-one.h>
#include <linbox/solutions/rank.h>
#include <linbox/util/matrix-stream.h>

#define SP_STOR SparseMatrixFormat::SparseSeq

using namespace LinBox;

/// rank or rank mod p
int main (int argc, char **argv)
{
	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::clog);

	if (argc < 2 || argc > 3) {
		std::cerr << "Usage: rank <matrix-file-in-supported-format> [<p>]" << std::endl;
		return -1;
	}

	std::ifstream input (argv[1]);
	if (!input) {
		std::cerr << "Error opening matrix file: " << argv[1] << std::endl;
		return -1;
	}

	size_t r;

	Givaro::QField<Givaro::Rational> QQ;
	LinBox::Timer tim ; tim.clear() ; tim.start();
	MatrixStream<Givaro::QField<Givaro::Rational>> ms (QQ, input);
	SparseMatrix<Givaro::QField<Givaro::Rational>, SP_STOR> A ( ms );

	tim.stop();
	std::clog << "matrix is " << A.rowdim() << " by " << A.coldim() << " (" << tim << ")" << std::endl;
	tim.clear() ; tim.start();

	if (argc == 2) { // rank over the rational numbers.
		/* We could pick a random prime and work mod that prime, But
		 * the point here is that the rank function in solutions/
		 * handles that issue.  Our matrix here is an integer or
		 * rational matrix and our concept is that we are getting the
		 * rank of that matrix by some blackbox magic inside linbox.
		 */
		LinBox::rank (r, A);
	}

	if (argc == 3) { // rank mod a prime
		uint32_t q = atoi(argv[2]);
		if (q == 0) {
				std::cerr << "second argument should be a non-zero integer or missing\n";
				return -1;
		}
		typedef Givaro::Modular<double> Field;

		Field F(q);
		if (q > F.maxCardinality()) {
			std::cerr << "your number is too big for this field" << std::endl;
			return -1 ;
		}

            // modular image of A
		SparseMatrix<Field, SP_STOR > B (F, A.rowdim(), A.coldim());
		MatrixHom::map(B, A);
        tim.stop();
		std::clog << "matrix is " << B.rowdim() << " by " << B.coldim()
                  <<" (time for map: "<< tim << ")" << std::endl;

        tim.clear();tim.start();
//         if (B.rowdim() <= 20 && B.coldim() <= 20)
//             B.write(std::clog) << std::endl;

		// Using the adaptive LinBox Solution
		LinBox::rank(r,B);
	}
	tim.stop();

	std::clog << "Rank is " << std::flush;
	std::cout << r << std::flush;
	std::clog << " (" << tim << " )" << std::endl;
	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
