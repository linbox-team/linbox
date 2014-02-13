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

#include <iostream>
#include <sstream>

#include <linbox/field/modular.h>
#include <linbox/field/gf2.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/blackbox/zero-one.h>
#include <linbox/solutions/rank.h>
#include <linbox/util/matrix-stream.h>
#include <linbox/field/givaro.h>

#define SP_STOR SparseMatrixFormat::SparseSeq

using namespace LinBox;

/// rank or rank mod p
int main (int argc, char **argv)
{
	commentator().setMaxDetailLevel (-1);
	commentator().setMaxDepth (-1);
	commentator().setReportStream (std::cerr);

	if (argc < 2 || argc > 3) {
		std::cerr << "Usage: rank <matrix-file-in-supported-format> [<p>]" << std::endl;
		return -1;
	}

	std::ifstream input (argv[1]);
	if (!input) {
		std::cerr << "Error opening matrix file: " << argv[1] << std::endl;
		return -1;
	}

	long unsigned int r;

	LinBox::GivaroRational ZZ;
	MatrixStream<GivaroRational> ms( ZZ, input );
	SparseMatrix<GivaroRational, SP_STOR> A ( ms );
	// SparseMatrix<GivaroRational, SparseMatrixFormat::CSR> A ( ms );
	LinBox::Timer tim ; tim.clear() ; tim.start();
	if (argc == 2) { // rank over the rational numbers.

		/* We could pick a random prime and work mod that prime, But
		 * the point here is that the rank function in solutions/
		 * handles that issue.  Our matrix here is an integer or
		 * rational matrix and our concept is that we are getting the
		 * rank of that matrix by some blackbox magic inside linbox.
		 */
		std::cout << "matrix is " << A.rowdim() << " by " << A.coldim() << std::endl;
		// A.write(std::cout) << std::endl;;

		LinBox::rank (r, A);
	}
	if (argc == 3) { // rank mod a prime
		//for prime greater than wordsize:
		//stringstream qstr(argv[2]);
		//qstr >> q;
		/*
		integer q = atoi(argv[2]);
		typedef Modular<integer> Field;
		*/
		/*
		//to use doubles, prime < 2^{23}
		double q = atof(argv[2]);
		typedef Modular<double> Field;
		*/
		//to use ints, prime < 2^{31}
		int32_t q = atoi(argv[2]);
                if (q == 0) {
                        std::cerr << "second argument should be a non-zero integer or missing\n";
                        return -1;
                }
		typedef Modular<int32_t> Field;

		Field F(q);
/*
		MatrixStream<Field> ms( F, input );
*/
		// SparseMatrix<Field, SparseMatrixFormat::CSR > B (F, A.rowdim(), A.coldim());// modular image of A
		SparseMatrix<Field, SP_STOR > B (F, A.rowdim(), A.coldim());// modular image of A
		MatrixHom::map(B, A);
		std::cout << "matrix is " << B.rowdim() << " by " << B.coldim() << std::endl;
		//if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(std::cout) << std::endl;

		// Using the adaptive LinBox Solution
		LinBox::rank(r,B);

		// using BlackBoxes
#if 0 /*  too bad */
		   Method::Blackbox MBB;
		   MBB.certificate(true);
		   Linbox::rank(r, B, MBB);
#endif

		// using in place Sparse Elimination with linear pivoting

#if 0 /*  too bad */
		   Method::SparseElimination SE;
		   SE.strategy(Specifier::PIVOT_LINEAR);
		   rankin (r, B, SE);
		   if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(std::cout) << std::endl;
#endif


	}
	tim.stop();

	std::cout << "Rank is " << r << " (" << tim << " )" << std::endl;
	return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
