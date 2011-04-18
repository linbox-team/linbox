/** 
 * examples/rank.C
 *
 * Copyright (C) 2005, 2010 D. Saunders, J-G Dumas
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

/**\file examples/rank.C 
\brief Rank of sparse matrix over Z or Zp.
\ingroup examples
*/


#include <iostream>
#include <sstream>

//#include "linbox/field/modular.h"
#include "linbox/field/modular-double.h"
#include "linbox/field/gf2.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/blackbox/zero-one.h"
#include "linbox/solutions/rank.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/field/givaro-rational.h"

using namespace LinBox;
using namespace std;

/// rank or rank mod p
int main (int argc, char **argv)
{
    commentator.setMaxDetailLevel (-1);
    commentator.setMaxDepth (-1);
    commentator.setReportStream (std::cerr);

	if (argc < 2 || argc > 3) 
	{	cerr << "Usage: rank <matrix-file-in-supported-format> [<p>]" << endl; return -1; }

	ifstream input (argv[1]);
	if (!input) { cerr << "Error opening matrix file: " << argv[1] << endl; return -1; }

	long unsigned int r;

	if (argc == 2) { // rank over the rational numbers.

	   /* We could pick a random prime and work mod that prime, But the point here 
	   is that the rank function in solutions/ handles that issue.  Our matrix here 
	   is an integer or rational matrix and our concept is that we are getting the rank of that 
	   matrix by some blackbox magic inside linbox.
	   */
            LinBox::GivaroRational ZZ;
            MatrixStream<GivaroRational> ms( ZZ, input );
            SparseMatrix<GivaroRational> A ( ms );
		cout << "A is " << A.rowdim() << " by " << A.coldim() << endl;

		rank (r, A);
	}
	if (argc == 3) { // rank mod a prime 
		/* 
		//for prime greater than wordsize:
	    stringstream qstr(argv[2]);
		integer q;
		qstr >> q;
		typedef Modular<integer> Field;
		*/
		/* 
		//to use doubles, prime < 2^{23}
		double q = atof(argv[2]);
		typedef Modular<double> Field;
		*/
		//to use ints, prime < 2^{31}
		LinBox::int32 q = atoi(argv[2]);
		typedef Modular<LinBox::int32> Field;

		Field F(q);
		MatrixStream<Field> ms( F, input );
		SparseMatrix<Field, Vector<Field>::SparseSeq > B (ms);
		cout << "B is " << B.rowdim() << " by " << B.coldim() << endl;
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;

		// Using the adaptive LinBox Solution
		rank(r,B);
                    
		// using BlackBoxes
		/*
		Method::Blackbox MBB;
		MBB.certificate(true);
		rank(r, B, MBB);
		*/

		// using in place Sparse Elimination with linear pivoting
		
		/*
		Method::SparseElimination SE;
		SE.strategy(Specifier::PIVOT_LINEAR);
		rankin (r, B, SE);                    
		if (B.rowdim() <= 20 && B.coldim() <= 20) B.write(cout) << endl;
		*/


	}

	cout << "Rank is " << r << endl;
	return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
