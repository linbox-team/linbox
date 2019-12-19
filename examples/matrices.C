/*
 * examples/matrices.C
 *
 * Copyright (C) 2017  D. Saunders, Z. Wang, J-G Dumas
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

/** \file examples/matrices.C
 * @example  examples/matrices.C
 \brief example matrices that were chosen for Smith form testing.
 \ingroup examples

 \author bds & zw

 Various Smith form algorithms may be used for matrices over the
 integers or over Z_m.  Moduli greater than 2^32 are not supported here.
 Several types of example matrices may be constructed or the matrix be read from a file.
 Run the program with no arguments for a synopsis of the command line parameters.

 For the "adaptive" method, the matrix must be over the integers.
 This is expected to work best for large matrices.

 For the "2local" method, the computation is done mod 2^32.

 For the "local" method, the modulus must be a prime power.

 For the "ilio" method, the modulus may be arbitrary composite.
 If the modulus is a multiple of the integer determinant, the integer Smith form is obtained.
 Determinant plus ilio may be best for smaller matrices.

 This example was used during the design process of the adaptive algorithm.
*/

#include <linbox/linbox-config.h>
#include <iostream>
#include <string>

using namespace std;

#include <linbox/util/timer.h>
#include <linbox/matrix/dense-matrix.h>

#include "matrices.h"

template <class PIR>
void Mat(LinBox::DenseMatrix<PIR>& M, PIR& R, int & n, string src) ;

int main(int argc, char* argv[])
{
	if (argc < 3 or argc > 4) {

		cout << "usage: " << argv[0] << " type n [filename]"  << endl;

		cout << " type = `random', `random-rough', `tref',"
             << " 'moler', 'redheffer', "
             << " or `fib',"
			 << " and n is the dimension" << endl;
		cout  << " If filename is present, matrix is written there, else to cout." << endl;

		return 0;
	}

	string type = argv[1];

	int n = atoi(argv[2]);

// place B: Edit here and at place A for ring change
	//typedef PIRModular<int32_t> PIR;
	typedef Givaro::ZRing<Givaro::Integer> PIR;
	PIR R;
	LinBox::DenseMatrix<PIR> M(R,n,n);

	Mat(M,R,n,type);

    if (M.rowdim() <= 20 && M.coldim() <= 20) {
        M.write(std::clog, LinBox::Tag::FileFormat::Maple) << std::endl;
    }

	if (argc == 4) {
		ofstream out(argv[3]);
		M.write(out) << endl;
	} else {
		M.write(cout) << endl;
	}
}// main


/** Output matrix is determined by src which may be:
  "random-rough"
  This mat will have s, near sqrt(n), distinct invariant factors,
  each repeated twice), involving the s primes 101, 103, ...
  "random"
  This mat will have the same nontrivial invariant factors as
  diag(1,2,3,5,8, ... 999, 0, 1, 2, ...).
  "fib"
  This mat will have the same nontrivial invariant factors as
  diag(1,2,3,5,8, ... fib(k)), where k is about sqrt(n).
  The basic matrix is block diagonal with i-th block of order i and
  being a tridiagonal {-1,0,1} matrix whose snf = diag(i-1 1's, fib(i)),
  where fib(1) = 1, fib(2) = 2.  But note that, depending on n,
  the last block may be truncated, thus repeating an earlier fibonacci number.
  "file" (or any other string)
  Also "tref" and file with format "kdense"
  */
template <class PIR>
void Mat(LinBox::DenseMatrix<PIR>& M, PIR& R, int & n,
         string src) {

	if (src == "random-rough") RandomRoughMat(M, R, n);

	else if (src == "random") RandomFromDiagMat(M, R, n);

	else if (src == "fib") RandomFibMat(M, R, n);

	else if (src == "tref") TrefMat(M, R, n);

	else if (src == "krat") KratMat(M, R, n);

	else if (src == "moler") MolerMat(M, R, n);

	else if (src == "redheffer") RedhefferMat(M, R, n);

	else { // from cin, mostly pointless, but may effect a file format change.

			M.read(cin);
		n = M.rowdim();
	}

} // Mat
