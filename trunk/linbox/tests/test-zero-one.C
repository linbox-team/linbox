/* Copyright (C) LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#include <iostream>
#include <vector>
#include <utility>

#include "linbox/blackbox/zero-one.h"
#include "linbox/field/modular.h" 

#include "test-common.h"
#include "test-generic.h"

int main(int argc, char **argv) 
{
  bool pass = true;
  uint32 prime = 31337;
  size_t *rows, *cols, i;
  static size_t n = 1000, iter = 1;

  static Argument args[] = {
    { 'n', "-n N", "Set dimension of test matrix to NxN.", TYPE_INT, &n },
    { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &prime },
    { 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iter},
	{ '\0' }
  };

  parseArguments(argc, argv, args);

  typedef LinBox::ZeroOne<LinBox::Modular<LinBox::uint32> > Matrix;

  LinBox::Modular<LinBox::uint32> afield(prime);

  rows = new size_t[3 * n];
  cols = new size_t[3 * n];

  // "arrow" matrix
  for(i = 0; i < n; i++) { rows[i] = 0; cols[i] = i; } // first row
  for(i = 0; i < n; i++) { rows[n+i] = i; cols[n+i] = 0; } // first col
  for(i = 0; i < n; i++) { rows[2*n+i] = i; cols[2*n+i] = i; } // diag
  Matrix testMatrix(afield, rows, cols, n, n, 3*n - 2);

  /*
  for(i = 0; i < n; i++) { rows[i] = i; cols[i] = i; } // diag
  Matrix testMatrix(afield, rows, cols, n, n, n);
  */
/*
  Matrix testMatrix(afield);
  ifstream mat_in("data/n4c6.b9.186558x198895.sms");
  testMatrix.read(mat_in);
  std::cout << testMatrix.rowdim() << " " << testMatrix.coldim() << " " << testMatrix.nnz() << std::endl;
*/


	commentator.start("ZeroOne matrix blackbox test suite", "ZeroOne");

  pass = pass && testBlackbox(testMatrix);
  
  delete [] rows;
  delete [] cols;

	commentator.stop("ZeroOne matrix blackbox test suite");
  return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
