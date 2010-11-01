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



#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <utility>

#include "linbox/blackbox/zo.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/field/modular.h" 
#include "linbox/util/commentator.h"


#include "test-common.h"
#include "test-blackbox.h"

int main(int argc, char **argv) 
{
  bool pass = true;
  uint32 prime = 31337;
  static size_t n = 100000;

  static Argument args[] = 
    {{ 'n', "-n N", "Set dimension of test matrix to NxN.", TYPE_INT, &n }, 
     { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INT, &prime }, 
	 { '\0'}
	};

  parseArguments(argc, argv, args);

  typedef LinBox::Modular<uint32> Field;
  //typedef LinBox::Modular<LinBox::uint32> Field;
  typedef LinBox::ZeroOne<Field> Matrix;

  Field afield(prime);

  /* 
  // "arrow" matrix
  size_t *rows, *cols, i;
  rows = new size_t[3 * n + 1 - 2];
  cols = new size_t[3 * n + 1 - 2];

  for(i = 0; i < n; i++) { rows[i] = 0; cols[i] = i; } // first row
  for(i = 0; i < n - 1; i++) 
    { rows[n+2*i] = i + 1; cols[n+2*i] = 0; rows[n+2*i+1] = i + 1; cols[n+2*i+1] = i + 1; } // first col and the diag
  Matrix testMatrix(afield, rows, cols, n, n, 3 * n - 2);
  */

// random 3 per row matrix
	size_t *rows, *cols, i;
	const size_t npr = n / 10000;
	rows = new size_t[npr * n];
	cols = new size_t[npr * n];

    for(i = 0; i < n; i++)
        {
            set<size_t> a;
            while( a.size() < npr )
                a.insert(rand()%n);
            size_t j = 0;
            for(set<size_t>::iterator iter = a.begin(); j < npr; ++j, ++iter)
                {
                    rows[npr*i+j] = i;
                    cols[npr*i+j] = *iter;
                    //std::cout << rows[npr*i+j] << ", ";
                }
            //std::cout << std::endl;
        }
    ZeroOne<Field> testMatrix(afield, rows, cols, n, n, npr * n );

  /*
  Matrix testMatrix(afield);
  //ifstream mat_in("data/m133.b3.200200x200200.sms");
  ifstream mat_in("data/n4c6.b9.186558x198895.sms");
  //ifstream mat_in("data/small21x21.sms");
  //testMatrix.read(mat_in);
  testMatrix.read(cin);
  //LinBox::Transpose<Matrix> testMat(testMatrix);
  */

  //print out the dimensions and the number of non-zero entries of the matrix
  //std::cout << testMatrix.rowdim() << " " << testMatrix.coldim() << " " << testMatrix.nnz() << std::endl;


  //std::cout << std::endl << "ZeroOne matrix blackbox test suite" << std::endl;

  pass = pass && testBlackbox(testMatrix);
  //bool pass2 = testBlackbox(testMat);
  
  //delete [] rows;
  //delete [] cols;

  //return pass&&pass2 ? 0 : -1;
  return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
