/*
 * examples/nullspacebasis_rank.C
 *
 * Copyright (C) 2014 J-G. Dumas
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

/**\file examples/nullspacebasis_rank.C
 * @example  examples/nullspacebasis_rank.C
  \brief NullSpace of sparse matrix over GFq.
  \brief nullspace is allocated rank \times n.
  \ingroup examples
  */

#include <iostream>
#include "linbox/matrix/dense-matrix.h"
#include "linbox/field/Givaro/givaro-gfq.h"
#include "linbox/algorithms/gauss.h"

using namespace LinBox;

int main (int argc, char **argv)
{
  if ( argc <  2 || argc > 4) {
    std::cerr << "Usage to get a random null space basis over GF(p,k):  <matrix-file-in-SMS-format> p [k]" << std::endl;
    return -1;
  }

  std::ifstream input (argv[1]);
  if (!input) {
    std::cerr << "Error opening matrix file " << argv[1] << std::endl;
    return -1;
  }
  int pVal = atoi(argv[2]);
  int kVal = atoi(argv[3]);


  //typedef Modular<int> Field;
  typedef GivaroGfq Field;
  Field F(pVal, argc>3?kVal:1);
  SparseMatrix<Field, SparseMatrixFormat::SparseSeq > A (F);
  A.read (input);
  std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

  GaussDomain<Field> GD(F);
  
  typename Field::Element Det;
  unsigned long Rank;
  size_t Ni(A.rowdim()),Nj(A.coldim());

  Permutation<Field> P((int)Nj,F);

  GD.InPlaceLinearPivoting(Rank, Det, A, P, Ni, Nj );

  for(size_t i=0; i< Ni; ++i) {
    if (A[i].size() == 0) {
      size_t j(i);
      if (nextnonzero(j,Ni,A)) {
	A[i] = A[j];
	A[j].resize(0);
      }
      else {
	break;
      }
    }
  }
  size_t nullity = A.coldim()-Rank;
  BlasMatrix<Field> NullSpace(F,A.coldim(),nullity);
  GD.nullspacebasis(NullSpace, Rank, A, P);

  NullSpace.write( std::cerr << "X:=", Tag::FileFormat::Maple ) << ';' << std::endl;
    
  std::cerr << "NullsSpace dimensions:" << NullSpace.rowdim() << 'x' << NullSpace.coldim() << std::endl;
    
  return 0;
}
