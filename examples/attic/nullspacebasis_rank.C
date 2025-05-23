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

#include <givaro/gfq.h>
#include <iostream>
#include "linbox/matrix/dense-matrix.h"
#include "linbox/algorithms/gauss.h"

using namespace LinBox;

int main (int argc, char **argv)
{
  if ( argc < 3 || argc > 4) {
    std::cerr << "Usage to get a random null space basis over GF(p,k):  <matrix-file-in-SMS-format> p [k]" << std::endl;
    return -1;
  }

  std::ifstream input (argv[1]);
  if (!input) {
    std::cerr << "Error opening matrix file " << argv[1] << std::endl;
    return -1;
  }
  int pVal = atoi(argv[2]);


  //typedef Givaro::Modular<int> Field;
  typedef Givaro::GFqDom<int64_t> Field;
  Field F(pVal, argc>3?atoi(argv[3]):1);
  SparseMatrix<Field, SparseMatrixFormat::SparseSeq > A (F);
  A.read (input);
  std::cout << "A is " << A.rowdim() << " by " << A.coldim() << std::endl;

  GaussDomain<Field> GD(F);
  
  typename Field::Element Det;
  size_t Rank;
  size_t Ni(A.rowdim()),Nj(A.coldim());

  Permutation<Field> P(F,(int)Nj);

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
  DenseMatrix<Field> NullSpace(F,A.coldim(),nullity);
  Givaro::Timer chrono; chrono.start();
  GD.nullspacebasis(NullSpace, Rank, A, P);
  chrono.stop();

  NullSpace.write( std::cout << "X:=", Tag::FileFormat::Maple ) << ';' << std::endl;
    
  std::clog << "NullSpace dimensions:" << NullSpace.rowdim() << 'x' << NullSpace.coldim()  << ' '
            << chrono << std::endl;
    
  return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
