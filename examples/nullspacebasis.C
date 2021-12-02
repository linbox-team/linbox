/*
 * examples/nullspacebasis.C
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

/**\file examples/nullspacebasis.C
 * @example  examples/nullspacebasis.C
  \brief NullSpace of sparse matrix over GFq.
  \brief nullspace is allocated m \times n.
  \ingroup examples
  */

#include <iostream>
#include "linbox/matrix/dense-matrix.h"
#include <givaro/gfq.h>
#include "linbox/algorithms/gauss.h"

using namespace LinBox;

int main (int argc, char **argv)
{
    if ( argc <  3 || argc > 4) {
        std::cerr << "Usage to get a random null space basis over GF(p,k):  <matrix-file-in-SMS-format> p [k]" << std::endl;
        return -1;
    }
    
    std::ifstream input (argv[1]);
    if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }
    
	//typedef Givaro::Modular<int> Field;
    typedef Givaro::GFqDom<int64_t> Field;
    Field F(atoi(argv[2]),argc>3?atoi(argv[3]):1);
    SparseMatrix<Field, SparseMatrixFormat::SparseSeq > B (F);
    B.read (input);
    std::cout << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;
    
    DenseMatrix<Field> NullSpace(F,B.coldim(),B.coldim());
    GaussDomain<Field> GD(F);
    
    GD.nullspacebasisin(NullSpace, B);
    
    NullSpace.write( std::cerr << "X:=", Tag::FileFormat::Maple ) << ';' << std::endl;

    std::cerr << "NullsSpace dimensions:" << NullSpace.rowdim() << 'x' << NullSpace.coldim() << std::endl;

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
