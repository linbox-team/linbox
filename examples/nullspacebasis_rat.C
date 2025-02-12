/*
 * examples/nullspacebasis_rat.C
 *
 * Copyright (C) The LinBox group
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
#include <givaro/givrational.h>
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/algorithms/gauss.h"

using namespace LinBox;

typedef Givaro::QField<Givaro::Rational> QRat;
typedef MatrixStream<QRat> QMstream;
typedef SparseMatrix<QRat, SparseMatrixFormat::SparseSeq > QSparseMat;

int main (int argc, char **argv)
{
    if ((argc != 2) && (argc != 3)) {
        std::cerr << "Usage to get a null space basis over Q:  <matrix-file-in-SMS-format> [output-file]" << std::endl;
        return -1;
    }

    std::ifstream input (argv[1]);
    if (!input) { std::cerr << "Error opening matrix file " << argv[1] << std::endl; return -1; }

    QRat QQ;
    QMstream ms(QQ, input);
    QSparseMat B (ms);
    std::cout << "B is " << B.rowdim() << " by " << B.coldim() << std::endl;

    SparseMatrix<QRat> NullSpace(QQ,B.coldim(),B.coldim());
    GaussDomain<QRat> GD(QQ);

    GD.nullspacebasisin(NullSpace, B);

    NullSpace.write( std::clog << "X:=", Tag::FileFormat::Maple ) << ';' << std::endl;

    if (argc>2) {
        std::ofstream output(argv[2]);
        NullSpace.write( output, Tag::FileFormat::SMS );
        std::clog << "Nullspace basis written in " << argv[2] << std::endl;
    }

    std::clog << "NullSpace dimensions:" << NullSpace.rowdim() << 'x' << NullSpace.coldim() << std::endl;

    return 0;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
