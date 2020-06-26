/*
 * examples/ratdet.C
 *
 * Copyright (C) 2020 The LinBox group
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

/**\file examples/ratdet.C
 * @example  examples/ratdet.C
 \brief Determinant of rational matrix
 \ingroup examples
*/

#include <linbox/linbox-config.h>
#include <iostream>

#include <givaro/givrational.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/det.h>

using namespace LinBox;
typedef SparseMatrixFormat::SparseSeq SparseStorage;
typedef Givaro::QField<Givaro::Rational> Rationals;

template<typename RMatrix>
int ratdet(int argc, char **argv) {

    std::ifstream input (argv[1]);
    if (!input) {
        std::cerr << "Error opening matrix file: " << argv[1] << std::endl;
        return -1;
    }

    Rationals QQ;
    MatrixStream<Rationals> ms (QQ, input);
    RMatrix A ( ms );

    Rationals::Element det_A;

    LinBox::Timer tim ; tim.clear() ; tim.start();
    det (det_A, A);
    tim.stop();

    std::cout << "Determinant is ";
    QQ.write(std::cout, det_A) << ':' << std::flush;
    std::clog << tim << std::endl;

    return 0;
}


int main (int argc, char **argv)
{
    commentator().setMaxDetailLevel (-1);
    commentator().setMaxDepth (-1);
    commentator().setReportStream (std::cerr);

    bool dense=true;
    if (argc < 2 || argc > 3) {
        std::cerr << "Usage: ratdet <matrix-file-in-supported-format> [d/s]" << std::endl;
        return -1;
    }

    if (argc == 3) {
        if (argv[2][0] != 'd') dense=false;
    }

    if (dense)
        return ratdet< DenseMatrix<Rationals> >(argc,argv);
    else
        return ratdet< SparseMatrix<Rationals, SparseStorage> >(argc, argv);
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
