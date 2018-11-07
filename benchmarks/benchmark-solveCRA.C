/* Copyright (C) 2018 The LinBox group
 * Written by Hongguang Zhu <zhuhongguang2014@gmail.com>
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file tests/test-solveCRA.C
 * @ingroup benchmarks
 * @brief Testing the MPI parallel/serial rational solver
 */

#define __LINBOX_HAVE_MPI
#define __Detailed_Time_Measurement

#include "givaro/modular.h"
#include "givaro/zring.h"
#include "linbox/linbox-config.h"
#include "linbox/matrix/sparse-matrix.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>

#include "linbox/matrix/random-matrix.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/mpicpp.h"

using namespace LinBox;

using Field = Givaro::ZRing<Integer>;

template <class Field, class Matrix>
static bool checkResult(const Field& ZZ, Matrix& A, BlasVector<Field>& B, BlasVector<Field>& X, Integer& d)
{
    BlasVector<Field> B2(ZZ, A.coldim());
    BlasVector<Field> B3(ZZ, A.coldim());
    A.apply(B2, X);

    Integer tmp;
    for (size_t j = 0; j < B.size(); ++j) {
        B3.setEntry(j, d * B.getEntry(j));
    }
    for (size_t j = 0; j < A.coldim(); ++j) {
        if (!ZZ.areEqual(B2[j], B3[j])) {
            std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
            std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
            return false;
        }
    }
    return true;
}

template <class Field, class Matrix>
void genData(Field& F, Matrix& A, size_t bits)
{
    typedef typename Field::RandIter RandIter;
    RandIter RI(F, bits, 5);
    LinBox::RandomDenseMatrix<RandIter, Field> RDM(F, RI);
    RDM.randomFullRank(A);
}

template <class Field>
void genData(Field& F, BlasVector<Field>& B, size_t bits)
{
    typedef typename Field::RandIter RandIter;
    RandIter RI(F, bits, 5);
    B.random(RI);
}

void run(BlasVector<Field>& X2, BlasMatrix<Field>& A, BlasVector<Field>& B, Communicator& communicator)
{
    Field ZZ;
    Field::Element d;

    double starttime = MPI_Wtime();
    solveCRA(X2, d, A, B, RingCategories::IntegerTag(),
             Method::BlasElimination(),
             // Method::Hybrid(communicator),
             &communicator);

    if (0 == communicator.rank()) {
        double endtime = MPI_Wtime();
        std::cout << "Total CPU time (seconds): " << endtime - starttime << std::endl;
        checkResult(ZZ, A, B, X2, d);
    }

    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char** argv)
{
    // @fixme To be able to not fill this file with __LINBOX_HAVE_MPI,
    // we should just have a dummy communicator at user level,
    // which wouldn't broadcast anything.
    Communicator communicator(&argc, &argv);

    size_t bits = 10;
    size_t n = 1;

    static Argument args[] = {{'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT, &n},
                              {'b', "-b B", "Set the mxaimum number of digits of integers to generate.", TYPE_INT, &bits},
                              END_OF_ARGUMENTS};
    parseArguments(argc, argv, args);

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);

    Field ZZ;
    DenseMatrix<Field> A(ZZ, n, n);
    BlasVector<Field> X(ZZ, A.coldim()), B(ZZ, A.coldim());

    // Generating data
    if (0 == communicator.rank()) {
        genData(ZZ, A, bits);
        genData(ZZ, B, bits);
    }

    communicator.bcast(A, 0);
    communicator.bcast(B, 0);

    run(X, A, B, communicator);

    return EXIT_SUCCESS;
}
