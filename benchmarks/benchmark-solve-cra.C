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
#include <unistd.h>
#include "givaro/modular.h"
#include "givaro/zring.h"
#include "linbox/linbox-config.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/mpicpp.h"

#include <iostream>
#include <stdio.h>
#include <stdlib.h>

using namespace LinBox;

using Field = Givaro::ZRing<Integer>;

#if defined(__LINBOX_HAVE_MPI)
#define getWTime(...) MPI_Wtime(__VA_ARGS__);
#else
#define getWTime(...) omp_get_wtime(__VA_ARGS__);
#endif

template <class Field, class Matrix>
static bool checkResult(Matrix& A, BlasVector<Field>& B, BlasVector<Field>& X, Integer& d)
{
    Field ZZ;
    BlasVector<Field> B2(ZZ, A.coldim());
    BlasVector<Field> B3(ZZ, A.coldim());
    A.apply(B2, X);
   
    for (size_t j = 0; j < B.size(); ++j) {
        ZZ.mul(B3[j],B[j],d);
    }

    for (size_t j = 0; j < A.coldim(); ++j) {
        if (!ZZ.isZero(B2[j]-B3[j])) {
            std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
            std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
            return false;
        }
    }
    return true;
}

template <class Field, class Matrix, class RandIterator>
void genData(Field& F, RandIterator& RandIter, Matrix& A,  BlasVector<Field>& B)
{

    LinBox::RandomDenseMatrix<typename Field::RandIter, Field> RDM(F, RandIter);
    RDM.randomFullRank(A);
    B.random(RandIter);
}

bool benchmark(size_t niter, BlasVector<Field>& x, BlasMatrix<Field>& A, BlasVector<Field>& B, Communicator* Cptr)
{
    Field::Element d;

    double startTime = getWTime();
    solveCRA(x, d, A, B, RingCategories::IntegerTag(), Method::DenseElimination(), Cptr);

    bool ok = false;
    if (Cptr->master()) {
        double endTime = getWTime();
        std::cout << "CPU time (seconds): " << (endTime - startTime)/double(niter) << std::endl;
        ok = checkResult(A, B, x, d);
    }

#ifdef __LINBOX_HAVE_MPI
    MPI_Bcast(&ok, 1, MPI_C_BOOL, 0, MPI_COMM_WORLD);
#endif
    return ok;
}

int main(int argc, char** argv)
{
    Communicator communicator(&argc, &argv);

    if (communicator.master()) {
        std::cout << "Communicator size: " << communicator.size() << std::endl;
    }

    int seed = -1;
    size_t bits = 10;
    size_t niter = 1;
    size_t n = 100;

    static Argument args[] = {{'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT, &n},
                              {'b', "-b B", "Set the maximum number of digits of integers to generate.", TYPE_INT, &bits},
                              {'i', "-i I", "Set the number of times to do the random unit tests.", TYPE_INT, &niter},
                              {'s', "-s SEED", "Set the seed for randomness (random if negative).", TYPE_INT, &seed},
                              END_OF_ARGUMENTS};
    parseArguments(argc, argv, args);

    if (seed < 0) {
        seed = time(NULL);
    }
    srand(seed);

/*    
#ifdef __LINBOX_HAVE_MPI
    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&bits, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&niter, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
*/

    Field ZZ;
    DenseMatrix<Field> A(ZZ, n, n);
    BlasVector<Field> b(ZZ, A.coldim());
    BlasVector<Field> x(ZZ, A.coldim());
    
    Field::RandIter RI(ZZ, bits, seed);
    
    bool ok = true;
    for (size_t j = 0u; j < niter; j++) {

        if (communicator.master()) {
            genData(ZZ, RI, A, b);
        }

#ifdef __LINBOX_HAVE_MPI
        communicator.bcast(A, 0);
        communicator.bcast(b, 0);
#endif

        ok = benchmark(niter, x, A, b, &communicator);
        if (!ok) break;

        ++seed;
    }

    if (!ok) {
        if (communicator.master()) {
            std::cerr << "Failed with seed: " << seed << std::endl;
        }
        return 1;
    }

    return 0;
}
