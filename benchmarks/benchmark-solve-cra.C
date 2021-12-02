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
#include <unistd.h>

using namespace LinBox;

using Ring = Givaro::ZRing<Integer>;

#if defined(__LINBOX_HAVE_MPI)
#define getWTime(...) MPI_Wtime(__VA_ARGS__)
#define Time2Seconds(start,end) static_cast<double>(end-start)
#elif defined(__LINBOX_HAVE_OMP)
#define getWTime(...) omp_get_wtime(__VA_ARGS__)
#define Time2Seconds(start,end) static_cast<double>(end-start)
#else
#include <chrono>
#define getWTime(...) std::chrono::high_resolution_clock::now()
#define Time2Seconds(start,end) std::chrono::duration_cast<std::chrono::duration<double>>(end-start).count()
#endif

template <class Matrix>
static bool checkResult(Matrix& A, BlasVector<Ring>& B, BlasVector<Ring>& X, Integer& d)
{
    const Ring& ZZ = A.field();
    BlasVector<Ring> B2(ZZ, A.coldim());
    BlasVector<Ring> B3(ZZ, A.coldim());
    A.apply(B2, X);

    VectorDomain<Ring> VD(ZZ);
    VD.mul(B3, B, d);

    if (!VD.areEqual(B2, B3)) {
        std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
        std::cerr << "               The solution of solveCRA is incorrect                " << std::endl;
        std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
        return false;
    }

    return true;
}

template <class Matrix, class RandIterator>
void genData(Ring& ZZ, RandIterator& randIter, Matrix& A, BlasVector<Ring>& B)
{

    LinBox::RandomDenseMatrix<typename Ring::RandIter, Ring> RDM(ZZ, randIter);
    RDM.randomFullRank(A);
    B.random(randIter);
}

bool benchmark(size_t niter, BlasVector<Ring>& x, BlasMatrix<Ring>& A, BlasVector<Ring>& B, Communicator& communicator)
{
    Ring::Element d;

    auto startTime = getWTime();
    Method::CRAAuto method;
    method.pCommunicator = &communicator;
    solve(x, d, A, B, method);

    bool ok = false;
    if (communicator.master()) {
        auto endTime = getWTime();
        std::cout << "CPU time (seconds): " << Time2Seconds(startTime,endTime) / double(niter) << std::endl;
        ok = checkResult(A, B, x, d);
    }

    communicator.bcast(ok, 0);
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
    bool loop = false;

    static Argument args[] = {{'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT, &n},
                              {'b', "-b B", "Set the maximum number of digits of integers to generate.", TYPE_INT, &bits},
                              {'i', "-i I", "Set the number of times to do the random unit tests.", TYPE_INT, &niter},
                              {'s', "-s SEED", "Set the seed for randomness (random if negative).", TYPE_INT, &seed},
                              {'l', "-l", "Infinite loop (ignoring -i).", TYPE_BOOL, &loop},
                              END_OF_ARGUMENTS};
    parseArguments(argc, argv, args);

    // As only master uses the seed, there is no need to broadcast it
    if (communicator.master()) {
        if (seed < 0) {
            seed = time(NULL);
        }
        srand(seed);
    }

    Ring ZZ;
    DenseMatrix<Ring> A(ZZ, n, n);
    BlasVector<Ring> b(ZZ, A.coldim());
    BlasVector<Ring> x(ZZ, A.coldim());

    bool ok = true;
    for (size_t j = 0u; loop | (j < niter); j++) {
        if (communicator.master()) {
            Ring::RandIter randIter(ZZ, seed);
            randIter.setBitsize(bits);
            genData(ZZ, randIter, A, b);
        }

        communicator.bcast(A, 0);
        communicator.bcast(b, 0);

        ok = benchmark(niter, x, A, b, communicator);
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
