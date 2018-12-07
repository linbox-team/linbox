/* Copyright (C) 2018 LinBox
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/solutions/hadamard-bound.h"
#include "linbox/util/commentator.h"
#include "linbox/matrix/random-matrix.h"

#include <iostream>

using namespace LinBox;

using Field = Givaro::ZRing<Integer>;

void test(size_t n)
{
    Field F;
    BlasMatrix<Field> A(F, n, n);

    Field::RandIter randIter(F, 10); // @fixme bits
    RandomDenseMatrix<Field::RandIter, Field> RDM(F, randIter);
    RDM.randomFullRank(A);

    std::cout << HadamardBound(A) << std::endl;
    std::cout << FastHadamardBound(A) << std::endl;

    // @fixme Generate full rank matrices and test that the bound does

    // ---- Rational solve

    BlasVector<Field> b(F, n);

    auto hadamardBound = RationalSolveHadamardBound(A, b);
    std::cout << hadamardBound.numBoundBitSize << " " << hadamardBound.denBoundBitSize << std::endl;
}

int main(int argc, char** argv)
{
    static size_t n = 10;
    static int iterations = 1;

    static Argument args[] = {{'n', "-n N", "Set dimension of test objects to NxN.", TYPE_INT, &n},
                              {'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations},
                              END_OF_ARGUMENTS};

    parseArguments(argc, argv, args);

    bool ok = true;

    for (auto i = 0; i < iterations; ++i) {
        test(n);
    }

    return ok ? 0 : -1;
}
