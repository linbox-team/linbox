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
#include "linbox/matrix/random-matrix.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/solve.h"
#include "linbox/solutions/hadamard-bound.h"
#include "linbox/util/commentator.h"

#include <iostream>

using namespace LinBox;

using Field = Givaro::ZRing<Integer>;

void test(size_t n)
{
    Field F;
    BlasMatrix<Field> A(F, n, n);

    // Generate a full rank matrix
    Field::RandIter randIter(F, 10); // @fixme bits
    RandomDenseMatrix<Field::RandIter, Field> RDM(F, randIter);
    RDM.randomFullRank(A);

    // ---- Determinant

    // Compute the bounds
    auto hb = HadamardBound(A);
    auto fastHB = HadamardBound(A);

    // Compute the effective determinant
    Integer detA;
    det(detA, A);
    std::cout << "Det: " << Givaro::logtwo(detA) << std::endl;

    if (fastHB < hb) {
        std::cerr << "Fast Hadamard bound is somehow better than the precise one." << std::endl;
        exit(-1);
    }

    if (Givaro::logtwo(detA) > hb) {
        std::cerr << "The Hadamard bound does not bound the determinant." << std::endl;
        exit(-2);
    }

    // ---- Rational solve

    BlasVector<Field> b(F, n);

    // Compute the bounds
    auto rationalSolveHB = RationalSolveHadamardBound(A, b);

    // Compute the effective solution
    BlasVector<Field> num(F, n);
    Field::Element den;
    solve(num, den, A, b);

    std::cout << num[0] << " " << den << std::endl;
    std::cout << Givaro::logtwo(num[0]) << " " << Givaro::logtwo(den) << std::endl;
    std::cout << rationalSolveHB.numBoundBitSize << " " << rationalSolveHB.denBoundBitSize << std::endl;
}

int main(int argc, char** argv)
{
    static size_t n = 10;
    static int iterations = 1;

    // @fixme seed
    // @fixme bitsize

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
