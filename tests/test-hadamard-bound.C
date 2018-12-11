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
    auto fastHB = FastHadamardBound(A);

    // Compute the effective determinant
    Integer detA;
    det(detA, A);

    if (fastHB < hb) {
        std::cerr << "Fast Hadamard bound is somehow better than the precise one." << std::endl;
        exit(-1);
    }

    if (Givaro::logtwo(Givaro::abs(detA)) > hb) {
        std::cerr << "The Hadamard bound does not bound the determinant." << std::endl;
        exit(-2);
    }

    // ---- Rational solve

    BlasVector<Field> b(F, n);
    b.random();

    // Compute the bounds
    auto rationalSolveHB = RationalSolveHadamardBound(A, b);

    // Compute the effective solution
    BlasVector<Field> num(F, n);
    Field::Element den;
    solve(num, den, A, b);

    for (size_t i = 0u; i < n; ++i) {
        if (Givaro::logtwo(Givaro::abs(num[i])) > rationalSolveHB.numBoundBitSize) {
            std::cerr << "The rational solve Hadamard bound does not bound the numerator." << std::endl;
            std::cout << "num[i]: " << Givaro::logtwo(Givaro::abs(num[i])) << " > " << rationalSolveHB.numBoundBitSize << std::endl;
            exit(-3);
        }
    }

    if (Givaro::logtwo(Givaro::abs(den)) > rationalSolveHB.denBoundBitSize) {
        std::cerr << "The rational solve Hadamard bound does not bound the denominator." << std::endl;
        std::cout << "den: " << Givaro::logtwo(den) << " > " << rationalSolveHB.denBoundBitSize << std::endl;
        exit(-4);
    }
}

int main(int argc, char** argv)
{
    static size_t n = 10;
    static int iterations = 1;

    // @fixme seed
    // @fixme bitsize

    // @fixme print seed on failure

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
