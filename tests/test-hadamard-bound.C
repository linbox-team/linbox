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
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/random-matrix.h"
#include "linbox/solutions/det.h"
#include "linbox/solutions/hadamard-bound.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/commentator.h"

#include <iostream>

using namespace LinBox;

using Field = Givaro::ZRing<Integer>;

template <class TMatrix, class TVector>
bool test_with_matrix_vector(size_t n, size_t bitSize, int* seed)
{
    Field F;
    TMatrix A(F, n, n);

    *seed += 1;
    srand(*seed);
    Field::RandIter randIter(F, bitSize, *seed);

    // Generate a full rank matrix
    RandomDenseMatrix<Field::RandIter, Field> RDM(F, randIter);
    RDM.randomFullRank(A);

    // ---- Determinant

    // Compute the bounds
    auto hb = HadamardBound(A);
    auto fastHb = FastHadamardBound(A);

    // Compute the effective determinant
    Integer detA;
    det(detA, A);

    // std::cout << "det hb fastHb : " << Givaro::logtwo(Givaro::abs(detA)) << " " << hb << " " << fastHb << std::endl;

    if (fastHb + 1 < hb) {
        std::cerr << "Fast Hadamard bound is somehow better than the precise one." << std::endl;
        return false;
    }

    if (Givaro::logtwo(Givaro::abs(detA)) > hb) {
        std::cerr << "The Hadamard bound does not bound the determinant." << std::endl;
        std::cerr << "den: " << Givaro::logtwo(Givaro::abs(detA)) << " > " << hb << std::endl;
        return false;
    }

    // ---- Rational solve

    if (detA > 0) {
        TVector b(F, n);
        b.random(randIter);

        // Compute the bounds
        auto rationalSolveHb = RationalSolveHadamardBound(A, b);

        // Compute the effective solution
        TVector num(F, n);
        Field::Element den;
        solve(num, den, A, b);

        for (size_t i = 0u; i < n; ++i) {
            if (Givaro::logtwo(Givaro::abs(num[i])) > rationalSolveHb.numBoundBitSize) {
                std::cerr << "The rational solve Hadamard bound does not bound the numerator." << std::endl;
                std::cerr << "num[i]: " << Givaro::logtwo(Givaro::abs(num[i])) << " > " << rationalSolveHb.numBoundBitSize
                        << std::endl;
                return false;
            }
        }

        if (Givaro::logtwo(Givaro::abs(den)) > rationalSolveHb.denBoundBitSize) {
            std::cerr << "The rational solve Hadamard bound does not bound the denominator." << std::endl;
            std::cerr << "den: " << Givaro::logtwo(den) << " > " << rationalSolveHb.denBoundBitSize << std::endl;
            return false;
        }
    }

    return true;
}

int main(int argc, char** argv)
{
    size_t n = 10;
    size_t bitSize = 10;
    int iterations = 1;
    int seed = time(NULL);
    bool loop = false;

    static Argument args[] = {{'n', "-n N", "Set dimension of test objects to NxN.", TYPE_INT, &n},
                              {'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations},
                              {'l', "-l", "Infinite testing.", TYPE_BOOL, &loop},
                              {'s', "-s", "Seed for randomness.", TYPE_INT, &seed},
                              {'b', "-b", "Bit size of integers.", TYPE_INT, &bitSize},
                              END_OF_ARGUMENTS};

    // @todo Test rational matrices,
    // and, if so, specialize vectorNormBitSize for it

    parseArguments(argc, argv, args);

    bool ok = true;

    int startingSeed;
    for (auto i = 0; ok && (loop || i < iterations); ++i) {
        startingSeed = seed;
        ok = ok && test_with_matrix_vector<BlasMatrix<Field>, BlasVector<Field>>(n, bitSize, &seed);
        ok = ok && test_with_matrix_vector<SparseMatrix<Field>, BlasVector<Field>>(n, bitSize, &seed);
    }

    if (!ok) {
        std::cerr << "Failed with seed: " << startingSeed << std::endl;
    }

    return ok ? 0 : -1;
}
