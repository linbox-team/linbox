/**
* Copyright (C) LinBox
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

/**
 * This is testing the serialization routines,
 * by serializing, then unserializing and checking equality.
 *
 * Basic types (integer and floating points) are checked first.
 * Custom LinBox classes (BlasMatrix, SparseMatrix, ...) are checked too.
 */

#include "linbox/matrix/random-matrix.h"
#include "linbox/util/serialization.h"

using namespace LinBox;

template <class T>
bool test(T& output, const T& input)
{
    // The vector of bytes that will store the serialized input
    // is created with a certain size of uninitialized data.
    // This is used to check that the unserialized does
    // understand the offset start.
    uint64_t randomOffset = rand() % 100;
    std::vector<uint8_t> bytes(randomOffset);

    auto bytesWritten = serialize(bytes, input);
    if (bytesWritten == 0u) {
        return false;
    }

    auto bytesRead = unserialize(output, bytes, randomOffset);
    if (bytesRead != bytesWritten) {
        return false;
    }

    return true;
}

template <class T>
bool check_basic_type(const T& input)
{
    T output;
    if (!test(output, input)) {
        return false;
    }

    if (output != input) {
        return false;
    }

    return true;
}

template <class T>
bool test_basic_type()
{
    // @note This little trick is to ensure that 64 bits type
    // have not just their lower bits written.
    T input = rand();
    input *= rand();
    input += rand();

    return check_basic_type(input);
}

bool test_integer()
{
    // Generating a big integer to ensure it is spreading over multiple limbs.
    Integer input(1);
    for (int i = 0; i < 20; ++i) {
        input *= (1 + rand());
    }
    input += rand();

    return check_basic_type(input);
}

// Check matrix equality after serialize/unserialize.
template <class Field, class Matrix>
bool check_matrix(const Field& F, const Matrix& input)
{
    Matrix output(F);

    if (!test(output, input)) {
        return false;
    }

    if (output.rowdim() != input.rowdim() || output.coldim() != input.coldim()) {
        return false;
    }

    for (auto i = 0u; i < input.rowdim(); i++) {
        for (auto j = 0u; j < input.coldim(); j++) {
            if (!F.areEqual(input.getEntry(i, j), output.getEntry(i, j))) {
                return false;
            }
        }
    }

    return true;
}

// Check vector equality after serialize/unserialize.
template <class Field, class Vector>
bool check_vector(const Field& F, const Vector& input)
{
    Vector output(F);

    if (!test(output, input)) {
        return false;
    }

    if (output.size() != input.size()) {
        return false;
    }

    for (auto i = 0u; i < input.size(); i++) {
        if (!F.areEqual(input[i], output[i])) {
            return false;
        }
    }

    return true;
}

// Tests serialibility of matrices and vectors of specified field elements.
template <class Field>
bool test_field(const Integer& q)
{
    // --- Test dense matrix

    Field F(q);
    BlasMatrix<Field> denseMatrix(F, 10 + rand() % 100, 10 + rand() % 100);

    // Fill denseMatrix with random values!
    typename Field::RandIter R(F);
    RandomDenseMatrix<typename Field::RandIter, Field> RandMat(F, R);
    RandMat.random(denseMatrix);

    check_matrix(F, denseMatrix);

    // --- Test sparse matrix

    SparseMatrix<Field> sparseMatrix(F, denseMatrix.rowdim(), denseMatrix.coldim());
    for (auto i = 0u; i < denseMatrix.rowdim(); i++) {
        for (auto j = 0u; j < denseMatrix.coldim(); j++) {
            if (rand() % 2 == 0) {
                sparseMatrix.setEntry(i, j, denseMatrix.getEntry(i, j));
            }
        }
    }

    check_matrix(F, sparseMatrix);

    // --- Test dense vector

    BlasVector<Field> denseVector(F, denseMatrix.rowdim());
    for (auto i = 0u; i < denseMatrix.rowdim(); i++) {
        denseVector[i] = denseMatrix.getEntry(i, 0);
    }

    check_vector(F, denseVector);

    return true;
}

int main(int argc, char** argv)
{
    Integer q = 101;
    uint64_t seed = time(nullptr);
    bool loop = false;

    Argument as[] = {{'q', "-q Q", "Set the field characteristic (-1 for random).", TYPE_INTEGER, &q},
                     {'s', "-s seed", "Set seed for the random generator", TYPE_UINT64, &seed},
                     {'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL, &loop},
                     END_OF_ARGUMENTS};

    FFLAS::parseArguments(argc, argv, as);

    srand(seed);

    bool ok = true;
    do {
        ok = ok && test_basic_type<int8_t>();
        ok = ok && test_basic_type<uint8_t>();
        ok = ok && test_basic_type<int16_t>();
        ok = ok && test_basic_type<uint16_t>();
        ok = ok && test_basic_type<int32_t>();
        ok = ok && test_basic_type<uint32_t>();
        ok = ok && test_basic_type<int64_t>();
        ok = ok && test_basic_type<uint64_t>();
        ok = ok && test_basic_type<float>();
        ok = ok && test_basic_type<double>();

        ok = ok && test_integer();

        ok = ok && test_field<Givaro::ZRing<Integer>>(q);
        ok = ok && test_field<Givaro::Modular<float>>(q);
        ok = ok && test_field<Givaro::Modular<double>>(q);
    } while (loop && ok);

    if (!ok) std::cerr << "Failed with seed: " << seed << std::endl;

    return !ok;
}
