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

/*! @file tests/test-mpi-comm.C
 * @ingroup benchmarks
 * @brief Check MPI communicator interface
 */

#include <givaro/modular.h>
#include <givaro/zring.h>
#include <linbox/linbox-config.h>
#include <linbox/matrix/sparse-matrix.h>

#include "linbox/util/mpicpp.h"

#include "linbox/blackbox/random-matrix.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;

template <class Field, class Matrix>
static bool ensureEqual(const Field& F, Matrix& A, Matrix& A2)
{
    for (auto i = 0u; i < A.rowdim(); ++i) {
        for (auto j = 0u; j < A.coldim(); ++j) {
            if (!F.areEqual(A.getEntry(i, j), A2.getEntry(i, j))) {
                std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
                std::cerr << "               The data communicated is inconsistent                " << std::endl;
                std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
                return false;
            }
        }
    }

    return true;
}

template <class Field>
static bool ensureEqual(const Field& F, BlasVector<Field>& B, BlasVector<Field>& B2)
{
    for (size_t j = 0; j < B.size(); ++j) {
        if (!F.areEqual(B[j], B2[j])) {
            std::cerr << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>" << std::endl;
            std::cerr << "               The data communicated is inconsistent                " << std::endl;
            std::cerr << "<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<" << std::endl;
            return false;
        }
    }
    return true;
}

template <class Field, class Matrix>
void genData(Field& F, Givaro::Integer q, Matrix& A, size_t seed, size_t bits)
{
    typedef typename Field::RandIter RandIter;
    RandIter RI(F, seed);
    LinBox::RandomDenseMatrix<RandIter, Field> RDM(F, RI);
    RDM.random(A);
}

template <class Field>
void genData(Field& F, Givaro::Integer q, BlasVector<Field>& B, size_t seed, size_t bits)
{
    typedef typename Field::RandIter RandIter;
    RandIter RI(F, seed);
    B.random(RI);
}

// 0 ssend B
// 1 recv B as B2
// 0 bcast B
// 1 bcast-recv as B
// 1 checks that B == B2
template <class Field, class Object>
bool test_ssend_recv_bcast(Field& F, Object& B, Object& B2, Communicator& comm)
{
    if (comm.rank() == 0) {
        comm.ssend(B, 1);
    }
    else if (comm.rank() == 1) {
        comm.recv(B2, 0);
    }

    bool ok = false;
    comm.bcast(B, 0);
    if (comm.rank() == 1) {
        ok = ensureEqual(F, B, B2);
    }
    MPI_Bcast(&ok, 1, MPI_CXX_BOOL, 1, MPI_COMM_WORLD);

    return ok;
}

// 0 send B
// 1 recv B as B2
// 1 send B2
// 0 recv B2 as B3
// 0 check that B == B3
template <class Field, class Object>
bool test_send_recv(Field& F, Object& B, Object& B2, Object& B3, Communicator& comm)
{
    if (comm.rank() == 0) {
        comm.send(B, 1);
        comm.recv(B3, 1);
    }
    else if (comm.rank() == 1) {
        comm.recv(B2, 0);
        comm.send(B2, 0);
    }

    bool ok = false;
    if (comm.rank() == 0) {
        ok = ensureEqual(F, B, B3);
    }
    MPI_Bcast(&ok, 1, MPI_CXX_BOOL, 0, MPI_COMM_WORLD);

    return ok;
}

template <class Field>
bool test_with_field(Givaro::Integer q, size_t bits, size_t ni, size_t nj, Communicator& comm, size_t& seed)
{
    Field ZZ(q);

    DenseMatrix<Field> denseMatrix(ZZ, ni, nj), denseMatrix2(ZZ, ni, nj), denseMatrix3(ZZ, ni, nj);
    SparseMatrix<Field> sparseMatrix(ZZ, ni, nj), sparseMatrix2(ZZ, ni, nj), sparseMatrix3(ZZ, ni, nj);
    BlasVector<Field> blasVector(ZZ, nj), blasVector2(ZZ, nj), blasVector3(ZZ, nj);

    // Generating random data for matrice and vector
    if (0 == comm.rank()) {
        genData(ZZ, q, denseMatrix, seed, bits);
        seed += 1;
        genData(ZZ, q, sparseMatrix, seed, bits);
        seed += 1;
        genData(ZZ, q, blasVector, seed, bits);
        seed += 1;
    }

    bool ok = true;
    ok = ok && test_ssend_recv_bcast(ZZ, blasVector, blasVector2, comm);
    ok = ok && test_ssend_recv_bcast(ZZ, denseMatrix, denseMatrix2, comm);
    ok = ok && test_ssend_recv_bcast(ZZ, sparseMatrix, sparseMatrix2, comm);

    ok = ok && test_send_recv(ZZ, blasVector, blasVector2, blasVector3, comm);
    ok = ok && test_send_recv(ZZ, denseMatrix, denseMatrix2, denseMatrix3, comm);
    ok = ok && test_send_recv(ZZ, sparseMatrix, sparseMatrix2, sparseMatrix3, comm);

    return ok;
}

int main(int argc, char** argv)
{
    Communicator comm(&argc, &argv);

    Givaro::Integer q = 101;
    size_t bits = 10;
    size_t niter = 1;
    size_t m = 10;
    size_t n = 10;
    size_t seed = time(nullptr);
    bool loop = false;

    Argument args[] = {{'b', "-b B", "Set the maximum number of digits of integers to generate.", TYPE_INT, &bits},
                       {'i', "-i I", "Set the number of iteration over unit test sets.", TYPE_INT, &niter},
                       {'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL, &loop},
                       {'m', "-m M", "Set column dimension of test matrices to M.", TYPE_INT, &m},
                       {'n', "-n N", "Set row dimension of test matrices to N.", TYPE_INT, &n},
                       {'s', "-s SEED", "Seed used for randomness.", TYPE_INT, &seed},
                       {'q', "-q Q", "Set the field cardinality.", TYPE_INTEGER, &q},
                       END_OF_ARGUMENTS};
    parseArguments(argc, argv, args);

    if (comm.size() < 2) {
        std::cerr << "This test requires at least 2 MPI nodes, but " << comm.size() << " provided." << std::endl;
        std::cerr << "Please run it with mpirun." << std::endl;
        return -2;
    }

    MPI_Bcast(&seed, 1, MPI_INT, 0, MPI_COMM_WORLD);

    bool ok = true;
    for (auto j = 0u; ok && (loop || j < niter); j++) {
        size_t startingSeed = seed;
        srand(seed);

        ok = ok && test_with_field<Givaro::Modular<float>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::Modular<double>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::Modular<int32_t>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::Modular<int64_t>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::Modular<Integer>>(q, bits, m, n, comm, seed);

        ok = ok && test_with_field<Givaro::ZRing<float>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::ZRing<double>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::ZRing<int32_t>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::ZRing<int64_t>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::ZRing<Integer>>(q, bits, m, n, comm, seed);

        ok = ok && test_with_field<Givaro::ModularBalanced<float>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::ModularBalanced<double>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::ModularBalanced<int32_t>>(q, bits, m, n, comm, seed);
        ok = ok && test_with_field<Givaro::ModularBalanced<int64_t>>(q, bits, m, n, comm, seed);

        if (!ok && comm.rank() == 0) {
            std::cerr << "Failed with seed " << startingSeed << std::endl;
            break;
        }
    }

    return ok ? 0 : -1;
}
