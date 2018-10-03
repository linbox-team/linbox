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
#define __LINBOX_HAVE_MPI //! Necessary to compile correctly for the communicator !
#include <givaro/modular.h>
#include <givaro/zring.h>
#include <iostream>
#include <linbox/linbox-config.h>
#include <linbox/matrix/sparse-matrix.h>
#include <stdio.h>
#include <stdlib.h>

#include "linbox/util/mpicpp.h"
#include <mpi.h>
//#include "linbox/util/mpi-gmp.inl"

#include "linbox/blackbox/random-matrix.h"
#include "linbox/matrix/random-matrix.h"

using namespace LinBox;

// @fixme Some of these look like duplicates
template <class Field, class Matrix>
static bool checkResult(const Field& F, Matrix& A, Matrix& A2)
{
    // A.write(std::cout << " A|A3: \n",Tag::FileFormat::Maple) << ';' << std::endl;
    // A2.write(std::cout << " A2|A4: \n",Tag::FileFormat::Maple) << ';' << std::endl;
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
static bool checkResult(const Field& F, BlasVector<Field>& B, BlasVector<Field>& B2)
{
    // B.write(std::cout << " B: \n",Tag::FileFormat::Maple) << ';' << std::endl;
    // B2.write(std::cout << " B2: \n",Tag::FileFormat::Maple) << ';' << std::endl;
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
void genData(Field& F, Givaro::Integer q, Matrix& A, size_t bits)
{
    typedef typename Field::RandIter RandIter;
    RandIter RI(F);
    LinBox::RandomDenseMatrix<RandIter, Field> RDM(F, RI);
    RDM.random(A);
}

template <class Field>
void genData(Field& F, Givaro::Integer q, BlasVector<Field>& B, size_t bits)
{
    typedef typename Field::RandIter RandIter;
    RandIter RI(F);
    B.random(RI);
}

// 0 send B
// 1 recv B as B2
// 0 bcast B
// 1 bcast-recv as B
// 1 checks that B == B2
template <class Field, class Object>
void test_ssend_recv_bcast(Field& F, Object& objectSend, Object& objectRecv, Communicator* Cptr)
{
    if (0 == Cptr->rank()) {
        // double starttime, endtime;
        // MPI_Barrier(MPI_COMM_WORLD);
        // starttime = MPI_Wtime();
        // MPI data distribution for Integer type value
        Cptr->ssend(objectSend, 1);
        // MPI_Barrier(MPI_COMM_WORLD);
        // endtime   = MPI_Wtime();
        // std::cout<<"MPI data distribution used CPU time (seconds): " <<endtime-starttime<<std::endl;
    }
    else if (1 == Cptr->rank()) {
        Cptr->recv(objectRecv, 0);
    }

    Cptr->bcast(objectSend, 0);
    if (1 == Cptr->rank()) {
        checkResult(F, objectSend, objectRecv);
    }
}

template <class Field>
void test_with_field(Givaro::Integer q, size_t bits, size_t ni, size_t nj, Communicator* Cptr)
{
    Field ZZ(q);

    DenseMatrix<Field> A(ZZ, ni, nj), A2(ZZ, ni, nj);
    BlasVector<Field> B2(ZZ, A.coldim()), B(ZZ, A.coldim());
    SparseMatrix<Field> sparseSend(ZZ, ni, nj), sparseRecv(ZZ, ni, nj);

    // Generating random data for matrice and vector
    if (0 == Cptr->rank()) {
        genData(ZZ, q, A, bits);
        genData(ZZ, q, sparseSend, bits);
        genData(ZZ, q, B, bits);
    }

    // @fixme Test also send() (not ssend)
    // Is send useful if we have ssend?

    test_ssend_recv_bcast(ZZ, B, B2, Cptr);
    test_ssend_recv_bcast(ZZ, A, A2, Cptr);
    test_ssend_recv_bcast(ZZ, sparseSend, sparseRecv, Cptr);

    // @fixme Test also
    // 0 send B
    // 1 recv B as B2
    // 1 send B2
    // 0 recv B2 as B3
    // 0 check that B == B3
}

int main(int argc, char** argv)
{
    Communicator Cptr(&argc, &argv);
    size_t bits, niter, n, ni;
    bool loop = false;

    bits = 10, niter = 1, n = 3, ni = 1;
    Givaro::Integer q = 101;
    static Argument args[] = {{'b', "-b B", "Set the maximum number of digits of integers to generate.", TYPE_INT, &bits},
                              {'i', "-i I", "Set the number of iteration over unit test sets.", TYPE_INT, &niter},
                              {'l', "-loop Y/N", "run the test in an infinite loop.", TYPE_BOOL, &loop},
                              {'n', "-n N", "Set column and row dimension of test matrices to N.", TYPE_INT, &n},
                              {'q', "-q Q", "Set the field cardinality.", TYPE_INTEGER, &q},
                              END_OF_ARGUMENTS};
    parseArguments(argc, argv, args);

    if (Cptr.rank() == 0) {
        std::cout << "Connected to " << Cptr.size() << " nodes." << std::endl;
    }

    MPI_Bcast(&n, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&niter, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&loop, 1, MPI::BOOL, 0, MPI_COMM_WORLD);

    bool peak = false;
    srand(time(NULL)); // @fixme Use user-provided seed if any

    for (auto j = 0u; loop || j < niter; j++) {
        if (0 == Cptr.rank()) {
            ni = rand() % n + 1;
            if (ni < n / 2 && ni % 2 == 0 && !peak) ni = 1;

            peak = !peak;
            std::cout << " Test with dimension: " << ni << " x " << ni << std::endl;
        }

        MPI_Bcast(&ni, 1, MPI_INT, 0, MPI_COMM_WORLD);

        // @fixme Why always square? Try rectangular matrices!
        // @fixme Pass reference to communicator, no more pointer
        test_with_field<Givaro::ZRing<Integer>>(q, bits, ni, ni, &Cptr);

        test_with_field<Givaro::Modular<float>>(q, bits, ni, ni, &Cptr);
        test_with_field<Givaro::Modular<double>>(q, bits, ni, ni, &Cptr);
        test_with_field<Givaro::Modular<int32_t>>(q, bits, ni, ni, &Cptr);
        test_with_field<Givaro::Modular<int64_t>>(q, bits, ni, ni, &Cptr);

        test_with_field<Givaro::ZRing<float>>(q, bits, ni, ni, &Cptr);
        test_with_field<Givaro::ZRing<double>>(q, bits, ni, ni, &Cptr);
        test_with_field<Givaro::ZRing<int32_t>>(q, bits, ni, ni, &Cptr);

        test_with_field<Givaro::ModularBalanced<float>>(q, bits, ni, ni, &Cptr);
        test_with_field<Givaro::ModularBalanced<double>>(q, bits, ni, ni, &Cptr);
        test_with_field<Givaro::ModularBalanced<int32_t>>(q, bits, ni, ni, &Cptr);
        test_with_field<Givaro::ModularBalanced<int64_t>>(q, bits, ni, ni, &Cptr);

        test_with_field<Givaro::Modular<Givaro::Integer>>(q, bits, ni, ni, &Cptr);
    }

    return 0;
}
