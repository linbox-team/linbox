/*
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

#include <linbox/matrix/random-matrix.h>
#include <linbox/solutions/solve.h>

using namespace LinBox;

namespace {
    template <class Matrix>
    struct TypeToString;

    template <class MatrixArgs>
    struct TypeToString<DenseMatrix<MatrixArgs>> {
        static constexpr const char* value = "DenseMatrix";
    };

    template <class... MatrixArgs>
    struct TypeToString<SparseMatrix<MatrixArgs...>> {
        static constexpr const char* value = "SparseMatrix";
    };

    template <class Matrix>
    inline const char* type_to_string(const Matrix& A)
    {
        return TypeToString<Matrix>::value;
    }

    template <class SolveMethod, class ResultVector, class Matrix, class Vector>
    void print_error(const ResultVector& x, const Matrix& A, const Vector& b, bool verbose, std::string reason)
    {
        if (verbose) {
            A.write(std::cerr << "A: ", Tag::FileFormat::Maple) << std::endl;
            std::cerr << "b: " << b << std::endl;
            std::cerr << "x: " << x << std::endl;
        }

        std::cerr << "/!\\ " << SolveMethod::name() << " on " << type_to_string(A) << " over ";
        A.field().write(std::cerr);
        std::cerr << " of size " << A.rowdim() << "x" << A.coldim() << " FAILS (" << reason << ")" << std::endl;
    }
}

template <class SolveMethod, class Matrix, class Domain, class ResultDomain>
bool test_solve(Domain& D, ResultDomain& RD, Communicator* pCommunicator, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using Vector = DenseVector<Domain>;
    using ResultVector = DenseVector<ResultDomain>;

    //
    // Generating data
    //

    Matrix A(D, m, n);
    Vector b(D, m);
    ResultVector x(RD, n);

    // @note RandomDenseMatrix can also fill a SparseMatrix, it has just a sparsity of 0.
    typename Domain::RandIter randIter(D, bitSize, seed);
    typename Domain::RandIter vectorRandIter(D, vectorBitSize, seed);
    LinBox::RandomDenseMatrix<typename Domain::RandIter, Domain> RDM(D, randIter);

    // @note Within our test m == n implies invertible,
    // but for the rectangular case, we pass whatever.
    if (m == n) {
        RDM.randomFullRank(A);
        if (bitSize != vectorBitSize) {
            b.random(vectorRandIter);
        } else {
            b.random(randIter);
        }
    }
    else {
        RDM.random(A);

        // We set b so that the system is not inconsistent
        Vector x(D, n);
        x.random(randIter);
        A.apply(b, x); // b <- Ax
    }

    if (verbose) {
        std::cout << "--- Testing " << SolveMethod::name() << " on " << type_to_string(A) << " over ";
        D.write(std::cout) << " of size " << m << "x" << n << std::endl;
    }

    //
    // Copying data for final check
    //

    using ResultMatrix = typename Matrix::template rebind<ResultDomain>::other;
    ResultMatrix RA(A, RD);
    ResultVector Rb(RD, b);

    //
    // Solving
    //

    SolveMethod method;
    method.pCommunicator = pCommunicator;

    try {
        solve(x, A, b, method);
        // @fixme Test this result too

        solveInPlace(x, A, b, method);
    } catch (...) {
        print_error<SolveMethod>(x, RA, b, verbose, "throws error");
        return false;
    }

    //
    // Checking result is correct
    //

    ResultVector RAx(RD, m);
    RA.apply(RAx, x);

    VectorDomain<ResultDomain> VD(RD);
    if (!VD.areEqual(RAx, Rb)) {
        print_error<SolveMethod>(x, RA, Rb, verbose, "Ax != b");
        return false;
    }

    return true;
}

//
// Rational solve
//

template <class SolveMethod, class Matrix>
bool test_rational_solve(Communicator& communicator, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using IntegerDomain = Givaro::ZRing<Integer>;
    using RationalDomain = Givaro::QField<Givaro::Rational>;

    IntegerDomain D;
    RationalDomain RD;

    bool ok = true;
    ok &= test_solve<SolveMethod, Matrix>(D, RD, &communicator, m, m, bitSize, vectorBitSize, seed, verbose);
    ok &= test_solve<SolveMethod, Matrix>(D, RD, &communicator, m, n, bitSize, vectorBitSize, seed, verbose);
    return ok;
}

template <class SolveMethod>
bool test_all_rational_solve(Communicator& communicator, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using IntegerDomain = Givaro::ZRing<Integer>;

    bool ok = true;
    ok &= test_rational_solve<SolveMethod, DenseMatrix<IntegerDomain>>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);
    ok &= test_rational_solve<SolveMethod, SparseMatrix<IntegerDomain>>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);
    return ok;
}

template <class SolveMethod>
bool test_dense_rational_solve(Communicator& communicator, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using IntegerDomain = Givaro::ZRing<Integer>;
    return test_rational_solve<SolveMethod, DenseMatrix<IntegerDomain>>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);
}

template <class SolveMethod>
bool test_sparse_rational_solve(Communicator& communicator, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using IntegerDomain = Givaro::ZRing<Integer>;
    return test_rational_solve<SolveMethod, SparseMatrix<IntegerDomain>>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);
}

//
// Modular solve
//

template <class SolveMethod, class Matrix>
bool test_modular_solve(Integer& q, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using ModularDomain = Givaro::Modular<double>;

    ModularDomain D(q);

    bool ok = true;
    ok &= test_solve<SolveMethod, Matrix>(D, D, nullptr, m, m, bitSize, vectorBitSize, seed, verbose);
    ok &= test_solve<SolveMethod, Matrix>(D, D, nullptr, m, n, bitSize, vectorBitSize, seed, verbose);
    return ok;
}

template <class SolveMethod>
bool test_all_modular_solve(Integer& q, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using ModularDomain = Givaro::Modular<double>;

    bool ok = true;
    ok &= test_modular_solve<SolveMethod, DenseMatrix<ModularDomain>>(q, m, n, bitSize, vectorBitSize, seed, verbose);
    ok &= test_modular_solve<SolveMethod, SparseMatrix<ModularDomain>>(q, m, n, bitSize, vectorBitSize, seed, verbose);
    return ok;
}

template <class SolveMethod>
bool test_dense_modular_solve(Integer& q, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using ModularDomain = Givaro::Modular<double>;
    return test_modular_solve<SolveMethod, DenseMatrix<ModularDomain>>(q, m, n, bitSize, vectorBitSize, seed, verbose);
}

template <class SolveMethod>
bool test_sparse_modular_solve(Integer& q, int m, int n, int bitSize, int vectorBitSize, int seed, bool verbose)
{
    using ModularDomain = Givaro::Modular<double>;
    return test_modular_solve<SolveMethod, SparseMatrix<ModularDomain>>(q, m, n, bitSize, vectorBitSize, seed, verbose);
}

int main(int argc, char** argv)
{
    Integer q = 101;
    bool verbose = false;
    bool loop = false;
    int seed = -1;
    int bitSize = 100;
    int vectorBitSize = -1;
    int m = 32;
    int n = 24;

    static Argument args[] = {{'q', "-q", "Field characteristic.", TYPE_INTEGER, &q},
                              {'v', "-v", "Enable verbose mode.", TYPE_BOOL, &verbose},
                              {'l', "-l", "Infinite loop of tests.", TYPE_BOOL, &loop},
                              {'s', "-s", "Seed for randomness.", TYPE_INT, &seed},
                              {'b', "-b", "Bit size for rational solve tests.", TYPE_INT, &bitSize},
                              {'B', "-B", "Vector bit size for rational solve tests (defaults to -b if not specified).", TYPE_INT, &vectorBitSize},
                              {'m', "-m", "Row dimension of matrices.", TYPE_INT, &m},
                              {'n', "-n", "Column dimension of matrices.", TYPE_INT, &n},
                              END_OF_ARGUMENTS};

    if (vectorBitSize < 0) {
        vectorBitSize = bitSize;
    }

    if (seed < 0) {
        seed = time(nullptr);
    }

    parseArguments(argc, argv, args);

    if (verbose) {
        commentator().setReportStream(std::cout);
    }

    bool ok = true;
    Communicator communicator(0, nullptr);

    do {
        ok &= test_all_rational_solve<Method::Auto>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);

        // @bug When bitSize = 5 and vectorBitSize = 50, CRA fails
        ok &= test_all_rational_solve<Method::CRAAuto>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);
        ok &= test_all_rational_solve<Method::Dixon>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);

        // @note SymbolicNumeric methods are only implemented on DenseMatrix
        // @fixme Singular case fails
        // ok &= test_dense_rational_solve<Method::SymbolicNumericOverlap>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);
        // @fixme Fails
        // ok &= test_sparse_rational_solve<Method::SymbolicNumericNorm>(communicator, m, n, bitSize, vectorBitSize, seed, verbose);

        ok &= test_all_modular_solve<Method::Auto>(q, m, n, bitSize, vectorBitSize, seed, verbose);
        ok &= test_all_modular_solve<Method::Auto>(q, m, n, bitSize, vectorBitSize, seed, verbose);
        ok &= test_all_modular_solve<Method::DenseElimination>(q, m, n, bitSize, vectorBitSize, seed, verbose);
        ok &= test_all_modular_solve<Method::SparseElimination>(q, m, n, bitSize, vectorBitSize, seed, verbose);

        // @fixme Dense does not compile
        // ok &= test_dense_modular_solve<Method::Wiedemann>(q, m, n, bitSize, vectorBitSize, seed, verbose);
        ok &= test_sparse_modular_solve<Method::Wiedemann>(q, m, n, bitSize, vectorBitSize, seed, verbose);

        // @fixme Dense is segfaulting
        // ok &= test_dense_modular_solve<Method::Lanczos>(q, m, n, bitSize, vectorBitSize, seed, verbose);
        // @fixme Singular case fails
        // ok &= test_sparse_modular_solve<Method::Lanczos>(q, m, n, bitSize, vectorBitSize, seed, verbose);

        // @fixme Dense does not compile
        // ok &= test_dense_modular_solve<Method::BlockLanczos>(q, m, n, bitSize, vectorBitSize, seed, verbose);
        // @fixme Sparse is segfaulting
        // ok &= test_sparse_modular_solve<Method::BlockLanczos>(q, m, n, bitSize, vectorBitSize, seed, verbose);

        // @deprecated These do not compile anymore
        // ok &= test_all_modular_solve<Method::BlockWiedemann>(q, m, n, bitSize, vectorBitSize, seed, verbose);
        // ok &= test_all_modular_solve<Method::Coppersmith>(q, m, n, bitSize, vectorBitSize, seed, verbose);

        if (!ok) {
            std::cerr << "Failed with seed: " << seed << std::endl;
        }

        seed += 1;
    } while (ok && loop);

    return 0;
}
