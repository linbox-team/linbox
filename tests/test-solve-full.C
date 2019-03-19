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
    struct matrixName;

    template <class MatrixArgs>
    struct matrixName<DenseMatrix<MatrixArgs>> {
        static constexpr const char* value = "DenseMatrix";
    };

    template <class... MatrixArgs>
    struct matrixName<SparseMatrix<MatrixArgs...>> {
        static constexpr const char* value = "SparseMatrix";
    };

    enum : uint8_t {
        TryDense = 0x1,
        TrySparse = 0x2,
        TryBoth = TryDense | TrySparse,
    };
}

template <class SolveMethod, class ResultVector, class Matrix, class Vector>
void print_error(const ResultVector& x, const Matrix& A, const Vector& b, bool verbose, std::string reason)
{
    if (verbose) {
        A.write(std::cerr << "A: ", Tag::FileFormat::Maple) << std::endl;
        std::cerr << "b: " << b << std::endl;
        std::cerr << "x: " << x << std::endl;
    }

    std::cerr << "/!\\ " << SolveMethod::name() << " on " << matrixName<Matrix>::value << " over ";
    A.field().write(std::cerr);
    std::cerr << " of size " << A.rowdim() << "x" << A.coldim() << " FAILS (" << reason << ")" << std::endl;
}

template <class SolveMethod, class Matrix, class Domain, class ResultDomain>
bool test_solve(Domain& D, ResultDomain& RD,Communicator* pCommunicator, bool verbose, size_t m, size_t n)
{
    using Vector = DenseVector<Domain>;
    using ResultVector = DenseVector<ResultDomain>;

    Matrix A(D, m, n);
    Vector b(D, m);
    ResultVector x(RD, n);

    // @note RandomDenseMatrix can also fill a SparseMatrix, it has just a sparsity of 0.
    typename Domain::RandIter randIter(D, 10, 0); // @fixme Set seed, bits
    LinBox::RandomDenseMatrix<typename Domain::RandIter, Domain> RDM(D, randIter);
    b.random();

    // @note Within our test m == n implies invertible,
    // but for the rectangular case, we pass whatever.
    if (m == n) {
        RDM.randomFullRank(A);
    }
    else {
        RDM.random(A);
    }

    if (verbose) {
        std::cout << "--- Testing " << SolveMethod::name() << " on " << matrixName<Matrix>::value << " over ";
        D.write(std::cout) << " of size " << m << "x" << n << std::endl;
    }

    SolveMethod method;
    method.pCommunicator = pCommunicator;

    try {
        solve(x, A, b, method);
        solveInPlace(x, A, b, method);
    } catch (...) {
        print_error<SolveMethod>(x, A, b, verbose, "throws error");
        return false;
    }

    // @fixme CHECK RESULT IS CORRECT
    // if (x[0] != 1 || x[1] != 52) {
    // print_error<SolveMethod>(x, A, b, verbose, "Ax != b");
    // return false;
    // }

    return true;
}

// Testing rational solve over the integers
template <class SolveMethod>
void test_rational_solve(Communicator& communicator, uint8_t tryFlags, bool verbose)
{
    using IntegerDomain = Givaro::ZRing<Integer>;
    using RationalDomain = Givaro::QField<Givaro::Rational>;

    IntegerDomain D;
    RationalDomain RD;

    if (tryFlags & TryDense) {
        // @fixme Forward n, m from arguments
        test_solve<SolveMethod, DenseMatrix<IntegerDomain>>(D, RD, &communicator, verbose, 2, 2);
        test_solve<SolveMethod, DenseMatrix<IntegerDomain>>(D, RD, &communicator, verbose, 2, 3);
    }

    if (tryFlags & TrySparse) {
        test_solve<SolveMethod, SparseMatrix<IntegerDomain>>(D, RD, &communicator, verbose, 2, 2);
        test_solve<SolveMethod, SparseMatrix<IntegerDomain>>(D, RD, &communicator, verbose, 2, 3);
    }
}

// Testing solve over a finite field
template <class SolveMethod>
void test_modular_solve(uint8_t tryFlags, bool verbose)
{
    using ModularDomain = Givaro::Modular<double>;

    ModularDomain D(101); // @fixme Use random or predefined as characteristic

    if (tryFlags & TryDense) {
        test_solve<SolveMethod, DenseMatrix<ModularDomain>>(D, D, nullptr, verbose, 2, 2);
        test_solve<SolveMethod, DenseMatrix<ModularDomain>>(D, D, nullptr, verbose, 2, 3);
    }

    if (tryFlags & TrySparse) {
        test_solve<SolveMethod, SparseMatrix<ModularDomain>>(D, D, nullptr, verbose, 2, 2);
        test_solve<SolveMethod, SparseMatrix<ModularDomain>>(D, D, nullptr, verbose, 2, 3);
    }
}

int main(int argc, char** argv)
{
    bool verbose = false;
    bool loop = false;

    static Argument args[] = {{'v', "-v", "Enable verbose mode.", TYPE_BOOL, &verbose},
                              {'l', "-l", "Infinite loop of tests.", TYPE_BOOL, &loop},
                              END_OF_ARGUMENTS};

    parseArguments(argc, argv, args);

    if (verbose) {
        commentator().setReportStream(std::cout);
    }

    Communicator communicator(0, nullptr);

    do {
        // test_rational_solve<Method::Auto>(communicator, TryBoth, verbose);
        // test_rational_solve<Method::CraAuto>(communicator, TryBoth, verbose);
        // test_rational_solve<Method::Dixon>(communicator, TryBoth, verbose);
        // test_rational_solve<Method::NumericSymbolicOverlap>(communicator, TryDense, verbose); // @fixme Singular case fails
        // test_rational_solve<Method::NumericSymbolicNorm>(communicator, TryDense, verbose);    // @fixme Fails

        test_modular_solve<Method::Auto>(TryBoth, verbose);
        test_modular_solve<Method::DenseElimination>(TryBoth, verbose);
        test_modular_solve<Method::SparseElimination>(TryBoth, verbose);
        // test_modular_solve<Method::Wiedemann>(TryBoth, verbose);
        // test_modular_solve<Method::Lanczos>(TryBoth, verbose);
        // test_modular_solve<Method::BlockLanczos>(TryBoth, verbose);

        // @deprecated These do not compile anymore
        // test_modular_solve<Method::BlockWiedemann>(TryBoth, verbose);
        // test_modular_solve<Method::Coppersmith>(TryBoth, verbose);
    } while (loop);

    return 0;
}
