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

#include <linbox/solutions/solve.h>

using namespace LinBox;

template <class SolveMethod, class ResultVector, class Matrix, class Vector>
void print_error(const ResultVector& x, const Matrix& A, const Vector& b, bool verboseEnabled, std::string reason)
{
    if (verboseEnabled) {
        A.write(std::cerr << "A: ", Tag::FileFormat::Maple) << std::endl;
        std::cerr << "b: " << b << std::endl;
        std::cerr << "x: " << x << std::endl;
    }

    // @note Within our tests A square implies non-singular
    bool singular = (A.rowdim() != A.coldim());
    std::cerr << "==> " << SolveMethod::name() << " on " << (singular ? "" : "non-") << "singular DenseMatrix over ";
    A.field().write(std::cerr);
    std::cerr << " FAILS (" << reason << ")" << std::endl;
}

template <class Ring, class SolveMethod>
void run_integer(Communicator& communicator, bool verboseEnabled, size_t dimension)
{

    // @fixme Test non-singular

    // @fixme Test singular (rectangular)

    Ring R;

    BlasVector<Ring> b(R, dimension);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    DenseMatrix<Ring> A(R, dimension, dimension);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    SolveMethod method;
    method.pCommunicator = &communicator;

    // Rational vector interface will call (num, den) one
    using RationalDomain = Givaro::QField<Givaro::Rational>;
    RationalDomain F;
    BlasVector<RationalDomain> x(F, dimension);

    if (verboseEnabled) {
        A.field().write(std::cout << "--------------- Testing " << SolveMethod::name() << " on DenseMatrix over ") << std::endl;
    }

    try {
        solve(x, A, b, method);
        solveInPlace(x, A, b, method);
    } catch (...) {
        print_error<SolveMethod>(x, A, b, verboseEnabled, "throws error");
        return;
    }

    if (x[0].nume() != 1 || x[0].deno() != 1 || x[1].nume() != 3 || x[1].deno() != 2) {
        print_error<SolveMethod>(x, A, b, verboseEnabled, "Ax != b");
    }
}

template <class Matrix, class Method>
void run_modular(bool verboseEnabled)
{
    using Field = typename Matrix::Field;

    Field F(101);

    BlasVector<Field> b(F, 2);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    Matrix A(F, 2, 2);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    Method method;
    method.blockingFactor = 1;

    if (verboseEnabled) {
        A.field().write(std::cout << "--------------- Testing " << Method::name() << " on DenseMatrix over ") << std::endl;
    }

    BlasVector<Field> x(F, 2);

    try {
        solve(x, A, b, method);
        solveInPlace(x, A, b, method);
    } catch (...) {
        print_error<SolveMethod>(x, A, b, verboseEnabled, "throws error");
        return;
    }

    if (x[0] != 1 || x[1] != 52) {
        print_error<SolveMethod>(x, A, b, verboseEnabled, "Ax != b");
    }
}

// Testing rational solve over the integers
template <class SolveMethod>
void test_rational_solve(Communicator& communicator, bool verboseEnabled)
{
    //
    // Testing invertible matrix
    //

    run_integer<Givaro::ZRing<Integer>, SolveMethod>(communicator, verboseEnabled, 2);


    //
    // Testing singular matrix
    //

    run_integer<Givaro::ZRing<Integer>, SolveMethod>(communicator, verboseEnabled, 3);
}

int main(int argc, char** argv)
{
    bool verboseEnabled = false;

    static Argument args[] = {{'v', "-v", "Enable verbose mode.", TYPE_BOOL, &verboseEnabled}, END_OF_ARGUMENTS};

    parseArguments(argc, argv, args);

    if (verboseEnabled) {
        commentator().setReportStream(std::cout);
    }

    Communicator communicator(0, nullptr);
    test_rational_solve<Method::Auto>(communicator, verboseEnabled);
    test_rational_solve<Method::CraAuto>(communicator, verboseEnabled);
    test_rational_solve<Method::Dixon>(communicator, verboseEnabled);
    test_rational_solve<Method::NumericSymbolicOverlap>(communicator, verboseEnabled); // @fixme Singular case fails
    test_rational_solve<Method::NumericSymbolicNorm>(communicator, verboseEnabled);    // @fixme Fails

    run_modular<DenseMatrix<Givaro::Modular<double>>, Method::Auto>(verboseEnabled);
    run_modular<SparseMatrix<Givaro::Modular<double>>, Method::Auto>(verboseEnabled);
    // run_modular<DenseMatrix<Givaro::Modular<double>>, Method::DenseElimination>();
    // run_modular<SparseMatrix<Givaro::Modular<double>>, Method::DenseElimination>();
    // run_modular<DenseMatrix<Givaro::Modular<double>>, Method::SparseElimination>();
    // run_modular<SparseMatrix<Givaro::Modular<double>>, Method::SparseElimination>();
    // // run_modular<DenseMatrix<Givaro::Modular<double>>, Method::Wiedemann>(); // @fixme Can't compile
    run_modular<SparseMatrix<Givaro::Modular<double>>, Method::Wiedemann>(verboseEnabled);
    // run_modular<DenseMatrix<Givaro::Modular<double>>, Method::Lanczos>(); // @fixme Segmentation fault
    // run_modular<SparseMatrix<Givaro::Modular<double>>, Method::Lanczos>(); // @fixme Segmentation fault
    // run_modular<DenseMatrix<Givaro::Modular<double>>, Method::BlockLanczos>(); // @fixme Can't compile
    // run_modular<SparseMatrix<Givaro::Modular<double>>, Method::BlockLanczos>(); // @fixme Segmentation fault

    // @deprecated These do not compile anymore
    // run_modular<DenseMatrix<Givaro::Modular<double>>, Method::BlockWiedemann>();
    // run_modular<DenseMatrix<Givaro::Modular<double>>, Method::Coppersmith>();

    // @fixme Test rectangular matrices...

    return 0;
}
