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

#include <linbox/matrix/dense-matrix.h>
#include <linbox/solutions/solve.h>

using namespace LinBox;

template <class Ring, class Method>
void run_integer(Communicator& communicator, bool verboseEnabled, size_t dimension) {
    Ring R;

    BlasVector<Ring> b(R, dimension);
    b.setEntry(0, 4);
    b.setEntry(1, 6);

    DenseMatrix<Ring> A(R, dimension, dimension);
    A.setEntry(0, 0, 4);
    A.setEntry(0, 1, 0);
    A.setEntry(1, 0, 0);
    A.setEntry(1, 1, 4);

    Method method;
    method.pCommunicator = &communicator;
    // @fixme Dixon fails with singular with SingularSolutionType::Determinist
    method.singularSolutionType = SingularSolutionType::Diophantine;

    // Rational vector interface will call (num, den) one
    Givaro::QField<Givaro::Rational> F;
    BlasVector<Givaro::QField<Givaro::Rational>> x(F, dimension);


    if (verboseEnabled) {
        std::cout << "--------------- Testing " << Method::name() << " on DenseMatrix<ZZ> (" << dimension << ")" << std::endl;
    }

    try {
        solve(x, A, b, method);
        solveInPlace(x, A, b, method);
    } catch (...) {
        std::cout << "===> OUCH: Throwing error" << std::endl;
        return;
    }

    if (x[0].nume() != 1 || x[0].deno() != 1 || x[1].nume() != 3 || x[1].deno() != 2) {
        A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
        std::cout << "b: " << b << std::endl;
        std::cout << "x: " << x << std::endl;
        std::cerr << "===> OUCH" << std::endl;
    }
}

template <class Matrix, class Method>
void run_modular(bool verboseEnabled) {
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
        std::cout << "--------------- Testing " << Method::name() << " on DenseMatrix<Modular<...>>" << std::endl;
    }

    BlasVector<Field> x(F, 2);

    try {
        solve(x, A, b, method);
        solveInPlace(x, A, b, method);
    } catch (...) {
        std::cout << "===> OUCH: Throwing error" << std::endl;
        return;
    }

    if (x[0] != 1 || x[1] != 52) {
        A.write(std::cout << "A: ", Tag::FileFormat::Maple) << std::endl;
        std::cout << "b: " << b << std::endl;
        std::cout << "x: " << x << std::endl;
        std::cerr << "===> OUCH" << std::endl;
    }
}

int main(int argc, char** argv)
{
    bool verboseEnabled = false;

	static Argument args[] = {
		{ 'v', "-v", "Enable verbose mode.", TYPE_BOOL,     &verboseEnabled },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

    if (verboseEnabled) {
        commentator().setReportStream(std::cout);
    }

    Communicator communicator(0, nullptr);

    run_integer<Givaro::ZRing<Integer>, Method::Auto>(communicator, verboseEnabled, 2);
    // run_integer<Givaro::ZRing<Integer>, Method::CraAuto>(communicator, 2);
    // run_integer<Givaro::ZRing<Integer>, Method::CraAuto>(communicator, 3);
    // run_integer<Givaro::ZRing<Integer>, Method::Dixon>(communicator, 2);
    // run_integer<Givaro::ZRing<Integer>, Method::Dixon>(communicator, 3);
    // run_integer<Givaro::ZRing<Integer>, Method::NumericSymbolicOverlap>(communicator, 2);
    // run_integer<Givaro::ZRing<Integer>, Method::NumericSymbolicOverlap>(communicator, 3); // @fixme Fails
    // run_integer<Givaro::ZRing<Integer>, Method::NumericSymbolicNorm>(communicator, 2); // @fixme Fails
    // run_integer<Givaro::ZRing<Integer>, Method::NumericSymbolicNorm>(communicator, 3); // @fixme Fails

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
