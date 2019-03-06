/*
 * Copyright(C) LinBox
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

#pragma once

#include <linbox/algorithms/diophantine-solver.h>
#include <linbox/algorithms/rational-solver.h> // @todo Rename dixon-rational-solver?
#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/matrix/sparse-matrix.h>

namespace LinBox {
    namespace {
        // @fixme Remove when RationalSolver is renamed
        template <class... Args>
        using DixonRationalSolver = RationalSolver<Args...>;

        template <class Matrix>
        struct MethodForMatrix {
        };

        template <class Ring>
        struct MethodForMatrix<DenseMatrix<Ring>> {
            using type = Method::Dixon;
        };

        template <class Ring>
        struct MethodForMatrix<SparseMatrix<Ring>> {
            using type = Method::SparseElimination;
        };
    }

    /**
     * \brief Solve specialisation for Dixon on dense matrices.
     */
    template <class Matrix, class Vector>
    void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b,
               const RingCategories::IntegerTag& tag, const MethodWIP::Dixon& m)
    {
        commentator().start("solve.dixon.integer.dense");
        linbox_check((A.coldim() != xNum.size()) || (A.rowdim() != b.size()));

        using Ring = typename Matrix::Field;
        using Field = Givaro::Modular<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()));

        // @note The template argument Method::Dixon means that we want to use a dense algorithm.
        using Solver = DixonRationalSolver<Ring, Field, PrimeGenerator, typename MethodForMatrix<Matrix>::type>;
        Solver dixonSolve(A.field(), primeGenerator);

        // Either A is known to be non-singular, or we just don't know yet.
        int maxTrials = m.trialsBeforeFailure;
        bool singular = (m.singularity == Singularity::Singular) || (A.rowdim() != A.coldim());
        SolverReturnStatus status = SS_OK;
        if (!singular) {
            status = dixonSolve.solveNonsingular(xNum, xDen, A, b, false, maxTrials);
            singular = (status != SS_OK);
        }

        // Either A is known to be singular, or we just failed trying to solve it as non-singular.
        if (singular) {
            SolverLevel level = (m.certifyInconsistency ? SL_LASVEGAS : SL_MONTECARLO);

            if (m.singularSolutionType == SingularSolutionType::Diophantine) {
                DiophantineSolver<Solver> diophantineSolve(dixonSolve);
                status = diophantineSolve.diophantineSolve(xNum, xDen, A, b, maxTrials, level);
            }
            else {
                bool randomSolutionType = (m.singularSolutionType == SingularSolutionType::Random);
                status = dixonSolve.monolithicSolve(xNum, xDen, A, b, false, randomSolutionType, maxTrials, level);
            }
        }

        commentator().stop("solve.dixon.integer.dense");

        if (status == SS_INCONSISTENT) {
            throw LinboxMathInconsistentSystem("From Dixon method.");
        }
    }
}
