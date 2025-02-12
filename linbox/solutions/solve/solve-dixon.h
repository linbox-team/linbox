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
#include <linbox/algorithms/rational-solver.h>
#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/matrix/sparse-matrix.h>

namespace LinBox {
    namespace {
        template <class Matrix>
        struct MethodForMatrix {
            using type = Method::Wiedemann; // For blackboxes
        };

        template <class Ring>
        struct MethodForMatrix<DenseMatrix<Ring>> {
            using type = Method::DenseElimination;
        };

        template <class Ring>
        struct MethodForMatrix<SparseMatrix<Ring>> {
            using type = Method::SparseElimination;
        };
    }

    /**
     * \brief Solve specialisation for Dixon on blackboxes matrices.
     */
    template <class IntVector, class Blackbox, class Vector>
    void solve(IntVector& xNum, typename IntVector::Element& xDen, const Blackbox& A, const Vector& b,
               const RingCategories::IntegerTag& tag, const Method::Dixon& m)
    {
        commentator().start("solve.dixon.integer.blackbox");
        linbox_check((A.coldim() == xNum.size()) && (A.rowdim() == b.size()));

        using Ring = typename Blackbox::Field;
        using Field = Givaro::Modular<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()));

        using Solver = DixonSolver<Ring, Field, PrimeGenerator, typename MethodForMatrix<Blackbox>::type>;
        Solver dixonSolve(A.field(), primeGenerator);

        // @fixme I'm still bit sad that we cannot use generically the function below,
        // just because RationalSolve<..., SparseElimination> has not the same
        // API (i.e. no solveNonSingular) than DixonSolver<..., DenseElimination>
        int maxTrials = m.trialsBeforeFailure;
        SolverReturnStatus status = dixonSolve.solve(xNum, xDen, A, b, maxTrials);

        commentator().stop("solve.dixon.integer.blackbox");

        if (status == SS_INCONSISTENT) {
            throw LinboxMathInconsistentSystem("From Dixon method.");
        }
        else if (status == SS_FAILED || status == SS_BAD_PRECONDITIONER) {
            throw LinboxError("From Dixon method.");
        }
    }

    /**
     * \brief Solve specialisation for Dixon on dense matrices.
     */
    template <class IntVector, class Ring, class Vector>
    void solve(IntVector& xNum, typename IntVector::Element& xDen, const DenseMatrix<Ring>& A, const Vector& b,
               const RingCategories::IntegerTag& tag, const Method::Dixon& m)
    {
        commentator().start("solve.dixon.integer.dense");
        linbox_check((A.coldim() == xNum.size()) && (A.rowdim() == b.size()));

        // @fixme Using Givaro::ModularBalanced<double> for the field makes Dixon fail...
        using Matrix = DenseMatrix<Ring>;
        using Field = Givaro::Modular<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()));

        using Solver = DixonSolver<Ring, Field, PrimeGenerator, typename MethodForMatrix<Matrix>::type>;
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
            if (m.singularSolutionType == SingularSolutionType::Diophantine) {
                SolverLevel level = (m.certifyInconsistency ? SL_LASVEGAS : SL_MONTECARLO);
                DiophantineSolver<Solver> diophantineSolve(dixonSolve);
                status = diophantineSolve.diophantineSolve(xNum, xDen, A, b, maxTrials, level);
            }
            else {
                status = dixonSolve.monolithicSolve(xNum, xDen, A, b, m);
            }
        }

        commentator().stop("solve.dixon.integer.dense");

        if (status == SS_INCONSISTENT) {
            throw LinboxMathInconsistentSystem("From Dixon method.");
        }
        else if (status == SS_FAILED || status == SS_BAD_PRECONDITIONER) {
            throw LinboxError("From Dixon method.");
        }
    }

    /**
     * \brief Solve specialisation for Dixon on sparse matrices.
     */
    template <class IntVector, class... MatrixArgs, class Vector>
    void solve(IntVector& xNum, typename IntVector::Element& xDen, const SparseMatrix<MatrixArgs...>& A, const Vector& b,
               const RingCategories::IntegerTag& tag, const Method::Dixon& m)
    {
        commentator().start("solve.dixon.integer.sparse");
        linbox_check((A.coldim() == xNum.size()) && (A.rowdim() == b.size()));

        using Matrix = SparseMatrix<MatrixArgs...>;
        using Ring = typename Matrix::Field;
        using Field = Givaro::Modular<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()));

        using Solver = DixonSolver<Ring, Field, PrimeGenerator, typename MethodForMatrix<Matrix>::type>;
        Solver dixonSolve(A.field(), primeGenerator);

        // @fixme I'm a bit sad that we cannot use generically the function above,
        // just because RationalSolve<..., SparseElimination> has not the same
        // API (i.e. no solveNonSingular) than DixonSolver<..., DenseElimination>
        int maxTrials = m.trialsBeforeFailure;
        SolverReturnStatus status = dixonSolve.solve(xNum, xDen, A, b, maxTrials);

        commentator().stop("solve.dixon.integer.sparse");

        if (status == SS_INCONSISTENT) {
            throw LinboxMathInconsistentSystem("From Dixon method.");
        }
        else if (status == SS_FAILED || status == SS_BAD_PRECONDITIONER) {
            throw LinboxError("From Dixon method.");
        }
    }
}
