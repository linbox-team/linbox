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
#include <linbox/algorithms/rational-solver.h> // @fixme Rename dixon-rational-solver?
#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/matrix/sparse-matrix.h>

namespace LinBox {
    // @fixme TBR when RationalSolver is renamed
    template <class... Args>
    using DixonRationalSolver = RationalSolver<Args...>;

    /**
     * \brief Solve specialisation for Dixon.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const MethodWIP::Dixon& m)
    {
        throw NotImplementedYet("Dixon solving.");
    }

    /**
     * \brief Solve specialisation for Dixon.
     */
    template <class ResultVector, class Matrix, class Vector, class IterationMethod>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag,
                        const MethodWIP::DixonCustom<IterationMethod>& m)
    {
        throw LinBoxFailure("Solve with MethodWIP::Dixon expects the rational result interface to be used.");
    }

    /**
     * \brief Solve specialisation for Dixon on dense matrices.
     */
    template <class Matrix, class Vector, class CategoryTag, class IterationMethod>
    void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b, const CategoryTag& tag,
               const MethodWIP::DixonCustom<IterationMethod>& m)
    {
        throw LinBoxFailure("Solve with MethodWIP::Dixon expects RingCategories::IntegerTag.");
    }

    /**
     * \brief Solve specialisation for Dixon on dense matrices.
     */
    // @fixme MethodWIP::Dixon should be templated with IterationMethod too!
    template <class MatrixField, class Vector, class IterationMethod>
    void solve(Vector& xNum, typename Vector::Field::Element& xDen, const BlasMatrix<MatrixField>& A, const Vector& b,
               const RingCategories::IntegerTag& tag, const MethodWIP::DixonCustom<IterationMethod>& m)
    {
        commentator().start("solve.dixon.integer.dense");
        linbox_check((A.coldim() != xNum.size()) || (A.rowdim() != b.size()));

        using Field = Givaro::Modular<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()));

        // @fixme Why can't I use Method::Auto here?
        using Solver = DixonRationalSolver<MatrixField, Field, PrimeGenerator, Method::Dixon>;
        Solver dixonSolve(A.field(), primeGenerator);

        // Either A is known to be non-singular, or we just don't know yet.
        int maxTrials = m.trialsBeforeThrowing;
        bool singular = (m.singularity == Singularity::Singular) || (A.rowdim() != A.coldim());
        if (!singular) {
            auto status = dixonSolve.solveNonsingular(xNum, xDen, A, b, false, maxTrials);
            singular = (status != SS_OK);
        }

        // Either A is known to be singular, or we just failed trying to solve it as non-singular.
        if (singular) {
            SolverLevel level = (m.findInconsistencyCertificate ? SL_LASVEGAS : SL_MONTECARLO);

            if (m.solutionType == SolutionType::Diophantine) {
                DiophantineSolver<Solver> diophantineSolve(dixonSolve);
                diophantineSolve.diophantineSolve(xNum, xDen, A, b, maxTrials, level);
            }
            else {
                bool randomSolutionType = (m.solutionType == SolutionType::Random);
                dixonSolve.monolithicSolve(xNum, xDen, A, b, false, randomSolutionType, maxTrials, level);
            }
        }

        commentator().stop("solve.dixon.integer.dense");

        // @fixme Tests should check conformance to that
        // if (status == SS_INCONSISTENT) {
        //     // @fixme The solve specification says we return a zero vector in that case...
        //     // But we throw instead, update spec.
        //     throw LinboxMathInconsistentSystem("Solve with Dixon method is impossible.");
        // }
    }
}
