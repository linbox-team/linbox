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

#include <linbox/algorithms/rational-solver.h>
#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/matrix/sparse-matrix.h>

namespace LinBox {
    /**
     * \brief Solve specialisation for Dixon.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Dixon& m)
    {
        throw NotImplementedYet("Dixon solving.");
    }

    /**
     * \brief Solve specialisation for Dixon.
     */
    template <class ResultVector, class Matrix, class Vector>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::ModularTag& tag,
                        const Method::Dixon& m)
    {
        throw LinBoxFailure("Cannot solve with Dixon method and Modular category tag.");
    }

    /**
     * \brief Solve specialisation for Dixon on dense matrices.
     */
    // @fixme Method::Dixon should be templated with IterationMethod too!
    template <class ResultVector, class MatrixField, class Vector>
    ResultVector& solve(ResultVector& x, const BlasMatrix<MatrixField>& A, const Vector& b, const RingCategories::IntegerTag& tag,
                        const Method::Dixon& m)
    {
        solve_precheck(x, A, b);

        // @fixme This is the original code solve.h:534 for this case... and it's ugly!
        commentator().start("Solve Integer Dixon for BlasMatrix", "solve.integer.dixon.dense");

        using Field = Givaro::Modular<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()));

        // @fixme This RationalSolver should be called something with Dixon in it...
        RationalSolver<MatrixField, Field, PrimeGenerator, DixonTraits> dixonSolve(A.field(), primeGenerator);

        // Do we know anything about A singularity?
        // @fixme implement

        // SolverReturnStatus status = SS_OK;
        // if singularity unknown and matrix is square, we try nonsingular solver
        // switch (m.singular()) {
        // case Specifier::SINGULARITY_UNKNOWN:
        //     switch (A.rowdim() == A.coldim() ? status = dixonSolve.solveNonsingular(x, d, A, b, false, (int)m.maxTries())
        //                                      : SS_SINGULAR) {

        //     case SS_OK: m.singular(Specifier::NONSINGULAR); break;
        //     case SS_SINGULAR:
        //         switch (m.solution()) {
        //         case DixonTraits::DETERMINIST:
        //             status = dixonSolve.monolithicSolve(x, d, A, b, false, false, (int)m.maxTries(),
        //                                                 (m.certificate() ? SL_LASVEGAS : SL_MONTECARLO));
        //             break;
        //         case DixonTraits::RANDOM:
        //             status = dixonSolve.monolithicSolve(x, d, A, b, false, true, (int)m.maxTries(),
        //                                                 (m.certificate() ? SL_LASVEGAS : SL_MONTECARLO));
        //             break;
        //         case DixonTraits::DIOPHANTINE: {
        //             DiophantineSolver<RationalSolver<Ring, Field, PrimeIterator<IteratorCategories::HeuristicTag>,
        //             DixonTraits>>
        //                 dsolve(dixonSolve);
        //             status =
        //                 dsolve.diophantineSolve(x, d, A, b, (int)m.maxTries(), (m.certificate() ? SL_LASVEGAS :
        //                 SL_MONTECARLO));
        //         } break;
        //         default: break;
        //         }
        //         break;
        //     default: break;

        //     }
        //     break;

        // case Specifier::NONSINGULAR: dixonSolve.solveNonsingular(x, d, A, b, false, (int)m.maxTries()); break;

        // case Specifier::SINGULAR:
        //     switch (m.solution()) {
        //     case DixonTraits::DETERMINIST:
        //         status = dixonSolve.monolithicSolve(x, d, A, b, false, false, (int)m.maxTries(),
        //                                             (m.certificate() ? SL_LASVEGAS : SL_MONTECARLO));
        //         break;

        //     case DixonTraits::RANDOM:
        //         status = dixonSolve.monolithicSolve(x, d, A, b, false, true, (int)m.maxTries(),
        //                                             (m.certificate() ? SL_LASVEGAS : SL_MONTECARLO));
        //         break;

        //     case DixonTraits::DIOPHANTINE: {
        //         DiophantineSolver<RationalSolver<Ring, Field, PrimeIterator<IteratorCategories::HeuristicTag>, DixonTraits>>
        //             dsolve(rsolve);
        //         status = dsolve.diophantineSolve(x, d, A, b, (int)m.maxTries(), (m.certificate() ? SL_LASVEGAS :
        //         SL_MONTECARLO));
        //     } break;

        //         // default:
        //         //  break;
        //     }
        // default: break;
        // }

        commentator().stop("solve.integer.dixon.dense");

        // if (status == SS_INCONSISTENT) {
        //     // @fixme The solve specification says we return a zero vector in that case...
        //     // But we throw instead, update spec.
        //     throw LinboxMathInconsistentSystem("Solve with Dixon method is impossible.");
        // }

        return x;
    }
}
