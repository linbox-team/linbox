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

#include <iostream> // @note Needed for givaro/ring-interface to compile

#include <linbox/algorithms/vector-fraction.h> // @fixme Why this is not inside vector/vector-fraction.h?
#include <linbox/field/field-traits.h>
#include <linbox/solutions/methods.h>

namespace LinBox {
    //
    // solve
    //

    /**
     * \brief Solve Ax = b, for x.
     *
     * Returns a vector x such that Ax = b.
     *
     * Specifically:
     * - A non singular: the unique solution is returned.
     * - A singular:
     *      - Consistent system: a random solution is returned.
     *          The method parameter can contain an hint that an
     *          arbitrary element of the solution space is acceptable instead,
     *          which can be faster to compute if one doesn't expect a result in that case.
     *      - Inconsistent system: LinboxMathInconsistentSystem is thrown.
     *      - Internal failure: LinboxError is thrown.
     *
     * CategoryTag is defaulted to `FieldTraits<Matrix::Field>::categoryTag()` when omitted.
     *
     * - Method::Auto
     *      - DenseMatrix   > Method::DenseElimination
     *      - SparseMatrix  > Method::SparseElimination
     *      - Otherwise
     *      |   - Row or column dimension < LINBOX_USE_BLACKBOX_THRESHOLD > Method::Elimination
     *      |   - Otherwise                                               > Method::Blackbox
     * - Method::Elimination
     *      - DenseMatrix   > Method::DenseElimination
     *      - SparseMatrix  > Method::SparseElimination
     *      - Otherwise     > Method::DenseElimination or Method::SparseElimination given matrix sparsity
     * - Method::DenseElimination
     *      - DenseMatrix
     *      |   - ModularTag > `LQUPMatrix<Field>::left_solve`
     *      |   - IntegerTag
     *      |   |   - RatVector > Method::Dixon
     *      |   |   - Otherwise > Error
     *      |   - Otherwise > Error (Not implemented)
     *      - Otherwise > Method::DenseElimination but copy to a DenseMatrix first
     * - Method::SparseElimination
     *      - SparseMatrix
     *      |   - IntegerTag > Method::Dixon
     *      |   - Otherwise > `GaussDomain<Field>::solveInPlace`
     *      - Otherwise > Method::SparseElimination but copy to SparseMatrix first
     * - Method::CRA
     *      - IntegerTag
     *      |   - Dispatch::Distributed > `ChineseRemainderDistributed`
     *      |   - Otherwise             > `RationalChineseRemainder`
     *      - Otherwise > Error
     * - Method::Dixon
     *      - IntegerTag
     *      |   - DenseMatrix   > `DixonSolver<..., Method::DenseElimination>`
     *      |   - SparseMatrix  > `DixonSolver<..., Method::SparseElimination>`
     *      |   - Otherwise     >  Error
     *      - Otherwise > Error
     * - Method::Blackbox > Method::Wiedemann
     * - Method::Wiedemann
     *      - ModularTag > `WiedemannSolver`
     *      - IntegerTag > Method::Dixon
     *      - Otherwise  > Error
     * - Method::BlockWiedemann [@deprecated, not tested]
     *      - ModularTag > `BlockWiedemannSolver`
     *      - Otherwise  > Error
     * - Method::Coppersmith [@deprecated, not tested]
     *      - ModularTag > `CoppersmithSolver`
     *      - Otherwise  > Error
     * - Method::Lanczos
     *      - ModularTag > `LanczosSolver`
     *      - Otherwise  > Error
     * - Method::BlockLanczos
     *      - ModularTag > `MGBlockLanczosSolver`
     *      - Otherwise  > Error
     * - Method::SymbolicNumericOverlap
     *      - IntegerTag
     *      |   - DenseMatrix > `RationalSolverSN<Ring, LPS<FMatrix>>`
     *      |   - Otherwise   > Error
     *      - Otherwise  > Error
     * - Method::SymbolicNumericNorm
     *      - IntegerTag > `DixonSolver<..., Method::SymbolicNumericNorm>`
     *      - Otherwise  > Error
     *
     * @param [out] x solution, can be a rational solution (vector of numerators and one denominator)
     * @param [in]  A matrix
     * @param [in]  b target
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return reference to \p x
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag, class SolveMethod>
    inline ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const SolveMethod& m)
    {
        throw LinBoxError("solve<" + CategoryTag::name() + ", " + SolveMethod::name() + "> does not exists.");
    }

    /**
     * \brief Solve dispatcher for automated category tag.
     */
    template <class ResultVector, class Matrix, class Vector, class SolveMethod>
    inline ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const SolveMethod& m)
    {
        return solve(x, A, b, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief Solve dispatcher for automated solve method.
     */
    template <class ResultVector, class Matrix, class Vector>
    inline ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b)
    {
        return solve(x, A, b, Method::Auto());
    }

    /**
     * \brief Solve specialisation on RationalTag, with a generic method.
     */
    template <class RatVector, class RatMatrix, class Vector, class SolveMethod>
    RatVector& solve(RatVector& x, const RatMatrix& A, const Vector& b, const RingCategories::RationalTag& tag, const SolveMethod& m)
    {
        return solve(x, A, b, tag, Method::CRA<SolveMethod>(m));
    }

    /**
     * \brief Solve specialisation on IntegerTag with Vector<QField> as result.
     *
     * This forward to the rational interface (num, den).
     * But will only work if the ResultVector if a vector of some Rational type.
     */
    template <class RatVector, class Matrix, class Vector, class SolveMethod>
    typename std::enable_if<std::is_same<typename SolveMethod::CategoryTag, RingCategories::IntegerTag>::value
                            && std::is_same<typename FieldTraits<typename RatVector::Field>::categoryTag, RingCategories::RationalTag>::value,
                            RatVector&>::type
    solve(RatVector& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag, const SolveMethod& m)
    {
        using Ring = typename Matrix::Field;
        using Element = typename Ring::Element;

        Vector xNum(b.field(), x.size());
        Element xDen;

        solve(xNum, xDen, A, b, tag, m);

        // The denominator being zero means computation failure
        if (b.field().isZero(xDen)) {
            throw LinboxError("Rational solve failed.");
        }

        // Copy result back to RatVector
        auto iXNum = xNum.begin();
        for (auto iX = x.begin(); iX != x.end(); ++iX) {
            *iX = typename RatVector::value_type(*iXNum, xDen);
            ++iXNum;
        }

        return x;
    }

    /**
     * \brief Solve specialisation on IntegerTag with VectorFraction as result.
     *
     * This forward to the rational interface (num, den).
     */
    template <class Matrix, class Vector, class SolveMethod>
    typename std::enable_if<std::is_same<typename SolveMethod::CategoryTag, RingCategories::IntegerTag>::value,
                            VectorFraction<typename Matrix::Field>&>::type
    solve(VectorFraction<typename Matrix::Field>& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag,
          const SolveMethod& m)
    {
        solve(x.numer, x.denom, A, b, tag, m);
        return x;
    }

    /**
     * \brief Rational solve Ax = b, for x expressed as xNum/xDen.
     *
     * Second interface for solving, only valid for RingCategories::IntegerTag.
     *
     * Solve with this interface will usually go for CRA or Dixon lifting,
     * as non-modular elimination would snowball elements to very big values.
     */
    template <class IntVector, class Matrix, class Vector, class CategoryTag, class SolveMethod>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const Matrix& A, const Vector& b, const CategoryTag& tag,
                      const SolveMethod& m)
    {
        throw LinBoxError("Rational solve is only valid for RingCategories::IntegerTag.");
    }

    /**
     * \brief Rational solve dispatcher for unimplemented methods.
     */
    template <class IntVector, class Matrix, class Vector, class SolveMethod>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const SolveMethod& m)
    {
        commentator().report(Commentator::LEVEL_UNIMPORTANT, ("Warning: Rational solve on RingCategories::IntegerTag with "
                                                              + SolveMethod::name() + " is forwarded to Method::Dixon instead.")
                                                                 .c_str());

        solve(xNum, xDen, A, b, tag, reinterpret_cast<const Method::Dixon&>(m));
    }

    /**
     * \brief Rational solve dispatcher for automated category tag.
     */
    template <class IntVector, class Matrix, class Vector, class SolveMethod>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const Matrix& A, const Vector& b, const SolveMethod& m)
    {
        solve(xNum, xDen, A, b, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief Rational solve dispatcher for automated solve method.
     */
    template <class IntVector, class Matrix, class Vector>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const Matrix& A, const Vector& b)
    {
        solve(xNum, xDen, A, b, Method::Auto());
    }

    //
    // solveInPlace
    //

    /**
     * \brief Solve Ax = b, for x.
     *
     * Returns a vector x such that Ax = b.
     * A can be modified.
     *
     * See documentation for `solve`.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag, class SolveMethod>
    inline ResultVector& solveInPlace(ResultVector& x, Matrix& A, const Vector& b, const CategoryTag& tag, const SolveMethod& m)
    {
        // @note This is called if the specialization has not been implemented,
        // which means there might not be any "in place" version of the solve.
        // So, we just forward to classic solve.

        return solve(x, A, b, tag, m);
    }

    /**
     * \brief Solve in place dispatcher for automated category tag.
     */
    template <class ResultVector, class Matrix, class Vector, class SolveMethod>
    inline ResultVector& solveInPlace(ResultVector& x, Matrix& A, const Vector& b, const SolveMethod& m)
    {
        return solveInPlace(x, A, b, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief Solve in place dispatcher for automated solve method.
     */
    template <class ResultVector, class Matrix, class Vector>
    inline ResultVector& solveInPlace(ResultVector& x, Matrix& A, const Vector& b)
    {
        return solveInPlace(x, A, b, Method::Auto());
    }

    /**
     * \brief Rational solve in place Ax = b, for x expressed as xNum/xDen.
     * The matrix A might be modified.
     *
     * Second interface for solving in place, only valid for RingCategories::IntegerTag.
     */
    template <class IntVector, class Matrix, class Vector, class SolveMethod, class CategoryTag>
    inline void solveInPlace(IntVector& xNum, typename IntVector::Element& xDen, Matrix& A, const Vector& b, const CategoryTag& tag,
                             const SolveMethod& m)
    {
        // @note This is called if the specialization has not been implemented,
        // which means there might not be any "in place" version of the solve.
        // So, we just forward to classic solve.

        solve(xNum, xDen, A, b, tag, m);
    }

    /**
     * \brief Rational solve in place dispatcher for automated category tag.
     */
    template <class IntVector, class Matrix, class Vector, class SolveMethod>
    inline void solveInPlace(IntVector& xNum, typename IntVector::Element& xDen, Matrix& A, const Vector& b, const SolveMethod& m)
    {
        solveInPlace(xNum, xDen, A, b, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief Rational solve in place dispatcher for automated solve method.
     */
    template <class IntVector, class Matrix, class Vector>
    inline void solveInPlace(IntVector& xNum, typename IntVector::Element& xDen, Matrix& A, const Vector& b)
    {
        solveInPlace(xNum, xDen, A, b, Method::Auto());
    }
}

//
// solve
//

#include "./solve/solve-auto.h"

// Elimination
#include "./solve/solve-dense-elimination.h"
#include "./solve/solve-elimination.h"
#include "./solve/solve-sparse-elimination.h"

// Integer-based
#include "./solve/solve-cra.h"
#include "./solve/solve-dixon.h"
#include "./solve/solve-numeric-symbolic.h"

// Blackbox
#include "./solve/solve-blackbox.h"
#include "./solve/solve-lanczos.h"
#include "./solve/solve-wiedemann.h"
