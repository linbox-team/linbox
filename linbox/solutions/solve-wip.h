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

#ifndef __LINBOX_solve_H
#define __LINBOX_solve_H

#include <iostream> // @note Needed for givaro/ring-interface to compile

#include <linbox/field/field-traits.h>
#include <linbox/solutions/methods-wip.h>
#include <linbox/util/debug.h> // NotImplementedYet

namespace LinBox {

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
     *      - Inconsistent system: zero vector is returned. @fixme Not always true...
     *
     * CategoryTag is defaulted to FieldTraits<Matrix::Field>::categoryTag().
     *
     * - Method::Auto
     *      - DenseMatrix   > Method::Elimination
     *      - Otherwise     > Method::Blackbox (or Elimination given the size of the matrix)
     *      - @fixme Why don't use Elimination if sparseMatrix?
     * - Method::Elimination
     *      - SparseMatrix  > Method::SparseElimination
     *      - Otherwise     > Method::DenseElimination
     * - Method::DenseElimination
     *      - DenseMatrix
     *      |   - ModularTag > `LQUPMatrix<Field>::left_solve`
     *      |   - IntegerTag
     *      |   |   - RatVector > Method::Dixon
     *      |   |   - Otherwise > Error
     *      |   - Otherwise > @fixme NIY
     *      - Otherwise > Method::DenseElimination but copy to a DenseMatrix first
     * - Method::SparseElimination
     *      - SparseMatrix
     *          - IntegerTag > Method::Dixon
     *          - Otherwise > `GaussDomain<GF2>::solveInPlace`
     *      - Otherwise > Method::SparseElimination but copy to SparseMatrix first
     * - Method::Cra
     *      - IntegerTag
     *      |   - Dispatch::Distributed > `MPIratChineseRemainder`
     *      |   - Otherwise             > `RationalRemainder`
     *      - Otherwise > Error
     * - Method::Dixon
     *      - IntegerTag
     *      |   - DenseMatrix   > `RationalSolver<..., Method::Dixon>`
     *      |   - SparseMatrix  > `RationalSolver<..., Method::SparseElimination>`
     *      |   - Otherwise     > @fixme NIY Does dixon need to read the matrix?
     *      - Otherwise > Error
     * - Method::Blackbox > Method::Wiedemann
     * - Method::Wiedemann
     *      - ModularTag > `WiedemannSolver`
     *      - Otherwise > Error
     * - Method::BlockWiedemann [@deprecated, not tested]
     *      - ModularTag > `BlockWiedemannSolver`
     *      - Otherwise > Error
     * - Method::Coppersmith [@deprecated, not tested]
     *      - ModularTag > `CoppersmithSolver`
     *      - Otherwise > Error
     * - @fixme Lanczos and others... (from bbsolve.h)
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
        return solve(x, A, b, MethodWIP::Auto());
    }

    /**
     * \brief Solve specialisation on IntegerTag.
     *
     * This forward to the rational interface (num, den).
     * But will only work if the ResultVector if a vector of some Rational type.
     */
    template <class ResultVector, class Matrix, class Vector, class SolveMethod>
    typename std::enable_if<std::is_same<typename SolveMethod::CategoryTag, RingCategories::IntegerTag>::value,
                            ResultVector&>::type
    solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag, const SolveMethod& m)
    {
        using Ring = typename Vector::Field;
        using Element = typename Ring::Element;

        BlasVector<Ring> xNum(x.field(), x.size());
        Element xDen;

        solve(xNum, xDen, A, b, tag, m);

        // Copy result back to ResultVector
        auto iXNum = xNum.begin();
        for (auto iX = x.begin(); iX != x.end(); ++iX) {
            *iX = typename ResultVector::value_type(*iXNum, xDen);
            ++iXNum;
        }

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
    template <class Matrix, class Vector, class CategoryTag, class SolveMethod>
    inline void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b,
                      const CategoryTag& tag, const SolveMethod& m)
    {
        throw LinBoxError("Rational solve is only valid for RingCategories::IntegerTag.");
    }

    template <class Matrix, class Vector, class SolveMethod>
    inline void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const SolveMethod& m)
    {
        throw LinBoxError("Rational solve<RingCategories::IntegerTag, " + SolveMethod::name() + "> does not exists.");
    }

    /**
     * \brief Solve dispatcher for automated category tag.
     */
    template <class Matrix, class Vector, class SolveMethod>
    inline void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b, const SolveMethod& m)
    {
        return solve(xNum, xDen, A, b, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief Solve dispatcher for automated solve method.
     */
    template <class Matrix, class Vector>
    inline void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b)
    {
        return solve(xNum, xDen, A, b, MethodWIP::Auto());
    }
}

#include "./solve/solve-auto.h"
// #include "./solve/solve-blackbox.h"

// Elimination
#include "./solve/solve-dense-elimination.h"
#include "./solve/solve-elimination.h"
#include "./solve/solve-sparse-elimination.h"

// Integer-based
#include "./solve/solve-cra.h"
#include "./solve/solve-dixon.h"

// @fixme What are those for?
#include "./solve/solve-wiedemann.h"

// #include "./solve/solvein.h" @fixme

#endif
