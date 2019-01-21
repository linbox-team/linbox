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

#include <iostream> // @fixme This is needed for givaro ring-interface to compile

#include <linbox/field/field-traits.h>
#include <linbox/solutions/methods.h>
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
     *      - Inconsistent system: zero vector is returned.
     *
     * CategoryTag is defaulted to FieldTraits<Matrix::Field>::categoryTag().
     *
     * SolveMethod determines which algorithm will be used.
     * Using Method::Auto chooses the Method according to Matrix and Ring information:
     * - RingCategories::IntegerTag         => Method::Dixon @fixme Requires correct interface
     * - RingCategories::<Other>            => Method::Elimination
     * - Method::Elimination
     *      - Matrix == DenseMatrix         => Method::DenseElimination
     *      - Matrix == SparseMatrix        => Method::SparseElimination
     *      - Matrix == <Other>             => Method::Blackbox
     * - Method::BlackBox
     *      - Matrix == <Other>             => Method::Wiedemann @fixme
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
        throw NotImplementedYet("Solve specialisation is not implemented yet.");
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
        throw LinBoxError("Rational solve is only valid for RingCategories::IntegerTag and Method::Cra or Method::Dixon.");
    }

    template <class Matrix, class Vector, class CategoryTag, class SolveMethod>
    inline void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const SolveMethod& m)
    {
        throw NotImplementedYet("Rational solve specialisation is not implemented yet.");
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
        return solve(xNum, xDen, A, b, Method::Auto());
    }
}

#include "./solve/solve-utils.h"

#include "./solve/solve-auto.h"
#include "./solve/solve-cra.h"
#include "./solve/solve-dixon.h"

#include "./solve/solve-dense-elimination.h"
#include "./solve/solve-elimination.h"
#include "./solve/solve-sparse-elimination.h"

#include "./solve/solvein.h" // @fixme

#endif
