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

#include <linbox/field/field-traits.h>
#include <linbox/solutions/methods-wip.h>

namespace LinBox {
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
        return solveInPlace(x, A, b, MethodWIP::Auto());
    }

    /**
     * \brief Rational solve in place Ax = b, for x expressed as xNum/xDen.
     * The matrix A might be modified.
     *
     * Second interface for solving in place, only valid for RingCategories::IntegerTag.
     */
    template <class Matrix, class Vector, class SolveMethod, class CategoryTag>
    inline void solveInPlace(Vector& xNum, typename Vector::Field::Element& xDen, Matrix& A, const Vector& b,
                             const CategoryTag& tag, const SolveMethod& m)
    {
        // @note This is called if the specialization has not been implemented,
        // which means there might not be any "in place" version of the solve.
        // So, we just forward to classic solve.

        solve(xNum, xDen, A, b, tag, m);
    }

    /**
     * \brief Rational solve in place dispatcher for automated category tag.
     */
    template <class Matrix, class Vector, class SolveMethod>
    inline void solveInPlace(Vector& xNum, typename Vector::Field::Element& xDen, Matrix& A, const Vector& b,
                             const SolveMethod& m)
    {
        solveInPlace(xNum, xDen, A, b, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief Rational solve in place dispatcher for automated solve method.
     */
    template <class Matrix, class Vector>
    inline void solveInPlace(Vector& xNum, typename Vector::Field::Element& xDen, Matrix& A, const Vector& b)
    {
        solveInPlace(xNum, xDen, A, b, MethodWIP::Auto());
    }
}

// @fixme Other solveInPlace
