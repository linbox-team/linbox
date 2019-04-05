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

#include <linbox/matrix/dense-matrix.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/methods.h>

namespace LinBox {
    //
    // solve
    //

    /**
     * \brief Solve specialisation for Auto.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Auto& m)
    {
        if (useBlackboxMethod(A)) {
            return solve(x, A, b, tag, reinterpret_cast<const Method::Blackbox&>(m));
        }
        else {
            return solve(x, A, b, tag, reinterpret_cast<const Method::Elimination&>(m));
        }
    }

    /**
     * \brief Solve specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class ResultVector, class Field, class Vector>
    ResultVector& solve(ResultVector& x, const DenseMatrix<Field>& A, const Vector& b,
                        const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return solve(x, A, b, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief Solve specialisation for Auto with SparseMatrix and ModularTag.
     */
    template <class ResultVector, class... MatrixArgs, class Vector>
    ResultVector& solve(ResultVector& x, const SparseMatrix<MatrixArgs...>& A, const Vector& b,
                        const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return solve(x, A, b, tag, reinterpret_cast<const Method::SparseElimination&>(m));
    }

    /**
     * \brief Solve specialisation for Auto and IntegerTag.
     */
    template <class ResultVector, class Matrix, class Vector>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag,
                        const Method::Auto& m)
    {
        return solve(x, A, b, tag, reinterpret_cast<const Method::Dixon&>(m));
    }

    /**
     * \brief Solve specialisation for Auto and RationalTag.
     */
    template <class ResultVector, class Matrix, class Vector>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::RationalTag& tag,
                        const Method::Auto& m)
    {
        return solve(x, A, b, tag, Method::CRAAuto(m));
    }

    //
    // solve Rational API
    //

    /**
     * \brief Solve specialization for Auto and IntegerTag.
     */
    template <class IntVector, class Matrix, class Vector>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const Method::Auto& m)
    {
        solve(xNum, xDen, A, b, tag, reinterpret_cast<const Method::Dixon&>(m));
    }

    //
    // solveInPlace
    //

    /**
     * \brief Solve in place specialisation for Auto.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solveInPlace(ResultVector& x, Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Auto& m)
    {
        if (useBlackboxMethod(A)) {
            return solve(x, A, b, tag, reinterpret_cast<const Method::Blackbox&>(m));
        }
        else {
            return solve(x, A, b, tag, reinterpret_cast<const Method::Elimination&>(m));
        }
    }

    /**
     * \brief Solve in place specialisation for Auto and IntegerTag.
     */
    template <class ResultVector, class Matrix, class Vector>
    ResultVector& solveInPlace(ResultVector& x, Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag,
                               const Method::Auto& m)
    {
        return solveInPlace(x, A, b, tag, reinterpret_cast<const Method::Dixon&>(m));
    }

    /**
     * \brief Solve in place specialisation for Auto with DenseMatrix and non-IntegerTag.
     */
    template <class ResultVector, class Field, class Vector, class CategoryTag>
    typename std::enable_if<!std::is_same<CategoryTag, RingCategories::IntegerTag>::value, ResultVector&>::type solveInPlace(
        ResultVector& x, DenseMatrix<Field>& A, const Vector& b, const CategoryTag& tag, const Method::Auto& m)
    {
        return solveInPlace(x, A, b, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief Solve in place specialisation for Auto with SparseMatrix and non-IntegerTag.
     */
    template <class ResultVector, class... MatrixArgs, class Vector, class CategoryTag>
    typename std::enable_if<!std::is_same<CategoryTag, RingCategories::IntegerTag>::value, ResultVector&>::type solveInPlace(
        ResultVector& x, SparseMatrix<MatrixArgs...>& A, const Vector& b, const CategoryTag& tag, const Method::Auto& m)
    {
        return solveInPlace(x, A, b, tag, reinterpret_cast<const Method::SparseElimination&>(m));
    }

    //
    // solveInPlace Rational API
    //

    /**
     * \brief Solve in place specialization for Auto and IntegerTag.
     */
    template <class Matrix, class Vector>
    inline void solveInPlace(Vector& xNum, typename Vector::Element& xDen, Matrix& A, const Vector& b,
                             const RingCategories::IntegerTag& tag, const Method::Auto& m)
    {
        solveInPlace(xNum, xDen, A, b, tag, reinterpret_cast<const Method::Dixon&>(m));
    }
}
