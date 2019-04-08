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
     * \brief Solve specialisation for Elimination.
     */
    template <class Matrix, class Vector, class CategoryTag>
    Vector& solve(Vector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Elimination& m)
    {
        // @fixme Make a DenseMatrix, compute sparsity,
        // and decide with magic whether to use Sparse or DenseElimination.
        // Ah, and do the very same within solveInPlace.
        return solve(x, A, b, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief Solve specialisation for Elimination with DenseMatrix.
     */
    template <class MatrixField, class Vector, class CategoryTag>
    Vector& solve(Vector& x, const DenseMatrix<MatrixField>& A, const Vector& b, const CategoryTag& tag,
                  const Method::Elimination& m)
    {
        return solve(x, A, b, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief Solve specialisation for Elimination with SparseMatrix.
     */
    template <class MatrixField, class Vector, class CategoryTag>
    Vector& solve(Vector& x, const SparseMatrix<MatrixField>& A, const Vector& b, const CategoryTag& tag,
                  const Method::Elimination& m)
    {
        return solve(x, A, b, tag, reinterpret_cast<const Method::SparseElimination&>(m));
    }

    //
    // solveInPlace
    //

    /**
     * \brief Solve in place specialisation for Elimination.
     */
    template <class Matrix, class Vector, class CategoryTag>
    Vector& solveInPlace(Vector& x, Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Elimination& m)
    {
        return solveInPlace(x, A, b, tag, reinterpret_cast<const Method::Blackbox&>(m));
    }

    /**
     * \brief Solve in place specialisation for Elimination with DenseMatrix.
     */
    template <class MatrixField, class Vector, class CategoryTag>
    Vector& solveInPlace(Vector& x, DenseMatrix<MatrixField>& A, const Vector& b, const CategoryTag& tag,
                         const Method::Elimination& m)
    {
        return solveInPlace(x, A, b, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief Solve in place specialisation for Elimination with SparseMatrix.
     */
    template <class MatrixField, class Vector, class CategoryTag>
    Vector& solveInPlace(Vector& x, SparseMatrix<MatrixField>& A, const Vector& b, const CategoryTag& tag,
                         const Method::Elimination& m)
    {
        return solveInPlace(x, A, b, tag, reinterpret_cast<const Method::SparseElimination&>(m));
    }
}
