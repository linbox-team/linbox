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

#include <linbox/algorithms/gauss.h>
#include <linbox/algorithms/matrix-hom.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/methods-wip.h>

namespace LinBox {
    /**
     * \brief Solve specialisation for SparseElimination.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag,
                        const MethodWIP::SparseElimination& m)
    {
        commentator().report(Commentator::LEVEL_UNIMPORTANT,
                             "Warning: Solve implicitly convert to a sparse matrix because of MethodWIP::SparseElimination.");

        using Field = typename Matrix::Field;
        SparseMatrix<Field> ASparse(A.field(), A.rowdim(), A.coldim());

        MatrixHom::map(ASparse, A);
        return solve(x, ASparse, b, tag, m);
    }

    /**
     * \brief Solve specialisation for SparseElimination with SparseMatrix.
     */
    template <class ResultVector, class... MatrixArgs, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const SparseMatrix<MatrixArgs...>& A, const Vector& b, const CategoryTag& tag,
                        const MethodWIP::SparseElimination& m)
    {
        commentator().start("solve.sparse-elimination.any.sparse");
        linbox_check((A.coldim() != x.size()) || (A.rowdim() != b.size()));

        // @fixme We should call solvein, that way the specialization above
        // would do the same and don't double copy the matrix.
        SparseMatrix<MatrixArgs...> ACopy(A);

        using Field = typename SparseMatrix<MatrixArgs...>::Field;
        GaussDomain<Field> gaussDomain(ACopy.field());
        gaussDomain.solveInPlace(x, ACopy, b);

        commentator().stop("solve.sparse-elimination.any.sparse");

        return x;
    }

    /**
     * \brief Solve specialisation for SparseElimination with SparseMatrix on IntegerTag.
     */
    template <class ResultVector, class... MatrixArgs, class Vector>
    ResultVector& solve(ResultVector& x, const SparseMatrix<MatrixArgs...>& A, const Vector& b,
                        const RingCategories::IntegerTag& tag, const MethodWIP::SparseElimination& m)
    {
        return solve(x, A, b, tag, reinterpret_cast<const MethodWIP::Dixon&>(m));
    }
}
