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

#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/methods.h>

namespace LinBox {
    /**
     * \brief Solve specialisation for SparseElimination.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag,
                        const Method::SparseElimination& m)
    {
        // @fixme Emit warning about conversion to sparse

        using Field = typename Matrix::Field;
        SparseMatrix<Field> ASparse(A.field(), A.rowdim(), A.coldim());

        // @fixme Isn't that old syntax?
        MatrixHom::map(ASparse, A);

        return solve(x, ASparse, b, tag, m);
    }

    /**
     * \brief Solve specialisation for SparseElimination with SparseMatrix.
     */
    template <class ResultVector, class... MatrixArgs, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const SparseMatrix<MatrixArgs...>& A, const Vector& b, const CategoryTag& tag,
                        const Method::SparseElimination& m)
    {
        solve_precheck(x, A, b);

        commentator().start("Solve SparseElimination for SparseMatrix", "solve.any.sparse-elimination.sparse");

        // @fixme...
        // using Field = typename SparseMatrix<MatrixArgs...>::Field;
        // GaussDomain<Field> gaussDomain(A.field());
        // gaussDomain.solvein(x, A, b, generator);

        commentator().stop("solve.any.sparse-elimination.sparse");

        return x;
    }
}
