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
#include <linbox/solutions/methods.h>

namespace LinBox {
    //
    // solve
    //

    /**
     * \brief Solve specialisation for SparseElimination.
     */
    template <class Matrix, class Vector>
    Vector& solve(Vector& x, const Matrix& A, const Vector& b, const RingCategories::ModularTag& tag, const Method::SparseElimination& m)
    {
        commentator().report(Commentator::LEVEL_UNIMPORTANT,
                             "Warning: Solve implicitly convert to a sparse matrix because of Method::SparseElimination.");

        using Field = typename Matrix::Field;
        SparseMatrix<Field> ASparse(A.field(), A.rowdim(), A.coldim());
        MatrixHom::map(ASparse, A);
        return solveInPlace(x, ASparse, b, tag, m);
    }

    /**
     * \brief Solve specialisation for SparseElimination with SparseMatrix.
     */
    template <class... MatrixArgs, class Vector>
    Vector& solve(Vector& x, const SparseMatrix<MatrixArgs...>& A, const Vector& b, const RingCategories::ModularTag& tag,
                  const Method::SparseElimination& m)
    {
        SparseMatrix<MatrixArgs...> ACopy(A);
        return solveInPlace(x, ACopy, b, tag, m);
    }

    //
    // solveInPlace
    //

    /**
     * \brief Solve in place specialisation for SparseElimination with SparseMatrix.
     */
    template <class... MatrixArgs, class Vector>
    Vector& solveInPlace(Vector& x, SparseMatrix<MatrixArgs...>& A, const Vector& b, const RingCategories::ModularTag& tag,
                         const Method::SparseElimination& m)
    {
        commentator().start("solve-in-place.sparse-elimination.any.sparse");
        linbox_check((A.coldim() == x.size()) && (A.rowdim() == b.size()));

        using Field = typename SparseMatrix<MatrixArgs...>::Field;
        GaussDomain<Field> gaussDomain(A.field());
        gaussDomain.solveInPlace(x, A, b);

        commentator().stop("solve-in-place.sparse-elimination.any.sparse");

        return x;
    }
}
