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

#include <linbox/matrix/densematrix/blas-matrix.h>
#include <linbox/matrix/sparse-matrix.h>

namespace LinBox {
    /**
     * \brief Solve specialisation for BlasElimination.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag,
                        const Method::BlasElimination& m)
    {
        // @fixme Original code would copy the sparse or other to a BlasMatrix
        throw NotImplementedYet("Blas eliminating.");
    }
    /**
     * \brief Solve specialisation for BlasElimination on dense matrices with ModularTag.
     */
    template <class Field, class Vector>
    Vector& solve(Vector& x, const BlasMatrix<Field>& A, const Vector& b, const RingCategories::ModularTag& tag,
                        const Method::BlasElimination& m)
    {
        solve_precheck(x, A, b);

        commentator().start("Solve Modular BlasElimination for BlasMatrix", "solve.modular.blas-elimination.dense");

        LQUPMatrix<Field> LQUP(A);
        LQUP.left_solve(x, b);

        commentator().stop("solve.modular.blas-elimination.dense");

        return x;
    }

    /**
     * \brief Solve specialisation for BlasElimination on dense matrices with IntegerTag.
     */
    template <class ResultVector, class Field, class Vector>
    ResultVector& solve(ResultVector& x, const BlasMatrix<Field>& A, const Vector& b, const RingCategories::IntegerTag& tag,
                        const Method::BlasElimination& m)
    {
        // @fixme This is the original code for this case... but it goes for a Dixon!
        // Is that really what we want?

        // @fixme By the way, what is Method::NonBlasElimination?
        return solve(x, A, b, tag, Method::Dixon(m));
    }
}
