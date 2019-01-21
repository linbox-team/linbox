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
    /**
     * \brief Solve specialisation for DenseElimination.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag,
                        const Method::DenseElimination& m)
    {
        // @fixme Original code would copy the sparse or other to a DenseMatrix
        throw NotImplementedYet("Dense eliminating.");
    }

    /**
     * \brief Solve specialisation for DenseElimination on dense matrices with ModularTag.
     */
    template <class Field, class Vector>
    Vector& solve(Vector& x, const DenseMatrix<Field>& A, const Vector& b, const RingCategories::ModularTag& tag,
                        const Method::DenseElimination& m)
    {
        solve_precheck(x, A, b);

        commentator().start("Solve Modular DenseElimination for DenseMatrix", "solve.modular.dense-elimination.dense");

        LQUPMatrix<Field> LQUP(A);
        LQUP.left_solve(x, b);

        commentator().stop("solve.modular.dense-elimination.dense");

        return x;
    }

    /**
     * \brief Solve specialisation for DenseElimination on dense matrices with IntegerTag.
     */
    template <class ResultVector, class Field, class Vector>
    ResultVector& solve(ResultVector& x, const DenseMatrix<Field>& A, const Vector& b, const RingCategories::IntegerTag& tag,
                        const Method::DenseElimination& m)
    {
        // @fixme This is the original code for this case... but it goes for a Dixon!
        // Is that really what we want?
        return solve(x, A, b, tag, Method::Dixon(m));
    }
}
