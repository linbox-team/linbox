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
    template <class Matrix, class Vector>
    Vector& solve(Vector& x, const Matrix& A, const Vector& b, const RingCategories::ModularTag& tag,
                  const Method::DenseElimination& m)
    {
        commentator().report(Commentator::LEVEL_UNIMPORTANT,
                             "Warning: Solve implicitly convert to a dense matrix because of Method::DenseElimination."
                             "This is usually expected behavior for small-size blackboxes.");

        // Copy the matrix into a dense one.
        DenseMatrix<typename Matrix::Field> ACopy(A);

        return solve(x, ACopy, b, tag, m);
    }

    /**
     * \brief Solve specialisation for DenseElimination on dense matrices with ModularTag.
     */
    template <class Field, class Vector>
    Vector& solve(Vector& x, const DenseMatrix<Field>& A, const Vector& b, const RingCategories::ModularTag& tag,
                  const Method::DenseElimination& m)
    {
        linbox_check((A.coldim() == x.size()) && (A.rowdim() == b.size()));

        commentator().start("solve.dense-elimination.modular.dense");

        PLUQMatrix<Field> PLUQ(A);
        PLUQ.left_solve(x, b);

        commentator().stop("solve.dense-elimination.modular.dense");

        return x;
    }
}
