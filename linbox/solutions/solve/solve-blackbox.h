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
     * \brief Solve specialisation for Blackbox.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Blackbox& m)
    {
        return solve(x, A, b, tag, reinterpret_cast<const Method::Wiedemann&>(m));
    }

    //
    // solveInPlace
    //

    /**
     * \brief Solve in place specialisation for Blackbox.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solveInPlace(ResultVector& x, Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Blackbox& m)
    {
        return solveInPlace(x, A, b, tag, reinterpret_cast<const Method::Wiedemann&>(m));
    }
}
