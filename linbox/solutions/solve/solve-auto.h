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

#include <linbox/solutions/methods.h>

namespace LinBox {
    /**
     * \brief Solve specialisation for Auto.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Auto& m)
    {
        // @fixme useBB(), as the auto should go according to sparse or so
        return solve(x, A, b, tag, Method::Elimination());
    }

    /**
     * \brief Solve specialisation for Auto and IntegerTag.
     */
    template <class ResultVector, class Matrix, class Vector>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag,
                        const Method::Auto& m)
    {
        return solve(x, A, b, tag, Method::Dixon());
    }
}
