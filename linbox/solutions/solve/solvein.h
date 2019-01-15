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

namespace LinBox {
    /**
     * \brief Solve Ax = b, for x.
     *
     * Returns a vector x such that Ax = b.
     * A will be modified.
     *
     * See documentation for `solve`.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag, class SolveMethod>
    inline ResultVector& solvein(ResultVector& x, Matrix& A, const Vector& b, const CategoryTag& tag, const SolveMethod& m)
    {
        // @fixme Should solvein be the entry point for everything and solve just copy the matrix before executing?

        // Here, we just fallback to non in-place solve
        return solve(x, A, b, tag, m);
    }
}
