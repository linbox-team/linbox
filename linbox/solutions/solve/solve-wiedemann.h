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

#include <linbox/algorithms/block-wiedemann.h>
#include <linbox/solutions/methods.h>

namespace LinBox {
    /**
     * \brief Solve specialisation for Wiedemann.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Wiedemann& m)
    {
        solve_precheck(x, A, b);

        // @fixme Well, it did not existed before...

        return x;
    }

    /**
     * \brief Solve specialisation for BlockWiedemann.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag,
                        const Method::BlockWiedemann& m)
    {
        solve_precheck(x, A, b);

        // @fixme This does not work

        // using Context = BlasMatrixDomain<typename Matrix::Field>;
        // Context BMD(A.field());
        // BlockWiedemannSolver<Context> BWS(BMD, m.blockingFactor(), m.blockingFactor() + 1);
        // BWS.solve(x, A, b);

        return x;
    }
}
