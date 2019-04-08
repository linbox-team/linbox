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

#include <linbox/algorithms/lanczos.h>
#include <linbox/algorithms/mg-block-lanczos.h>

namespace LinBox {
    //
    // Lanczos
    //

    /**
     * \brief Solve specialisation for Lanczos.
     */
    template <class Matrix, class Vector, class CategoryTag>
    Vector& solve(Vector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const Method::Lanczos& m)
    {
        throw LinboxError("Method::Lanczos can only be used with RingCategories::ModularTag.");
    }

    /**
     * \brief Solve specialisation for Lanczos with ModularTag.
     */
    template <class Matrix, class Vector>
    Vector& solve(Vector& x, const Matrix& A, const Vector& b, const RingCategories::ModularTag& tag, const Method::Lanczos& m)
    {
        using Solver = LanczosSolver<typename Matrix::Field, Vector>;
        Solver solver(A.field(), m);

        return solver.solve(A, x, b);
    }

    //
    // BlockLanczos
    //

    /**
     * \brief Solve specialisation for BlockLanczos.
     */
    template <class Matrix, class Vector, class CategoryTag>
    Vector& solve(Vector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const Method::BlockLanczos& m)
    {
        throw LinboxError("Method::BlockLanczos can only be used with RingCategories::ModularTag.");
    }

    /**
     * \brief Solve specialisation for BlockLanczos with ModularTag.
     */
    template <class Matrix, class Vector>
    Vector& solve(Vector& x, const Matrix& A, const Vector& b, const RingCategories::ModularTag& tag,
                  const Method::BlockLanczos& m)
    {
        using Solver = MGBlockLanczosSolver<typename Matrix::Field>;
        Solver solver(A.field(), m);

        bool solveResult = solver.solve(A, x, b);

        if (!solveResult) {
            throw LinboxMathInconsistentSystem("From BlockLanczos solve.");
        }

        return x;
    }
}
