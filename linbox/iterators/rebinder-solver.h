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

#include <linbox/field/rebind.h>
#include <linbox/vector/vector-traits.h>

namespace LinBox {
    /**
     * Initialized with a matrix A, vector b and a solve method,
     * returns x such that Ax = b (if any).
     *
     * The parentheses operator excepts the computation to take place
     * in a rebinded environment, with another field.
     *
     * This is used as a proxy within solve integer CRA
     * to provide each solve operations mod p.
     */
    template <class Matrix, class Vector, class SolveMethod>
    struct RebinderSolver {
        const Matrix& A;
        const Vector& b;
        const SolveMethod& m;

        RebinderSolver(const Matrix& _A, const Vector& _b, const SolveMethod& _m)
            : A(_A)
            , b(_b)
            , m(_m)
        {
        }

        template <typename Field>
        typename Rebind<Vector, Field>::other& operator()(typename Rebind<Vector, Field>::other& x, const Field& F) const
        {
            using FMatrix = typename Matrix::template rebind<Field>::other;
            using FVector = typename Rebind<Vector, Field>::other;

            FMatrix FA(A, F);

            // @fixme Why is the convention reversed here?
            // Above it's (-, Field) whereas here is (Field, -).
            FVector Fb(F, b);

            VectorWrapper::ensureDim(x, A.coldim());
            return solve(x, FA, Fb, m);
        }
    };
}
