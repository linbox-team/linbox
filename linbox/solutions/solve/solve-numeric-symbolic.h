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

#include <linbox/algorithms/numeric-solver-lapack.h>
#include <linbox/algorithms/rational-solver-sn.h>
#include <linbox/field/field-traits.h>
#include <linbox/field/param-fuzzy.h>
#include <linbox/solutions/methods-wip.h>
// #include <linbox/algorithms/numsym.h>

namespace LinBox {
    //
    // Numeric Symbolic Overlap
    // Youse's variant of Numeric Symbolic rational solve.
    //

    template <class Matrix, class Vector>
    inline void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const MethodWIP::NumericSymbolicOverlap& m)
    {
        throw LinBoxError("Rational solve with MethodWIP::NumericSymbolicOverlap only works with DenseMatrix.");
    }

    /**
     * \brief Solve specialisation for NumericSymbolicOverlap with IntegerTag on DenseMatrix.
     */
    template <class Ring, class Vector>
    inline void solve(Vector& xNum, typename Ring::Element& xDen, const DenseMatrix<Ring>& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const MethodWIP::NumericSymbolicOverlap& m)
    {
        commentator().start("solve.numeric-symbolic-overlap.integer");
        linbox_check((A.coldim() != x.size()) || (A.rowdim() != b.size()));

        using Field = ParamFuzzy;
        using FMatrix = BlasMatrix<Field>;
        using NumericSolver = LPS<FMatrix>;

        bool e = false; // @fixme Useless?
        NumericSolver numSolver;
        RationalSolverSN<Ring, NumericSolver> rsolver(b.field(), numSolver, e);

        int status = rsolver.solve(xNum, xDen, A, b);
        if (status) { // @fixme SNSolverReturnStatus is more precise than that
            throw LinboxMathInconsistentSystem("From NumericSymbolicOverlap solve.");
        }

        commentator().stop("solve.numeric-symbolic-overlap.integer");
    }

    //
    // Numeric Symbolic Norm
    // Wan's variant of Numeric Symbolic rational solve.
    //

    template <class Matrix, class Vector>
    inline void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const MethodWIP::NumericSymbolicNorm& m)
    {
        throw LinBoxError("Rational solve with MethodWIP::NumericSymbolicNorm only works with DenseMatrix.");
    }

    /**
     * \brief Solve specialisation for NumericSymbolicNorm with IntegerTag on DenseMatrix.
     */
    template <class Ring, class Vector>
    inline void solve(Vector& xNum, typename Ring::Element& xDen, const DenseMatrix<Ring>& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const MethodWIP::NumericSymbolicNorm& m)
    {
        commentator().start("solve.numeric-symbolic-norm.integer");
        linbox_check((A.coldim() != x.size()) || (A.rowdim() != b.size()));

        using Field = Givaro::Modular<int32_t>; // @fixme Why not double?
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;

        RationalSolver<Ring, Field, PrimeGenerator, Method::NumSymNorm> rsolver(b.field());

        int status = rsolver.solve(xNum, xDen, A, b);
        if (status) {
            throw LinboxMathInconsistentSystem("From NumericSymbolicNorm solve.");
        }

        commentator().stop("solve.numeric-symbolic-norm.integer");
    }
}
