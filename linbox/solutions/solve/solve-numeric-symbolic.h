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
#include <linbox/solutions/methods.h>

namespace LinBox {
    //
    // Numeric Symbolic Overlap
    // Youse's variant of Numeric Symbolic rational solve.
    //

// @note LPS<...> requires LAPACK
#if defined(__FFLASFFPACK_HAVE_LAPACK)

    template <class IntVector, class Matrix, class Vector>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const Method::SymbolicNumericOverlap& m)
    {
        throw LinBoxError("Rational solve with Method::SymbolicNumericOverlap only works with DenseMatrix.");
    }

    /**
     * \brief Solve specialisation for SymbolicNumericOverlap with IntegerTag on DenseMatrix.
     */
    template <class IntVector, class Ring, class Vector>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const DenseMatrix<Ring>& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const Method::SymbolicNumericOverlap& m)
    {
        commentator().start("solve.numeric-symbolic-overlap.integer");
        linbox_check((A.coldim() == xNum.size()) && (A.rowdim() == b.size()));

        using Field = ParamFuzzy;
        using FMatrix = BlasMatrix<Field>;
        using NumericSolver = LPS<FMatrix>;

        NumericSolver numSolver;
        RationalSolverSN<Ring, NumericSolver> rsolver(b.field(), numSolver, false);

        SolverReturnStatus status = rsolver.solve(xNum, xDen, A, b);

        commentator().stop("solve.numeric-symbolic-overlap.integer");

        if (status == SS_INCONSISTENT) {
            throw LinboxMathInconsistentSystem("From SymbolicNumericOverlap solve.");
        }
        else if (status != SS_OK) {
            throw LinboxError("Failed to solve with SymbolicNumericOverlap.");
        }
    }

#endif

    //
    // Numeric Symbolic Norm
    // Wan's variant of Numeric Symbolic rational solve.
    //

    template <class IntVector, class Matrix, class Vector>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const Method::SymbolicNumericNorm& m)
    {
        throw LinBoxError("Rational solve with Method::SymbolicNumericNorm only works with DenseMatrix.");
    }

    /**
     * \brief Solve specialisation for SymbolicNumericNorm with IntegerTag on DenseMatrix.
     */
    template <class IntVector, class Ring, class Vector>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const DenseMatrix<Ring>& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const Method::SymbolicNumericNorm& m)
    {
        commentator().start("solve.numeric-symbolic-norm.integer");
        linbox_check((A.coldim() == xNum.size()) && (A.rowdim() == b.size()));

        using Field = Givaro::Modular<int32_t>; // @fixme Why not double?
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;

        DixonSolver<Ring, Field, PrimeGenerator, Method::SymbolicNumericNorm> rsolver(b.field());

        SolverReturnStatus status = rsolver.solve(xNum, xDen, A, b);

        commentator().stop("solve.numeric-symbolic-norm.integer");

        if (status == SS_INCONSISTENT) {
            throw LinboxMathInconsistentSystem("From SymbolicNumericNorm solve.");
        }
        else if (status != SS_OK) {
            throw LinboxError("Failed to solve with SymbolicNumericNorm.");
        }
    }
}
