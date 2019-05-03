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

#include <linbox/algorithms/dixon-rns-solver.h>

namespace LinBox {
    /**
     * \brief Solve specialisation for DixonRNS on dense matrices.
     */
    template <class IntVector, class Ring, class Vector>
    void solve(IntVector& xNum, typename IntVector::Element& xDen, const DenseMatrix<Ring>& A, const Vector& b,
               const RingCategories::IntegerTag& tag, const Method::DixonRNS& m)
    {
        commentator().start("solve.dixon.integer.dense");

        using Field = Givaro::ModularBalanced<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()));

        DixonRNSSolver<Field, Ring, PrimeGenerator> solver(A.field(), primeGenerator);
        solver.solve(xNum, xDen, A, b, m);

        commentator().stop("solve.dixon.integer.dense");

        // @fixme Implement something like that
        // if (status == SS_INCONSISTENT) {
        //     throw LinboxMathInconsistentSystem("From DixonRNS method.");
        // } else if (status == SS_FAILED || status == SS_BAD_PRECONDITIONER) {
        //     throw LinboxError("From DixonRNS method.");
        // }
    }
}