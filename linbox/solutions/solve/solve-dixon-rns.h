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

#include <linbox/algorithms/multi-mod-lifting-container.h>
#include <linbox/algorithms/multi-mod-rational-reconstruction.h>

namespace LinBox {
    template <class Field, class Ring, class PrimeGenerator>
    class DixonRNSSolver {
    public:
        DixonRNSSolver(const Ring& ring, PrimeGenerator& primeGenerator)
            : _ring(ring)
            , _primeGenerator(primeGenerator)
        {
            /* @todo */
        }

        /**
         * Dense solving.
         */
        template <class RVector, class Vector>
        void solve(RVector& xNum, typename RVector::Element& xDen, const DenseMatrix<Ring>& A, const Vector& b,
                   const Method::DixonRNS& m)
        {
            // @fixme We should use some code from DixonSolver...
            // But that's hard so we just assume that A is square and invertible.
            linbox_check(A.rowdim() == A.coldim());

            using LiftingContainer = MultiModLiftingContainer<Field, Ring, PrimeGenerator>;

            commentator().start("[MultiModLifting] Init");
            LiftingContainer lc(_ring, _primeGenerator, A, b, m);
            MultiModRationalReconstruction<LiftingContainer> re(lc);
            commentator().stop("[MultiModLifting] Init");

            if (!re.getRational(xNum, xDen)) {
                std::cerr << "OUCH!" << std::endl;
            }
        }

    private:
        const Ring& _ring;
        PrimeGenerator& _primeGenerator;
    };

    /**
     * \brief Solve specialisation for DixonRNS on dense matrices.
     */
    template <class RVector, class Ring, class Vector>
    void solve(RVector& xNum, typename RVector::Element& xDen, const DenseMatrix<Ring>& A, const Vector& b,
               const RingCategories::IntegerTag& tag, const Method::DixonRNS& m)
    {
        commentator().start("solve.dixon-rns.integer.dense");

        // @fixme We don't know if we can use ModularBalanced<double>,
        // because of the rational reconstruction which might be
        // implicitly requiring 0-{p-1} representation of the p-adic sequence elements.
        using Field = Givaro::Modular<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()));

        // @fixme TO BE REMOVED
        DixonRNSSolver<Field, Ring, PrimeGenerator> solver(A.field(), primeGenerator);
        solver.solve(xNum, xDen, A, b, m);

        DixonSolver<Ring, Field, PrimeGenerator, Method::DenseElimination> classicSolver(A.field(), primeGenerator);
        classicSolver.solveNonsingular(xNum, xDen, A, b, false, m.trialsBeforeFailure);

        commentator().stop("solve.dixon-rns.integer.dense");

        // @fixme Implement something like that
        // if (status == SS_INCONSISTENT) {
        //     throw LinboxMathInconsistentSystem("From DixonRNS method.");
        // } else if (status == SS_FAILED || status == SS_BAD_PRECONDITIONER) {
        //     throw LinboxError("From DixonRNS method.");
        // }
    }
}