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

namespace LinBox {
    /**
     * From a MultiModLiftingContainer, will build
     * the solution on each prime, then will do a CRT reconstruction,
     * before reconstructing the rational.
     *
     * This does not do early termination.
     */
    template <class LiftingContainer>
    class MultiModRationalReconstruction {
        using Ring = typename LiftingContainer::Ring;
        using IElement = typename LiftingContainer::IElement;
        using IVector = typename LiftingContainer::IVector;

    public:
        MultiModRationalReconstruction(LiftingContainer& lc)
            : _lc(lc)
        {
        }

        bool getRational(IVector& xNum, IElement& xDen)
        {
            // Early out when the numerator is bounded by zero.
            if (_lc.numBound() == 0) {
                for (auto i = 0u; i < _lc.length(); ++i) {
                    _lc.ring().assign(xNum[i], _lc.ring().zero);
                }
                _lc.ring().assign(xDen, _lc.ring().one);
                return true;
            }

            VectorDomain<Ring> IVD(_lc.ring());

            // Stores each c0 + c1 pj + ... + ck pj^k for each pj
            std::vector<IVector> padicAccumulations(_lc.primesCount(), _lc.ring());
            // Temporary structure to store a ci for each pj
            std::vector<IVector> digits(_lc.primesCount(),
                                        _lc.ring()); // @fixme Could be a Field Element?
            // The pj^i for each pj
            std::vector<IElement> radices(_lc.primesCount(), 1);

            for (auto j = 0u; j < _lc.primesCount(); ++j) {
                padicAccumulations[j].resize(_lc.size());
                digits[j].resize(_lc.size());
            }

            for (auto i = 0u; i < _lc.length(); ++i) {
                _lc.next(digits);

                // @fixme Better use PolEval (except memory explosion?)
                for (auto j = 0u; j < _lc.primesCount(); ++j) {
                    IVD.axpyin(padicAccumulations[j], radices[j], digits[j]); // y <- y + p^i * ci
                    _lc.ring().mulin(radices[j], _lc.prime(j));
                    auto xxx = (_lc._A.getEntry(0, 0) * padicAccumulations[j][0] - _lc._b[0]) % radices[j];
                    // std::cout << "xxx " << j << "." << i << " " << _lc._A.getEntry(0, 0) << " * " << padicAccumulations[j][0] << " - " << _lc._b[0] << " mod " << radices[j] << std::endl;
                    std::cout << "xxx " << j << "." << i << " " << xxx << std::endl;
                }
            }

            // CRT reconstruction from paddicAccumulations
            using CRAField = Givaro::Modular<Integer>;
            RationalCRABuilderFullMultip<CRAField> craBuilder(_lc.log2Bound()
                                                              / 1.4427); // 1.4427 = 1 / log(2)

            {
                CRAField field(radices[0]);
                craBuilder.initialize(field, padicAccumulations[0]);
            }

            for (auto j = 1u; j < _lc.primesCount(); ++j) {
                CRAField field(radices[j]);
                craBuilder.progress(field, padicAccumulations[j]);
            }

            // Rational reconstruction
            // @note RR expects the bounds to be strict, this is why we add a + 1
            craBuilder.result(xNum, xDen, _lc.numBound() + 1, _lc.denBound() + 1);

            return true;
        }

    private:
        LiftingContainer& _lc;
    };

    // @fixme Move that to a file - and make it be a RationalSolver<Method::DixonRNS>
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
        void solve(RVector& xNum, typename RVector::Element& xDen, const DenseMatrix<Ring>& A,
                   const Vector& b, const Method::DixonRNS& m)
        {
            // @fixme We should use some code from DixonSolver...
            // But that's hard so we just assume that A is square and invertible.
            linbox_check(A.rowdim() == A.coldim());

            using LiftingContainer = MultiModLiftingContainer<Field, Ring, PrimeGenerator>;
            LiftingContainer lc(_ring, _primeGenerator, A, b, m);
            MultiModRationalReconstruction<LiftingContainer> re(lc);

            if (!re.getRational(xNum, xDen)) {
                std::cerr << "OUCH!" << std::endl;
            }

            // #ifdef DEBUG_HADAMARD_BOUND
            std::clog << "numLog " << Givaro::logtwo(Givaro::abs(xNum[0])) << " " << xNum[0] << ';' << std::endl;
            std::clog << "denLog " << Givaro::logtwo(xDen) << " " << xDen << ';' << std::endl;
            // #endif
        }

    private:
        const Ring& _ring;
        PrimeGenerator& _primeGenerator;
    };

    /**
     * \brief Solve specialisation for DixonRNS on dense matrices.
     */
    template <class RVector, class Ring, class Vector>
    void solve(RVector& xNum, typename RVector::Element& xDen, const DenseMatrix<Ring>& A,
               const Vector& b, const RingCategories::IntegerTag& tag, const Method::DixonRNS& m)
    {
        commentator().start("solve.dixon-rns.integer.dense");

        // @fixme We don't know if we can use ModularBalanced<double>,
        // because of the rational reconstruction which might be
        // implicitly requiring 0-{p-1} representation of the p-adic sequence elements.
        using Field = Givaro::Modular<double>;
        using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
        PrimeGenerator primeGenerator(FieldTraits<Field>::bestBitSize(A.coldim()), 12); // @fixme REMOVE SEED

        DixonRNSSolver<Field, Ring, PrimeGenerator> solver(A.field(), primeGenerator);
        solver.solve(xNum, xDen, A, b, m);

        commentator().stop("solve.dixon-rns.integer.dense");

        // @fixme Implement something like that
        // if (status == SS_INCONSISTENT) {
        //     throw LinboxMathInconsistentSystem("From DixonRNS method.");
        // } else if (status == SS_FAILED || status == SS_BAD_PRECONDITIONER) {
        //     throw LinboxError("From DixonRNS method.");
        // }
    }
}