/*
 * Copyright (C) 2019 LinBox Team
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#pragma once

#include "./rational-cra-builder-full-multip.h"

namespace LinBox {
    /**
     * From a MultiModLiftingContainer, will build
     * the solution on each prime, then will do a CRT reconstruction,
     * before reconstructing the rational.95
     *
     * This does not do early termination.
     */
    template <class LiftingContainer>
    class MultiModRationalReconstruction {
        using Ring = typename LiftingContainer::Ring;
        using IElement = typename LiftingContainer::IElement;
        using IVector = typename LiftingContainer::IVector;
        using FElement = typename LiftingContainer::FElement;
        using FVector = typename LiftingContainer::FVector;

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

            commentator().start("[MultiModLifting] Lifting");

            // Temporary structure to store a ci for each pj
            std::vector<FVector> digits;
            digits.reserve(_lc.primesCount());
            for (auto& F : _lc.primesFields()) {
                digits.emplace_back(F, _lc.size());
            }

            // The pj^i for each pj
            std::vector<IElement> radices(_lc.primesCount(), 1);

            // Stores each c0 + c1 pj + ... + ck pj^k for each pj
            std::vector<IVector> padicAccumulations(_lc.primesCount(), _lc.ring());
            for (auto j = 0u; j < _lc.primesCount(); ++j) {
                padicAccumulations[j].resize(_lc.size());
            }

            // @fixme Better use PolEval (or will it cause memory explosion?)
            VectorDomain<Ring> IVD(_lc.ring());
            for (auto i = 0u; i < _lc.length(); ++i) {
                _lc.next(digits);

				#pragma omp parallel for
                for (auto j = 0u; j < _lc.primesCount(); ++j) {
                    // @fixme @cpernet digits being a field vector, this will implicitly cast
                    // each of its elements to a Integer, is there something better?
                    // Or else, we just need an overload of Givaro::ZRing().axpyin() with a double as last parameter
                    IVD.axpyin(padicAccumulations[j], radices[j], digits[j]); // y <- y + p^i * ci
                    _lc.ring().mulin(radices[j], _lc.prime(j));
                }
            }
            commentator().stop("[MultiModLifting] Lifting");

            // CRT reconstruction from paddicAccumulations
            commentator().start("[MultiModLifting] CRT Reconstruction Progress");
            using CRAField = Givaro::Modular<Integer>;
            RationalCRABuilderFullMultip<CRAField> craBuilder(_lc.log2Bound() / 1.4427); // 1.4427 = 1 / log(2)

            {
                CRAField field(radices[0]);
                craBuilder.initialize(field, padicAccumulations[0]);
            }

            for (auto j = 1u; j < _lc.primesCount(); ++j) {
                CRAField field(radices[j]);
                craBuilder.progress(field, padicAccumulations[j]);
            }
            commentator().stop("[MultiModLifting] CRT Reconstruction Progress");

            // Rational reconstruction
            craBuilder.result(xNum, xDen, _lc.numBound());

            return true;
        }

    private:
        LiftingContainer& _lc;
    };
}
