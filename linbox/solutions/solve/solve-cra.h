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

/**
 * This file contains everything concerning solve using CRA method.
 *
 * CRA is structured to be able to dispatch each solve mod p
 * between multiple nodes or threads. This is one of the template argument of CRATraits.
 *
 * The other argument called IterationMethod is the method to use for
 * solving each mod p.
 */

#pragma once

#include <linbox/algorithms/cra-distributed.h>
#include <linbox/algorithms/rational-cra-builder-early-multip.h>
#include <linbox/algorithms/rational-cra-builder-full-multip.h>
#include <linbox/algorithms/rational-cra.h>
#include <linbox/field/rebind.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/solutions/hadamard-bound.h>
#include <linbox/util/commentator.h>
#include <linbox/util/debug.h> // NotImplementedYet
#include <linbox/vector/vector-traits.h>

namespace {
    /**
     * Initialized with a matrix A, vector b and a solve method,
     * returns x such that Ax = b (if any).
     *
     * The parentheses operator expects the computation to take place
     * in a rebinded environment, with another field.
     *
     * This is used as a proxy within solve integer CRA
     * to provide each solve operations mod p.
     */
    template <class Matrix, class Vector, class SolveMethod>
    struct CRASolveIteration {
        const Matrix& A;
        const Vector& b;
        const SolveMethod& m;

        CRASolveIteration(const Matrix& _A, const Vector& _b, const SolveMethod& _m)
            : A(_A)
            , b(_b)
            , m(_m)
        {
        }

        template <typename Field>
        typename LinBox::Rebind<Vector, Field>::other& operator()(typename LinBox::Rebind<Vector, Field>::other& x,
                                                                  const Field& F) const
        {
            using FMatrix = typename LinBox::Rebind<Matrix, Field>::other;
            using FVector = typename LinBox::Rebind<Vector, Field>::other;

            FMatrix FA(A, F);
            FVector Fb(F, b);

            LinBox::VectorWrapper::ensureDim(x, A.coldim());
            return solve(x, FA, Fb, m);
        }
    };

    template <class CRAField, class MatrixCategoryTag>
    struct BestCRABuilder {
        using type = LinBox::RationalCRABuilderFullMultip<CRAField>;
    };

    template <class CRAField>
    struct BestCRABuilder<CRAField, LinBox::RingCategories::RationalTag> {
        using type = LinBox::RationalCRABuilderEarlyMultip<CRAField>;
    };
}

namespace LinBox {
    /**
     * \brief Solve specialization with Chinese Remainder Algorithm method for an Integer or Rational tags.
     *
     * If a Dispatch::Distributed is used, please note that the result will only be set on the master node.
     */
    template <class IntVector, class Matrix, class Vector, class IterationMethod>
    inline void solve(IntVector& xNum, typename IntVector::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const Method::CRA<IterationMethod>& m)
    {
        //
        // Handle auto-dispatch.
        //

        Dispatch dispatch = m.dispatch;
        if (dispatch == Dispatch::Auto) {
            Method::CRA<IterationMethod> newM(m);

#if __LINBOX_HAVE_MPI
            // User has MPI enabled in config, but not specified if it wanted to use it,
            // we enable it with default communicator if needed.
            newM.dispatch = Dispatch::Distributed;
#else
            // @note Currently, Sequential = SMP if OpenMP is active on
            // the machine, as CRA -> ChineseRemainderSequential or ChineseRemainderOMP in cra-domain.h
            newM.dispatch = Dispatch::Sequential;
#endif

            return solve(xNum, xDen, A, b, tag, newM);
        }

        //
        // Declare communicator if none was yet.
        //

        if (m.dispatch == Dispatch::Distributed && m.pCommunicator == nullptr) {
            Method::CRA<IterationMethod> newM(m);
            Communicator communicator(nullptr, 0);
            newM.pCommunicator = &communicator;
            return solve(xNum, xDen, A, b, tag, newM);
        }

        //
        // Init all (Hadamard bound, prime generator).
        //

        if (m.master()) {
            commentator().start("solve.cra.integer");
            linbox_check((A.coldim() == xNum.size()) && (A.rowdim() == b.size()));
        }

        using CRAField = Givaro::ModularBalanced<double>;
        unsigned int bits = FieldTraits<CRAField>::bestBitSize(A.coldim());
        PrimeIterator<LinBox::IteratorCategories::HeuristicTag> primeGenerator(bits);
        CRASolveIteration<Matrix, Vector, IterationMethod> iteration(A, b, m.iterationMethod);

        // @note The result is stored to Integers, and will be converted
        // later back.
        using Ring = Givaro::ZRing<Integer>;
        Ring R;
        BlasVector<Ring> num(R, A.coldim());
        Integer den(1);

        //
        // Calling the right solver
        //

        // @fixme We would like to be able to choose between Early or Full at runtime,
        // but we can't yet.
        // @fixme Moreover, RationalSolveHadamardBound(A, b) for a rational matrix is not implemented yet...
        using MatrixCategoryTag = typename FieldTraits<typename Matrix::Field>::categoryTag;
        double hadamardLogBound = LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD;
        if (!std::is_same<MatrixCategoryTag, RingCategories::RationalTag>::value) {
            hadamardLogBound = RationalSolveHadamardBound(A, b).solutionLogBound;
        }

        using CRAAlgorithm = typename BestCRABuilder<CRAField, MatrixCategoryTag>::type;
        if (dispatch == Dispatch::Sequential) {
            LinBox::RationalChineseRemainder<CRAAlgorithm> cra(hadamardLogBound);
            cra(num, den, iteration, primeGenerator);
        }
#if defined(__LINBOX_HAVE_MPI)
        else if (dispatch == Dispatch::Distributed) {
            LinBox::ChineseRemainderDistributed<CRAAlgorithm> cra(hadamardLogBound, m.pCommunicator);
            cra(num, den, iteration, primeGenerator);
        }
#endif
        else {
            throw LinBox::NotImplementedYet("Integer CRA Solve with specified dispatch type is not implemented yet.");
        }

        //
        // Post-solve conversion.
        //

        // @note We need to convert because the storage might be some fixed-size integer types.
        if (m.master()) {
            typename IntVector::Field R = xNum.field();

            auto it_x = xNum.begin();
            auto it_num = num.begin();

            // Convert the result back
            for (; it_x != xNum.end(); ++it_x, ++it_num) {
                R.init(*it_x, *it_num);
            }
            R.init(xDen, den);

            // @note During Dispatch::Distributed, we do not dispatch the result to all other nodes,
            // to prevent unnecessary broadcast, as the doc says.

            commentator().stop("solve.cra.integer");
        }
    }

    /**
     * \brief Solve specialization with Chinese Remainder Algorithm method for a Rational matrix.
     */
    template <class RatVector, class RatMatrix, class IterationMethod>
    inline RatVector& solve(RatVector& x, const RatMatrix& A, const RatVector& b,
                            const RingCategories::RationalTag& tag, const Method::CRA<IterationMethod>& m)
    {
        using Ring = Givaro::ZRing<Integer>;

        Ring ZZ;
        Ring::Element xDen;
        DenseVector<Ring> xNum(ZZ, A.coldim());

        solve(xNum, xDen, A, b, RingCategories::IntegerTag(), m);

        // Copy result back to RatVector
        auto iXNum = xNum.begin();
        for (auto iX = x.begin(); iX != x.end(); ++iX) {
            *iX = typename RatVector::value_type(*iXNum, xDen);
            ++iXNum;
        }

        return x;
    }
}
