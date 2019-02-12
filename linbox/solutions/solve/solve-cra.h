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

#include <linbox/algorithms/cra-mpi.h> // @fixme Rename this file cra-rational-solver-distributed
#include <linbox/algorithms/rational-cra-full-multip.h>
#include <linbox/algorithms/rational-cra.h> // @fixme Rename this file cra-rational-solver
#include <linbox/field/rebind.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/solutions/hadamard-bound.h>
#include <linbox/util/commentator.h>
#include <linbox/vector/vector-traits.h>

namespace {
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
    struct CraRebinderSolver {
        const Matrix& A;
        const Vector& b;
        const SolveMethod& m;

        CraRebinderSolver(const Matrix& _A, const Vector& _b, const SolveMethod& _m)
            : A(_A)
            , b(_b)
            , m(_m)
        {
        }

        template <typename Field>
        typename LinBox::Rebind<Vector, Field>::other& operator()(typename LinBox::Rebind<Vector, Field>::other& x,
                                                                  const Field& F) const
        {
            using FMatrix = typename Matrix::template rebind<Field>::other;
            using FVector = typename LinBox::Rebind<Vector, Field>::other;

            FMatrix FA(A, F);

            // @fixme Why is the convention reversed here?
            // Above it's (-, Field) whereas here is (Field, -).
            FVector Fb(F, b);

            LinBox::VectorWrapper::ensureDim(x, A.coldim());
            return solve(x, FA, Fb, m);
        }
    };
}

namespace LinBox {
    /**
     * \brief Solve specialization with Chinese Remainder Algorithm method for an Integer ring.
     *
     * This one has a different signature (num, den) because of the IntegerTag implying a rational result.
     */
    template <class Matrix, class Vector, class IterationMethod>
    inline void solve(Vector& xNum, typename Vector::Field::Element& xDen, const Matrix& A, const Vector& b,
                      const RingCategories::IntegerTag& tag, const MethodWIP::Cra<IterationMethod>& m)
    {
        //
        // Handle auto-dispatch.
        //

        Dispatch dispatch = m.dispatch;
        if (dispatch == Dispatch::Auto) {
            MethodWIP::Cra<IterationMethod> newM(m);

            // User has MPI enabled in config, but not specified if it wanted to use it,
            // we enable it with default communicator if needed.
#if __LINBOX_HAVE_MPI
            newM.dispatch = Dispatch::Distributed;
#else
            newM.dispatch = Dispatch::Sequential; // @fixme Should we use Dispatch::Smp by default?
#endif

            return solve(xNum, xDen, A, b, tag, newM);
        }

        //
        // Declare communicator if none was yet.
        //

        if (m.dispatch == Dispatch::Distributed && m.pCommunicator == nullptr) {
            MethodWIP::Cra<IterationMethod> newM(m);
            Communicator communicator(nullptr, 0);
            newM.pCommunicator = &communicator;
            return solve(xNum, xDen, A, b, tag, newM);
        }

        //
        // Init all (Hadamard bound, prime generator).
        //

        if (m.master()) {
            commentator().start("solve.cra.integer");
            linbox_check((A.coldim() != xNum.size()) || (A.rowdim() != b.size()));
        }

        using CraField = Givaro::ModularBalanced<double>;
        const typename Matrix::Field& F = A.field();
        unsigned int bits = FieldTraits<CraField>::bestBitSize(A.coldim());
        PrimeIterator<LinBox::IteratorCategories::HeuristicTag> primeGenerator(bits);
        CraRebinderSolver<Matrix, Vector, IterationMethod> iteration(A, b, m.iterationMethod);

        // @note The result is stored to Integers, and will be converted
        // later back.
        BlasVector<Givaro::ZRing<Integer>> num(F, A.coldim());
        Integer den(1);

        auto hb = RationalSolveHadamardBound(A, b);
        double hadamardLogBound = 1.0 + hb.numLogBound + hb.denLogBound; // = ln2(2 * N * D)

        //
        // Calling the right solver
        //

        using CraAlgorithm = LinBox::FullMultipRatCRA<CraField>;
        if (dispatch == Dispatch::Sequential) {
            LinBox::RationalRemainder<CraAlgorithm> cra(hadamardLogBound);
            cra(num, den, iteration, primeGenerator);
        }
#if defined(__LINBOX_HAVE_MPI)
        else if (dispatch == Dispatch::Distributed) {
            // @fixme Rename RationalRemainderDistributed
            LinBox::MPIratChineseRemainder<CraAlgorithm> cra(hadamardLogBound, m.pCommunicator);
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
        // @fixme Specialize so that this copy isn't needed
        if (m.master()) {
            auto it_x = xNum.begin();
            auto it_num = num.begin();

            // Convert the result back to, what??
            for (; it_x != xNum.end(); ++it_x, ++it_num) {
                F.init(*it_x, *it_num);
            }
            F.init(xDen, den);

            // @fixme Should we synchronize x on all nodes?

            commentator().stop("solve.cra.integer");
        }
    }
}
