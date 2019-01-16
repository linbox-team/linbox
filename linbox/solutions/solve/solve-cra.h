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

#include <linbox/algorithms/cra-mpi.h> // @fixme Rename this file rational-remainder-distributed
#include <linbox/algorithms/rational-cra-full-multip.h>
#include <linbox/algorithms/rational-cra.h> // @fixme Rename this file rational-remainder
#include <linbox/iterators/rebinder-solver.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/solutions/hadamard-bound.h>
#include <linbox/util/commentator.h>

namespace {
    /**
     * Simple function to switch to a different algorithm according to the dispatch type.
     */
    template <class Vector, class Element, class Iteration, class PrimeGenerator, class DispatchType>
    inline void solve_from_dispatch(Vector& num, Element& den, double hadamardLogBound, Iteration& iteration,
                                    PrimeGenerator& primeGenerator, const DispatchType& dispatch)
    {
        throw LinBox::NotImplementedYet("Integer CRA Solve with specified dispatch type is not implemented yet.");
    }

    template <class Vector, class Element, class Iteration, class PrimeGenerator>
    inline void solve_from_dispatch(Vector& num, Element& den, double hadamardLogBound, Iteration& iteration,
                                    PrimeGenerator& primeGenerator, const LinBox::Dispatch::None& dispatch)
    {
        using CraAlgorithm = LinBox::FullMultipRatCRA<Givaro::ModularBalanced<double>>;
        LinBox::RationalRemainder<CraAlgorithm> cra(hadamardLogBound);
        cra(num, den, iteration, primeGenerator);
    }

#if defined(__LINBOX_HAVE_MPI) // @fixme Is this useful?
    template <class Vector, class Element, class Iteration, class PrimeGenerator>
    inline void solve_from_dispatch(Vector& num, Element& den, double hadamardLogBound, Iteration& iteration,
                                    PrimeGenerator& primeGenerator, const LinBox::Dispatch::Distributed& dispatch)
    {
        using CraAlgorithm = LinBox::FullMultipRatCRA<Givaro::ModularBalanced<double>>;
        // @fixme Rename RationalRemainderDistributed
        LinBox::MPIratChineseRemainder<CraAlgorithm> cra(hadamardLogBound, dispatch.communicator);
        cra(num, den, iteration, primeGenerator);
    }
#endif
}

namespace LinBox {
    /**
     * \brief Solve specialization with Chinese Remainder Algorithm method for an Integer ring.
     */
    template <class Matrix, class Vector, class IterationMethod, class DispatchType>
    inline MdrVector<Vector>& solve(MdrVector<Vector>& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag,
                                    const Method::CRAWIP<IterationMethod, DispatchType>& m)
    {
        if (m.dispatch.master()) {
            commentator().start("Solve Integer CRA", "solve.integer.cra");
            solve_precheck(x, A, b);
        }

        const typename Matrix::Field& F = A.field();
        unsigned int bits = 26 - (int)ceil(log(A.rowdim() * 0.7213475205));
        PrimeIterator<LinBox::IteratorCategories::HeuristicTag> genprime(bits);
        BlasVector<Givaro::ZRing<Integer>> num(F, A.coldim());
        Integer den(1);
        RebinderSolver<Matrix, Vector, IterationMethod> iteration(A, b, m.iterationMethod);

        auto hb = RationalSolveHadamardBound(A, b);
        double hadamardLogBound = 1.0 + hb.numLogBound + hb.denLogBound; // = ln2(2 * N * D)

        // @fixme This could take a MdrVector, right?
        solve_from_dispatch(num, den, hadamardLogBound, iteration, genprime, m.dispatch);

        // @note We need to convert because the storage might be some fixed-size integer types.
        if (m.dispatch.master()) {
            auto it_x = x.num.begin();
            auto it_num = num.begin();

            // Convert the result back to, what??
            for (; it_x != x.num.end(); ++it_x, ++it_num) {
                F.init(*it_x, *it_num);
            }
            F.init(x.den, den);

            // @fixme Should we synchronize x on all nodes?

            commentator().stop("solve.integer.cra");
        }

        return x;
    }

    /**
     * \brief Solve specialization with Chinese Remainder Algorithm method for an Integer ring, when dispatch is automated.
     */
    template <class Matrix, class Vector, class IterationMethod>
    inline MdrVector<Vector>& solve(MdrVector<Vector>& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag,
                                    const Method::CRAWIP<IterationMethod, Dispatch::Auto>& m)
    {
#if __LINBOX_HAVE_MPI
        // User has MPI enabled in config, but not specified if it wanted to use it,
        // we enable it with default communicator.
        Communicator communicator(nullptr, 0);

        Method::CRAWIP<IterationMethod, Dispatch::Distributed> newM;
        newM.iterationMethod = m.iterationMethod;
        newM.dispatch.communicator = &communicator;
#else
        // @fixme Should we use Dispatch::Threaded by default?
        Method::CRAWIP<IterationMethod, Dispatch::None> newM;
        newM.iterationMethod = m.iterationMethod;
#endif

        return solve(x, A, b, tag, newM);
    }
}
