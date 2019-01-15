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
 *.
 */

#ifndef __LINBOX_solve_H
#define __LINBOX_solve_H

#include <iostream> // @fixme This is needed for givaro ring-interface to compile

#include <linbox/algorithms/cra-mpi.h> // @fixme Rename this file rational-remainder-distributed
#include <linbox/algorithms/rational-cra-full-multip.h>
#include <linbox/algorithms/rational-cra.h> // @fixme Rename this file rational-remainder
#include <linbox/field/field-traits.h>
#include <linbox/iterators/rebinder-solver.h>
#include <linbox/randiter/random-prime.h>
#include <linbox/solutions/hadamard-bound.h>
#include <linbox/solutions/methods.h>
#include <linbox/util/debug.h> // NotImplementedYet
#include <linbox/vector/mdr-vector.h>

namespace LinBox {

    /**
     * \brief Solve Ax = b, for x.
     *
     * Returns a vector x such that Ax = b.
     *
     * Specifically:
     * - A non singular: the unique solution is returned.
     * - A singular:
     *      - Consistent system: a random solution is returned.
     *          The method parameter can contain an hint that an
     *          arbitrary element of the solution space is acceptable instead,
     *          which can be faster to compute if one doesn't expect a result in that case.
     *      - Inconsistent system: zero vector is returned.
     *
     * CategoryTag is defaulted to FieldTraits<Matrix::Field>::categoryTag().
     *
     * SolveMethod is expected to be one of the following:
     * - Method::CRA
     * - Method::Hybrid
     * - Method::Elimination
     * - Method::SparseElimination
     *
     * @param [out] x solution, can be a rational solution (vector of numerators and one denominator)
     * @param [in]  A matrix
     * @param [in]  b target
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return reference to \p x
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag, class SolveMethod>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag, const SolveMethod& m)
    {
        throw NotImplementedYet("Solve specialisation is not implemented yet.");
    }

    /**
     * \brief Solve dispatcher for automated category tag.
     */
    template <class ResultVector, class Matrix, class Vector, class SolveMethod>
    inline ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const SolveMethod& m)
    {
        return solve(x, A, b, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief Solve dispatcher for automated solve method.
     */
    template <class ResultVector, class Matrix, class Vector>
    inline ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b)
    {
        return solve(x, A, b, Method::Hybrid());
    }

    // @fixme Split to files

    //
    // solve_cra
    //

    namespace {
        template <class Vector, class Element, class Iteration, class PrimeGenerator, class DispatchType>
        inline void solve_cra(Vector& num, Element& den, double hadamardLogBound, Iteration& iteration,
                              PrimeGenerator& primeGenerator, const DispatchType& dispatch)
        {
            throw NotImplementedYet("Integer CRA Solve with specified dispatch type is not implemented yet.");
        }

        template <class Vector, class Element, class Iteration, class PrimeGenerator>
        inline void solve_cra(Vector& num, Element& den, double hadamardLogBound, Iteration& iteration,
                              PrimeGenerator& primeGenerator, const Dispatch::None& dispatch)
        {
            using CraAlgorithm = FullMultipRatCRA<Givaro::ModularBalanced<double>>;
            RationalRemainder<CraAlgorithm> cra(hadamardLogBound);
            cra(num, den, iteration, primeGenerator);
        }

        template <class Vector, class Element, class Iteration, class PrimeGenerator>
        inline void solve_cra(Vector& num, Element& den, double hadamardLogBound, Iteration& iteration,
                              PrimeGenerator& primeGenerator, const Dispatch::Distributed& dispatch)
        {
            using CraAlgorithm = FullMultipRatCRA<Givaro::ModularBalanced<double>>;
            // @fixme Rename RationalRemainderDistributed
            MPIratChineseRemainder<CraAlgorithm> cra(hadamardLogBound, dispatch.communicator);
            cra(num, den, iteration, primeGenerator);
        }
    }

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
        solve_cra(num, den, hadamardLogBound, iteration, genprime, m.dispatch);

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

    //
    // solvein
    //
    // Solve that can modify the input matrix
    //

    /**
     * \brief Solve Ax = b, for x.
     *
     * Returns a vector x such that Ax = b.
     * A will be modified.
     *
     * See comment for `solve`.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag, class SolveMethod>
    inline ResultVector& solvein(ResultVector& x, Matrix& A, const Vector& b, const CategoryTag& tag, const SolveMethod& m)
    {
        // @fixme Should solvein be the entry point for everything and solve just copy the matrix before executing?

        // Here, we just fallback to non in-place solve
        return solve(x, A, b, tag, m);
    }

    //
    // Precheck
    //

    template <class ResultVector, class Matrix, class Vector>
    inline void solve_precheck(const ResultVector& x, const Matrix& A, const Vector& b)
    {
        if ((A.coldim() != x.size()) || (A.rowdim() != b.size())) {
            throw LinboxError("Incompatible dimensions of data in 'solve'.");
        }
    }
}

#endif
