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

    /**
     * \brief Solve specialization with Chinese Remainder Algorithm method for an Integer ring.
     */
    template <class Matrix, class Vector, class IterationMethod>
    inline MdrVector<Vector>& solve(MdrVector<Vector>& x, const Matrix& A, const Vector& b, const RingCategories::IntegerTag& tag,
                                    const Method::CRAWIP<IterationMethod>& m)
    {
#if __LINBOX_HAVE_MPI
        // User has MPI enabled in config, but not specified if it wanted to use it,
        // we enable it with default communicator.
        Communicator communicator;
        if (m.dispatch.type == DispatchType::Unknown) {
            m.dispatch.type = DispatchType::Distributed;
            m.dispatch.communicator = &communicator;
        }
#endif

        if (!m.dispatch.communicator || m.dispatch.communicator->rank() == 0) {
            commentator().start("Solve Integer CRA", "solve.integer.cra");
            solve_precheck(x, A, b);
        }

        // Unspecified method, use default dispatch
        if (m.dispatch.type == DispatchType::Unknown) {
            // @fixme Should we use threaded by default?
            m.dispatch.type = DispatchType::None;
        }

        // @fixme Implement
        const typename Matrix::Field& F = A.field();
        unsigned int bits = 26 - (int)ceil(log(A.rowdim() * 0.7213475205));
        PrimeIterator<LinBox::IteratorCategories::HeuristicTag> genprime(bits);
        BlasVector<Givaro::ZRing<Integer>> num(F, A.coldim());
        Integer den(1);
        RebinderSolver<Matrix, Vector, IterationMethod> iteration(A, b, m.iterationMethod);

        auto hb = RationalSolveHadamardBound(A, b);
        double hadamard = 1.0 + hb.numLogBound + hb.denLogBound; // = ln2(2 * N * D)

        // @fixme Use templates to prevent that if?
        using CraAlgorithm = FullMultipRatCRA<Givaro::ModularBalanced<double>>;
        if (m.dispatch.type == DispatchType::None) {
            RationalRemainder<CraAlgorithm> cra(hadamard);
            // @fixme This could take a MdrVector, right?
            cra(num, den, iteration, genprime);
        }
        else if (m.dispatch.type == DispatchType::Distributed) {
            // @fixme Rename RationalRemainderDistributed
            MPIChineseRemainder<CraAlgorithm> cra(hadamard, m.dispatch.communicator);
            cra(num, den, iteration, genprime);
        }
        else {
            throw NotImplementedYet("Integer CRA Solve with specified dispatch type is not implemented yet.");
        }

        // @note We need to convert because the storage might be some fixed-size integer types.
        if (!m.dispatch.communicator || m.dispatch.communicator->rank() == 0) {
            auto it_x = x.num.begin();
            auto it_num = num.begin();

            // Convert the result back to, what??
            for (; it_x != x.end(); ++it_x, ++it_num) {
                F.init(*it_x, *it_num);
            }
            F.init(x.den, den);

            commentator().stop("solve.integer.cra");
        }
    }

    //
    // solve_in
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
    inline ResultVector& solve_in(ResultVector& x, Matrix& A, const Vector& b, const CategoryTag& tag, const SolveMethod& m)
    {
        // @fixme Should solve_in be the entry point for everything and solve just copy the matrix before executing?

        // Here, we just fallback to non in-place solve
        return solve(x, A, b, tag, m);
    }

    //
    // Precheck
    //

    template <class ResultVector, class Matrix, class Vector, class CategoryTag, class SolveMethod>
    inline bool solve_precheck(const ResultVector& x, const Matrix& A, const Vector& b)
    {
        if ((A.coldim() != x.size()) || (A.rowdim() != b.size())) {
            throw LinboxError("Incompatible dimensions of data in 'solve'.");
        }
    }
}

#endif
