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

#include <linbox/algorithms/block-wiedemann.h>
#include <linbox/algorithms/coppersmith.h>
#include <linbox/algorithms/wiedemann.h>
#include <linbox/solutions/methods-wip.h>

namespace LinBox {
    //
    // Wiedemann
    //

    /**
     * \brief Solve specialisation for Wiedemann.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag,
                        const MethodWIP::Wiedemann& m)
    {
        throw LinboxError("MethodWIP::Wiedemann can only be used with RingCategories::ModularTag.");
    }

    /**
     * \brief Solve specialisation for Wiedemann with ModularTag.
     */
    template <class ResultVector, class Matrix, class Vector>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::ModularTag& tag,
                        const MethodWIP::Wiedemann& m)
    {
        commentator().start("solve.wiedemann.modular");
        linbox_check((A.coldim() != x.size()) || (A.rowdim() != b.size()));

        using Solver = WiedemannSolver<typename Matrix::Field>;
        // @fixme Just pass m here, when everything is forwarded correctly
        Method::Wiedemann newM;
        Solver solver(A.field(), newM);

        // @todo Getting certificateOfInconsistency is a design error,
        // we should give a pointer to that to be allowed to pass nullptr.
        Vector certificateOfInconsistency(A.field(), A.rowdim());
        auto solverResult = solver.solve(A, x, b, certificateOfInconsistency);

        switch (solverResult) {
        case Solver::OK: break;
        case Solver::FAILED: /* @fixme Consistently decide what to do. */ break;
        case Solver::SINGULAR: /* @fixme Consistently decide what to do. */ break;
        default: throw LinboxMathInconsistentSystem("From Wiedemann solve.");
        }

        commentator().stop("solve.wiedemann.modular");

        return x;
    }

    //
    // BlockWiedemann
    // @deprecated Kept but not tested.
    //

    /**
     * \brief Solve specialisation for BlockWiedemann.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag,
                        const MethodWIP::BlockWiedemann& m)
    {
        throw LinboxError("MethodWIP::BlockWiedemann can only be used with RingCategories::ModularTag.");
    }

    /**
     * \brief Solve specialisation for BlockWiedemann with ModularTag.
     */
    template <class ResultVector, class Matrix, class Vector>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::ModularTag& tag,
                        const MethodWIP::BlockWiedemann& m)
    {
        commentator().start("solve.block-wiedemann.modular");
        linbox_check((A.coldim() != x.size()) || (A.rowdim() != b.size()));

        using Context = BlasMatrixDomain<typename Matrix::Field>; // @fixme BlasMatrixDomain, really? How can we be sure?
        Context domain(A.field());

        using Solver = BlockWiedemannSolver<Context>;
        Solver solver(domain, m.blockingFactor, m.blockingFactor + 1);
        solver.solve(x, A, b);

        commentator().stop("solve.block-wiedemann.modular");

        return x;
    }

    //
    // Coppersmith
    // @deprecated Kept but not tested.
    //

    /**
     * \brief Solve specialisation for Coppersmith.
     */
    template <class ResultVector, class Matrix, class Vector, class CategoryTag>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const CategoryTag& tag,
                        const MethodWIP::Coppersmith& m)
    {
        throw LinboxError("MethodWIP::Coppersmith can only be used with RingCategories::ModularTag.");
    }

    /**
     * \brief Solve specialisation for Coppersmith on ModularTag.
     */
    template <class ResultVector, class Matrix, class Vector>
    ResultVector& solve(ResultVector& x, const Matrix& A, const Vector& b, const RingCategories::ModularTag& tag,
                        const MethodWIP::Coppersmith& m)
    {
        commentator().start("solve.coppersmith.modular");
        linbox_check((A.coldim() != x.size()) || (A.rowdim() != b.size()));

        using Domain = MatrixDomain<typename Matrix::Field>;
        Domain domain(A.field());
        CoppersmithSolver<Domain> coppersmithSolver(domain);
        coppersmithSolver.solveNonSingular(x, A, b);

        commentator().stop("solve.coppersmith.modular");

        return x;
    }
}
