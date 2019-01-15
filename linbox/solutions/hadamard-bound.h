/*
 * Copyright (C) 2018
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

#include <linbox/field/field-traits.h>
#include <linbox/integer.h>
#include <linbox/matrix/matrix-category.h>
#include <linbox/matrix/matrix-traits.h>

#include <limits>

namespace LinBox {

    // ----- Vector norm

    template <class ConstIterator>
    double vectorLogNorm(const ConstIterator& begin, const ConstIterator& end)
    {
        Integer normSquared = 0;
        for (ConstIterator it = begin; it != end; ++it) {
            // Whatever field element it is,
            // it should be able to store the square without
            // loss of information.
            normSquared += (*it) * (*it);
        }

        if (normSquared == 0) {
            return 0.0;
        }

        return Givaro::logtwo(normSquared) / 2.0;
    }

    // ----- Detailed Hadamard bound

    struct HadamardLogBoundDetails {
        /**
         * Bit size of the minimal hadamard bound
         * between the row-wise and the col-wise ones.
         *
         * min { HadamardRow(A), HadamardCol(A) }
         */
        double logBound;
        /**
         * Bit size of the minimal hadamard bound
         * divided by the min of the norm vectors
         * between the row-wise and the col-wise ones.
         *
         * min { HadamardRow(A) / min || Ai,* ||,
         *       HadamardCol(A) / min || A*,j ||  }
         */
        double logBoundOverMinNorm;
    };

    /**
     * Precise Hadamard bound (bound on determinant) by taking
     * the row-wise euclidean norm.
     *
     * The result is expressed as bit size.
     */
    template <class IMatrix>
    void HadamardRowLogBound(double& logBound, double& minLogNorm, const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        HadamardRowLogBound(logBound, minLogNorm, A, tag);
    }

    template <class IMatrix>
    void HadamardRowLogBound(double& logBound, double& minLogNorm, const IMatrix& A, const MatrixCategories::RowColMatrixTag& tag)
    {
        logBound = 0.0;
        minLogNorm = std::numeric_limits<double>::infinity();

        for (auto rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            double rowLogNorm = vectorLogNorm(rowIt->begin(), rowIt->end());
            if (rowLogNorm == 0.0) {
                logBound = 0.0;
                minLogNorm = 0.0;
                return;
            }
            else if (rowLogNorm < minLogNorm) {
                minLogNorm = rowLogNorm;
            }
            logBound += rowLogNorm;
        }
    }

    template <class IMatrix>
    void HadamardRowLogBound(double& logBound, double& minLogNorm, const IMatrix& A, const MatrixCategories::RowMatrixTag& tag)
    {
        logBound = 0.0;
        minLogNorm = std::numeric_limits<double>::infinity();

        for (auto rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            Integer normSquared = 0;
            for (const auto& pair : *rowIt) {
                normSquared += (pair.second) * (pair.second);
            }
            if (normSquared == 0) {
                logBound = 0.0;
                minLogNorm = 0.0;
                return;
            }

            double logNormSquared = Givaro::logtwo(normSquared);
            if (logNormSquared < minLogNorm) {
                minLogNorm = logNormSquared;
            }
            logBound += logNormSquared;
        }

        // Square-root
        logBound /= 2.0;
        minLogNorm /= 2.0;
    }

    /**
     * Precise Hadamard bound (bound on determinant) by taking
     * the column-wise euclidean norm.
     *
     * The result is expressed as bit size.
     */
    template <class IMatrix>
    void HadamardColLogBound(double& logBound, double& minLogNorm, const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        HadamardColLogBound(logBound, minLogNorm, A, tag);
    }

    template <class IMatrix>
    void HadamardColLogBound(double& logBound, double& minLogNorm, const IMatrix& A, const MatrixCategories::RowColMatrixTag& tag)
    {
        logBound = 0.0;
        minLogNorm = std::numeric_limits<double>::infinity();

        typename IMatrix::ConstColIterator colIt;
        for (colIt = A.colBegin(); colIt != A.colEnd(); ++colIt) {
            double colLogNorm = vectorLogNorm(colIt->begin(), colIt->end());
            if (colLogNorm < minLogNorm) {
                minLogNorm = colLogNorm;
            }
            logBound += colLogNorm;
        }
    }

    template <class IMatrix>
    void HadamardColLogBound(double& logBound, double& minLogNorm, const IMatrix& A, const MatrixCategories::RowMatrixTag& tag)
    {
        logBound = 0.0;
        minLogNorm = std::numeric_limits<double>::infinity();

        // This vector contains the norm squared for each columns.
        std::vector<Integer> columnsNormsSquared(A.coldim());
        for (auto rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            for (const auto& pair : *rowIt) {
                columnsNormsSquared[pair.first] += (pair.second) * (pair.second);
            }
        }

        // All the norms have been computed, we check which one is the smallest
        // and compute the product (aka sum bitsize-wise) of them to make the logBound.
        for (const Integer& normSquared : columnsNormsSquared) {
            if (normSquared == 0) {
                logBound = 0.0;
                minLogNorm = 0.0;
                return;
            }
            double logNormSquared = Givaro::logtwo(normSquared);
            if (logNormSquared < minLogNorm) {
                minLogNorm = logNormSquared;
            }
            logBound += logNormSquared;
        }

        // Square-root
        logBound /= 2.0;
        minLogNorm /= 2.0;
    }

    /**
     * Precise Hadamard bound (bound on determinant) by taking the minimum
     * of the column-wise and the row-wise euclidean norm.
     *
     * The results are expressed as bit size.
     */
    template <class IMatrix>
    HadamardLogBoundDetails DetailedHadamardBound(const IMatrix& A)
    {
        double rowLogBound = 0;
        double rowMinLogNorm = 0;
        HadamardRowLogBound(rowLogBound, rowMinLogNorm, A);
        double rowLogBoundOverMinNorm = rowLogBound - rowMinLogNorm;

        double colLogBound = 0;
        double colMinLogNorm = 0;
        HadamardColLogBound(colLogBound, colMinLogNorm, A);
        double colLogBoundOverMinNorm = colLogBound - colMinLogNorm;

        HadamardLogBoundDetails data;
        data.logBound = std::min(rowLogBound, colLogBound) + 0.01;
        data.logBoundOverMinNorm = std::min(rowLogBoundOverMinNorm, colLogBoundOverMinNorm) + 0.01;
        return data;
    }

    // ----- Hadamard bound

    /**
     * Precise Hadamard bound (bound on determinant) by taking the minimum
     * of the column-wise and the row-wise euclidean norm.
     *
     * The result is expressed as bit size.
     */
    template <class IMatrix>
    double HadamardBound(const IMatrix& A)
    {
        return DetailedHadamardBound(A).logBound;
    }

    // ----- Fast Hadamard bound

    /**
     * Returns the bit size of the Hadamard bound.
     * This is a larger estimation but faster to compute.
     */
    template <class IMatrix>
    double FastHadamardBound(const IMatrix& A)
    {
        Integer max = 0;
        for (auto it = A.Begin(); it != A.End(); ++it) {
            const Integer& ai = *it;
            if (max < ai)
                max = ai;
            else if (max < -ai)
                max = -ai;
        }

        if (max == 0) {
            return 0;
        }

        auto n = std::max(A.rowdim(), A.coldim());
        double logBound = static_cast<double>(n) * (Givaro::logtwo(n) / 2.0 + Givaro::logtwo(max));

        return logBound;
    }

    // ----- Rational solve bound

    struct RationalSolveHadamardBoundData {
        double numLogBound;
        double denLogBound;
    };

    /**
     * Bound on the rational solution of a linear system Ax = b.
     *
     * Return bounds on the bit sizes of both denominator and numerator of the solution x.
     *
     * @note Matrix and Vector should be over Integer.
     */
    template <class IMatrix, class Vector>
    RationalSolveHadamardBoundData RationalSolveHadamardBound(const IMatrix& A, const Vector& b)
    {
        RationalSolveHadamardBoundData data;

        auto hadamardBound = DetailedHadamardBound(A);
        double bLogNorm = vectorLogNorm(b.begin(), b.end());

        data.numLogBound = hadamardBound.logBoundOverMinNorm + bLogNorm + 1;
        data.denLogBound = hadamardBound.logBound;
        return data;
    }
}
