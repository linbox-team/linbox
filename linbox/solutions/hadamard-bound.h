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

#include <linbox/matrix/matrix-category.h>
#include <linbox/field/field-traits.h>
#include <linbox/integer.h>

namespace LinBox {

    // ----- Vector norm

    // @todo Specialize for Rationals
    template <class ConstIterator>
    size_t vectorNormBitSize(const ConstIterator& begin, const ConstIterator& end)
    {
        Integer norm = 0;
        for (ConstIterator it = begin; it != end; ++it) {
            // Whatever field element it is,
            // it should be able to store the square without
            // loss of information.
            norm += (*it) * (*it);
        }

        if (norm == 0) {
            return 0;
        }

        return std::ceil(Givaro::logtwo(norm) / 2.0);
    }

    // ----- Detailed Hadamard bound

    struct DetailedHadamardBoundData {
        /**
         * Bit size of the minimal hadamard bound
         * between the row-wise and the col-wise ones.
         *
         * min { HadamardRow(A), HadamardCol(A) }
         */
        size_t boundBitSize;
        /**
         * Bit size of the minimal hadamard bound
         * divided by the min of the norm vectors
         * between the row-wise and the col-wise ones.
         *
         * min { HadamardRow(A) / min || Ai,* ||,
         *       HadamardCol(A) / min || A*,j ||  }
         */
        size_t boundOnMinNormBitSize;
    };

    /**
     * Precise Hadamard bound (bound on determinant) by taking
     * the row-wise euclidean norm.
     *
     * The result is expressed as bit size.
     */
    template <class IMatrix>
    void DetailedHadamardRowBound(size_t& boundBitSize, size_t& minNormBitSize, const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        DetailedHadamardRowBound(boundBitSize, minNormBitSize, A, tag);
    }

    template <class IMatrix>
    void DetailedHadamardRowBound(size_t& boundBitSize, size_t& minNormBitSize, const IMatrix& A, const MatrixCategories::RowColMatrixTag& tag)
    {
        boundBitSize = 0;
        minNormBitSize = -1u;

        for (auto rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            size_t rowNormBitSize = vectorNormBitSize(rowIt->begin(), rowIt->end());
            if (rowNormBitSize < minNormBitSize) {
                minNormBitSize = rowNormBitSize;
            }
            boundBitSize += rowNormBitSize;
        }
    }

    template <class IMatrix>
    void DetailedHadamardRowBound(size_t& boundBitSize, size_t& minNormBitSize, const IMatrix& A, const MatrixCategories::RowMatrixTag& tag)
    {
        boundBitSize = 0;
        minNormBitSize = -1u;

        // @fixme Could be mergeable with above if vectorNormBitSize takes a tag
        for (auto rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            Integer norm = 0;
            for (const auto& pair : *rowIt) {
                norm += (pair.second) * (pair.second);
            }
            if (boundBitSize < minNormBitSize) {
                minNormBitSize = boundBitSize;
            }
            boundBitSize += boundBitSize;
        }
    }

    /**
     * Precise Hadamard bound (bound on determinant) by taking
     * the column-wise euclidean norm.
     *
     * The result is expressed as bit size.
     */
    template <class IMatrix>
    size_t DetailedHadamardColBound(const IMatrix& A, size_t& minNormBitSize)
    {
        minNormBitSize = -1u;
        size_t normBitSize = 0;

        typename IMatrix::ConstColIterator colIt;
        for (colIt = A.colBegin(); colIt != A.colEnd(); ++colIt) {
            size_t colNormBitSize = vectorNormBitSize(colIt->begin(), colIt->end());
            if (colNormBitSize < minNormBitSize) {
                minNormBitSize = colNormBitSize;
            }
            normBitSize += colNormBitSize;
        }

        return normBitSize;
    }

    template <class Field>
    size_t DetailedHadamardColBound(const SparseMatrix<Field>& A, size_t& minNormBitSize)
    {
        return 0;
    }

    /**
     * Precise Hadamard bound (bound on determinant) by taking the minimum
     * of the column-wise and the row-wise euclidean norm.
     *
     * The results are expressed as bit size.
     */
    template <class IMatrix>
    DetailedHadamardBoundData DetailedHadamardBound(const IMatrix& A)
    {
        size_t rowBoundBitSize = 0;
        size_t minRowNormBitSize = 0;
        DetailedHadamardRowBound(rowBoundBitSize, minRowNormBitSize, A);
        size_t boundOnRowNormBitSize = rowBoundBitSize - minRowNormBitSize;

        size_t minColNormBitSize = 0;
        size_t colBoundBitSize = DetailedHadamardColBound(A, minColNormBitSize);
        size_t boundOnColNormBitSize = colBoundBitSize - minColNormBitSize;

        DetailedHadamardBoundData data;
        data.boundBitSize = std::min(rowBoundBitSize, colBoundBitSize);
        data.boundOnMinNormBitSize = std::min(boundOnRowNormBitSize, boundOnColNormBitSize);
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
    size_t HadamardBound(const IMatrix& A)
    {
        return DetailedHadamardBound(A).boundBitSize;
    }

    // ----- Fast Hadamard bound

    /**
     * Returns the bit size of the Hadamard bound.
     * This is a larger estimation but faster to compute.
     */
    template <class IMatrix>
    size_t FastHadamardBound(const IMatrix& A)
    {
        size_t hadamardBitSize;

        // @fixme With SparseMatrix, do we iterate over all zeros?
        // Find a better way.
        Integer max = 0;
        for (auto it = A.Begin(); it != A.End(); ++it) {
            Integer ai = *it;
            if (max < ai)
                max = ai;
            else if (max < -ai)
                max = -ai;
        }

        auto n = std::max(A.rowdim(), A.coldim());
        hadamardBitSize = 1 + n * (Givaro::logtwo(n) / 2 + Givaro::logtwo(max));

        return hadamardBitSize;
    }

    // ----- Rational solve bound

    struct RationalSolveHadamardBoundData {
        size_t numBoundBitSize;
        size_t denBoundBitSize;
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
        size_t bNormBitSize = vectorNormBitSize(b.begin(), b.end());

        data.numBoundBitSize = hadamardBound.boundOnMinNormBitSize + bNormBitSize;
        data.denBoundBitSize = hadamardBound.boundBitSize;
        return data;
    }
}