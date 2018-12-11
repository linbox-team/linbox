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

namespace LinBox {

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
    template <class BlackBox>
    size_t DetailedHadamardRowBound(const BlackBox& A, size_t& minNormBitSize)
    {
        minNormBitSize = -1u;
        size_t normBitSize = 0;

        typename BlackBox::ConstRowIterator rowIt;
        for (rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            Integer norm = 0;
            typename BlackBox::ConstRow::const_iterator col;
            for (col = rowIt->begin(); col != rowIt->end(); ++col) {
                norm += static_cast<Integer>((*col)) * (*col);
            }

            size_t rowNormBitSize = std::ceil(Givaro::logtwo(norm) / 2.0);
            if (rowNormBitSize < minNormBitSize) {
                minNormBitSize = rowNormBitSize;
            }
            normBitSize += rowNormBitSize;
        }

        return normBitSize;
    }

    /**
     * Precise Hadamard bound (bound on determinant) by taking
     * the column-wise euclidean norm.
     *
     * The result is expressed as bit size.
     */
    template <class BlackBox>
    size_t DetailedHadamardColBound(const BlackBox& A, size_t& minNormBitSize)
    {
        minNormBitSize = -1u;
        size_t normBitSize = 0;

        typename BlackBox::ConstColIterator colIt;
        for (colIt = A.colBegin(); colIt != A.colEnd(); ++colIt) {
            Integer norm = 0;
            typename BlackBox::ConstCol::const_iterator row;
            for (row = colIt->begin(); row != colIt->end(); ++row) {
                norm += static_cast<Integer>((*row)) * (*row);
            }

            size_t colNormBitSize = std::ceil(Givaro::logtwo(norm) / 2.0);
            if (colNormBitSize < minNormBitSize) {
                minNormBitSize = colNormBitSize;
            }
            normBitSize += colNormBitSize;
        }

        return normBitSize;
    }

    /**
     * Precise Hadamard bound (bound on determinant) by taking the minimum
     * of the column-wise and the row-wise euclidean norm.
     *
     * The results are expressed as bit size.
     */
    template <class BlackBox>
    DetailedHadamardBoundData DetailedHadamardBound(const BlackBox& A)
    {
        size_t minRowNormBitSize = 0;
        size_t rowBoundBitSize = DetailedHadamardRowBound(A, minRowNormBitSize);
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
    template <class BlackBox>
    size_t HadamardBound(const BlackBox& A)
    {
        return DetailedHadamardBound(A).boundBitSize;
    }

    // ----- Fast Hadamard bound

    /**
     * Returns the bit size of the Hadamard bound.
     * This is a larger estimation but faster to compute.
     */
    template <class BlackBox>
    size_t FastHadamardBound(const BlackBox& A)
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
    template <class BlackBox, class Vector>
    RationalSolveHadamardBoundData RationalSolveHadamardBound(const BlackBox& A, const Vector& b)
    {
        RationalSolveHadamardBoundData data;

        auto hadamardBound = DetailedHadamardBound(A);

        // Compute the vector norm
        Integer bNorm = 0;
        for (auto bIt = b.begin(); bIt != b.end(); ++bIt) {
            bNorm += static_cast<Integer>((*bIt)) * (*bIt);
        }
        size_t bNormBitSize = std::ceil(Givaro::logtwo(bNorm) / 2.0);

        data.numBoundBitSize = hadamardBound.boundOnMinNormBitSize + bNormBitSize;
        data.denBoundBitSize = hadamardBound.boundBitSize;
        return data;
    }
}
