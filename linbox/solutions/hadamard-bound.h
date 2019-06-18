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
    void vectorNormSquared(Integer& normSquared, const ConstIterator& begin,
                           const ConstIterator& end)
    {
        normSquared = 0;
        for (ConstIterator it = begin; it != end; ++it) {
            // Whatever field element it is,
            // it should be able to store the square without
            // loss of information.
            normSquared += (*it) * (*it);
        }
    }

    // ----- Detailed Hadamard bound

    struct HadamarBoundDetails {
        /**
         * The minimal hadamard bound
         * between the row-wise and the col-wise ones.
         *
         * min { HadamardRow(A), HadamardCol(A) }
         */
        Integer bound;

        /**
         * The minimal hadamard bound
         * divided by the min of the norm vectors
         * between the row-wise and the col-wise ones.
         *
         * min { HadamardRow(A) / min || Ai,* ||,
         *       HadamardCol(A) / min || A*,j ||  }
         */
        Integer boundOverMinNorm;
    };

    /**
     * Precise Hadamard bound (bound on determinant) by taking
     * the row-wise euclidean norm.
     */
    template <class IMatrix>
    void HadamardRowBound(Integer& bound, Integer& minNormSquared, const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        HadamardRowBound(bound, minNormSquared, A, tag);
    }

    template <class IMatrix>
    void HadamardRowBound(Integer& bound, Integer& minNormSquared, const IMatrix& A,
                          const MatrixCategories::RowColMatrixTag& tag)
    {
        bound = 1;
        minNormSquared = -1;

        for (auto rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            Integer rowNormSquared;
            vectorNormSquared(rowNormSquared, rowIt->begin(), rowIt->end());

            if (rowNormSquared == 0) {
                bound = 0;
                minNormSquared = 0;
                return;
            }

            if (minNormSquared < 0 || rowNormSquared < minNormSquared) {
                minNormSquared = rowNormSquared;
            }

            bound *= rowNormSquared;
        }

        bound = Givaro::sqrt(bound);
    }

    template <class IMatrix>
    void HadamardRowBound(Integer& bound, Integer& minNormSquared, const IMatrix& A,
                          const MatrixCategories::RowMatrixTag& tag)
    {
        bound = 1;
        minNormSquared = -1;

        for (auto rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            Integer normSquared = 0;
            for (const auto& pair : *rowIt) {
                normSquared += (pair.second) * (pair.second);
            }

            if (normSquared == 0) {
                bound = 0;
                minNormSquared = 0;
                return;
            }

            if (minNormSquared < 0 || normSquared < minNormSquared) {
                minNormSquared = normSquared;
            }

            bound *= normSquared;
        }

        bound = Givaro::sqrt(bound);
    }

    template <class IMatrix>
    void HadamardRowBound(Integer& bound, Integer& minNormSquared, const IMatrix& A,
                          const MatrixCategories::BlackboxTag& tag)
    {
        DenseMatrix<typename IMatrix::Field> ACopy(A);
        HadamardRowBound(bound, minNormSquared, ACopy);
    }

    /**
     * Precise Hadamard bound (bound on determinant) by taking
     * the column-wise euclidean norm.
     *
     * The result is expressed as bit size.
     */
    template <class IMatrix>
    void HadamardColBound(Integer& bound, Integer& minNormSquared, const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        HadamardColBound(bound, minNormSquared, A, tag);
    }

    template <class IMatrix>
    void HadamardColBound(Integer& bound, Integer& minNormSquared, const IMatrix& A,
                          const MatrixCategories::RowColMatrixTag& tag)
    {
        bound = 1;
        minNormSquared = -1;

        typename IMatrix::ConstColIterator colIt;
        for (colIt = A.colBegin(); colIt != A.colEnd(); ++colIt) {
            Integer colNormSquared;
            vectorNormSquared(colNormSquared, colIt->begin(), colIt->end());

            if (colNormSquared == 0) {
                bound = 0;
                minNormSquared = 0;
                return;
            }

            if (minNormSquared < 0 || colNormSquared < minNormSquared) {
                minNormSquared = colNormSquared;
            }

            bound *= colNormSquared;
        }

        bound = Givaro::sqrt(bound);
    }

    template <class IMatrix>
    void HadamardColBound(Integer& bound, Integer& minNormSquared, const IMatrix& A,
                          const MatrixCategories::RowMatrixTag& tag)
    {
        bound = 1;
        minNormSquared = -1;

        // This vector contains the norm squared for each columns.
        std::vector<Integer> columnsNormsSquared(A.coldim());
        for (auto rowIt = A.rowBegin(); rowIt != A.rowEnd(); ++rowIt) {
            for (const auto& pair : *rowIt) {
                columnsNormsSquared[pair.first] += (pair.second) * (pair.second);
            }
        }

        // All the norms have been computed, we check which one is the smallest
        // and compute the product (aka sum bitsize-wise) of them to make the bound.
        for (const Integer& normSquared : columnsNormsSquared) {
            if (normSquared == 0) {
                bound = 0;
                minNormSquared = 0;
                return;
            }

            if (minNormSquared < 0 || normSquared < minNormSquared) {
                minNormSquared = normSquared;
            }

            bound *= normSquared;
        }

        // Square-root
        bound = Givaro::sqrt(bound);
    }

    template <class IMatrix>
    void HadamardColBound(Integer& bound, Integer& minNormSquared, const IMatrix& A,
                          const MatrixCategories::BlackboxTag& tag)
    {
        DenseMatrix<typename IMatrix::Field> ACopy(A);
        HadamardColBound(bound, minNormSquared, ACopy);
    }

    /**
     * Precise Hadamard bound (bound on determinant) by taking the minimum
     * of the column-wise and the row-wise euclidean norm.
     *
     * The results are expressed as bit size.
     */
    template <class IMatrix>
    HadamarBoundDetails DetailedHadamardBound(const IMatrix& A)
    {
        Integer rowBound;
        Integer rowMinNormSquared;
        HadamardRowBound(rowBound, rowMinNormSquared, A);
        Integer rowBoundOverMinNorm = rowBound / Givaro::sqrt(rowMinNormSquared);
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "rowBound:=" << rowBound << ';' << std::endl;
        std::clog << "rowMinNormSquared:=" << rowMinNormSquared << ';' << std::endl;
        std::clog << "rowBoundOverMinNorm:=" << rowBoundOverMinNorm << ';' << std::endl;
#endif

        Integer colBound;
        Integer colMinNormSquared;
        HadamardColBound(colBound, colMinNormSquared, A);
        Integer colBoundOverMinNorm = colBound / Givaro::sqrt(colMinNormSquared);
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "colBound:=" << colBound << ';' << std::endl;
        std::clog << "colMinNormSquared:=" << colMinNormSquared << ';' << std::endl;
        std::clog << "colBoundOverMinNorm:=" << colBoundOverMinNorm << ';' << std::endl;
#endif

        HadamarBoundDetails data;
        data.bound = (rowBound < colBound) ? rowBound : colBound;
        data.boundOverMinNorm =
            (rowBoundOverMinNorm < colBoundOverMinNorm) ? rowBoundOverMinNorm : colBoundOverMinNorm;
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "bound:=" << data.bound << ';' << std::endl;
        std::clog << "boundOverMinNorm:=" << data.boundOverMinNorm << ';' << std::endl;
#endif

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
        return DetailedHadamardBound(A).bound;
    }

    // ----- Fast Hadamard bound

    template <class IMatrix>
    inline Integer& InfinityNorm(Integer& max, const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        return InfinityNorm(max, A, tag);
    }

    /**
     * Returns the maximal absolute value.
     */
    template <class IMatrix, class MTag>
    inline Integer& InfinityNorm(Integer& max, const IMatrix& A, const MTag& tag)
    {
        DenseMatrix<typename IMatrix::Field> ACopy(A);
        return InfinityNorm(max, ACopy, MatrixCategories::RowColMatrixTag());
    }

    template <class IMatrix>
    inline Integer& InfinityNorm(Integer& max, const IMatrix& A,
                                 const MatrixCategories::RowColMatrixTag& tag)
    {
        max = 0;
        for (auto it = A.Begin(); it != A.End(); ++it) {
            const Integer& ai = *it;
            if (max < ai)
                max = ai;
            else if (max < -ai)
                max = -ai;
        }

        return max;
    }

    /**
     * Returns the bit size of the Hadamard bound.
     * This is a larger estimation but faster to compute.
     */
    template <class IMatrix>
    inline double FastHadamardLogBound(const IMatrix& A, const Integer& infinityNorm)
    {
        if (infinityNorm == 0) {
            return 0.0;
        }

        uint64_t n = std::max(A.rowdim(), A.coldim());
        double bound =
            static_cast<double>(n) * (Givaro::logtwo(n) / 2.0 + Givaro::logtwo(infinityNorm));
        return bound;
    }

    template <class IMatrix>
    inline double FastHadamardLogBound(const IMatrix& A,
                                       const MatrixCategories::RowColMatrixTag& tag)
    {
        Integer infinityNorm;
        InfinityNorm(infinityNorm, A, tag);
        return FastHadamardLogBound(A, infinityNorm);
    }

    template <class IMatrix>
    inline double FastHadamardLogBound(const IMatrix& A, const MatrixCategories::BlackboxTag& tag)
    {
        DenseMatrix<typename IMatrix::Field> ACopy(A);
        return FastHadamardLogBound(ACopy);
    }

    template <class IMatrix>
    inline double FastHadamardLogBound(const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        return FastHadamardLogBound(A, tag);
    }

    /**
     * Bound on the coefficients of the characteristic polynomial
     * @bib "Efficient Computation of the Characteristic Polynomial". Dumas Pernet Wan ISSAC'05.
     *
     */
    template <class IMatrix>
    inline double FastCharPolyDumasPernetWanBound(const IMatrix& A, const Integer& infinityNorm)
    {
        // .105815875 = 0.21163275 / 2
        return FastHadamardLogBound(A, infinityNorm) + A.coldim() * .105815875;
    }

    /**
     * A.J. Goldstein et R.L. Graham.
     * A Hadamard-type bound on the coefficients of
     * a determinant of polynomials.
     * SIAM Review, volume 15, 1973, pages 657-658.
     *
     */
    template <class IMatrix>
    inline double FastCharPolyGoldsteinGrahamBound(const IMatrix& A, const Integer& infinityNorm)
    {
        Integer ggb(infinityNorm);
        ggb *= static_cast<uint64_t>(A.coldim());
        ggb += 2;
        ggb *= infinityNorm;
        ++ggb;
        return Givaro::logtwo(ggb) * A.coldim() / 2.0;
    }

    template <class IMatrix>
    inline double FastCharPolyHadamardBound(const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        Integer infinityNorm;
        InfinityNorm(infinityNorm, A, tag);
        const double DPWbound = FastCharPolyDumasPernetWanBound(A, infinityNorm);
        const double GGbound = FastCharPolyGoldsteinGrahamBound(A, infinityNorm);
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "DPWbound: " << DPWbound << std::endl;
        std::clog << "GGbound : " << GGbound << std::endl;
#endif
        return std::min(DPWbound, GGbound);
    }

    // ----- Rational solve bound

    struct RationalSolveHadamardBoundData {
        Integer numBound;        // N
        Integer denBound;        // D
        double solutionLogBound; // log2(2 * N * D)
    };

    /**
     * Bound on the rational solution of a linear system Ax = b.
     *
     * Return bounds on the bit sizes of both denominator and numerator of the solution x.
     *
     * @note Matrix and Vector should be over Integer.
     */
    template <class Matrix, class Vector>
    typename std::enable_if<std::is_same<typename FieldTraits<typename Matrix::Field>::categoryTag,
                                         RingCategories::IntegerTag>::value,
                            RationalSolveHadamardBoundData>::type
    RationalSolveHadamardBound(const Matrix& A, const Vector& b)
    {
        RationalSolveHadamardBoundData data;

        auto hadamardBound = DetailedHadamardBound(A);
        Integer bNormSquared;
        vectorNormSquared(bNormSquared, b.begin(), b.end());

        data.denBound = hadamardBound.bound;
        data.numBound = hadamardBound.boundOverMinNorm * Givaro::sqrt(bNormSquared);
        if (data.denBound == 0 || data.numBound == 0) {
            data.solutionLogBound = 0.0;
        }
        else {
            data.solutionLogBound = 1.0 + Givaro::logtwo(data.numBound)
                                    + Givaro::logtwo(data.denBound); // log2(2 * N * D)
        }

#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "numBound:=" << data.numBound << ';' << std::endl;
        std::clog << "denBound:=" << data.denBound << ';' << std::endl;
#endif
        return data;
    }

    /// @fixme Needed to solve-cra.h, but can't be used yet.
    template <class Matrix, class Vector>
    typename std::enable_if<std::is_same<typename FieldTraits<typename Matrix::Field>::categoryTag,
                                         RingCategories::RationalTag>::value,
                            RationalSolveHadamardBoundData>::type
    RationalSolveHadamardBound(const Matrix& A, const Vector& b)
    {
        throw NotImplementedYet("Hadamard bound on Rational matrices is not implemented yet.");
    }
}
