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

    // Returns false if the vector is null, true otherwise
    template <class ConstIterator>
    bool vectorLogNorm(double& logNorm, const ConstIterator& begin, const ConstIterator& end)
    {
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "BEG vectorLogNorm\n" ;
#endif
        Integer normSquared = 0;
        for (ConstIterator it = begin; it != end; ++it) {
            // normSquared += Integer(*it)*Integer(*it);
            Integer rit(*it);
            Integer::axpyin(normSquared, rit, rit);
        }

        if (normSquared == 0) {
            logNorm = 0.0;
            return false; // Vector is zero
        }

        logNorm = Givaro::logtwo(normSquared) / 2.0;

#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "normSquared:=" << normSquared << ';' << std::endl;
        std::clog << "vectorLogNorm:=" << logNorm << ';' << std::endl;
        std::clog << "END vectorLogNorm\n" ;
#endif
        return true;
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
            double rowLogNorm;
            if (vectorLogNorm(rowLogNorm, rowIt->begin(), rowIt->end())) {
                if (rowLogNorm < minLogNorm) {
                    minLogNorm = rowLogNorm;
                }
            }
            else {
                logBound = 0.0;
                minLogNorm = 0.0;
                return;
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

    template <class IMatrix>
    void HadamardRowLogBound(double& logBound, double& minLogNorm, const IMatrix& A, const MatrixCategories::BlackboxTag& tag)
    {
        DenseMatrix<typename IMatrix::Field> ACopy(A);
        HadamardRowLogBound(logBound, minLogNorm, ACopy);
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
            double colLogNorm;
            if (vectorLogNorm(colLogNorm, colIt->begin(), colIt->end())) {
                if (colLogNorm < minLogNorm) {
                    minLogNorm = colLogNorm;
                }
            }
            else {
                logBound = 0.0;
                minLogNorm = 0.0;
                return;
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

    template <class IMatrix>
    void HadamardColLogBound(double& logBound, double& minLogNorm, const IMatrix& A, const MatrixCategories::BlackboxTag& tag)
    {
        DenseMatrix<typename IMatrix::Field> ACopy(A);
        HadamardColLogBound(logBound, minLogNorm, ACopy);
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
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "BEG DetailedHadamardBound\n" ;
#endif
        double rowLogBound = 0.0;
        double rowMinLogNorm = 0.0;
        HadamardRowLogBound(rowLogBound, rowMinLogNorm, A);
        double rowLogBoundOverMinNorm = rowLogBound - rowMinLogNorm;
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "rowLogBound:=" << rowLogBound << ';' << std::endl;
        std::clog << "rowMinLogNorm:=" << rowMinLogNorm << ';' << std::endl;
        std::clog << "rowLogBoundOverMinNorm:=" << rowLogBoundOverMinNorm << ';' << std::endl;
#endif

        double colLogBound = 0.0;
        double colMinLogNorm = 0.0;
        HadamardColLogBound(colLogBound, colMinLogNorm, A);
        double colLogBoundOverMinNorm = colLogBound - colMinLogNorm;
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "colLogBound:=" << colLogBound << ';' << std::endl;
        std::clog << "colMinLogNorm:=" << colMinLogNorm << ';' << std::endl;
        std::clog << "colLogBoundOverMinNorm:=" << colLogBoundOverMinNorm << ';' << std::endl;
#endif

        HadamardLogBoundDetails data;
        data.logBound = std::min(rowLogBound, colLogBound);
        data.logBoundOverMinNorm = std::min(rowLogBoundOverMinNorm, colLogBoundOverMinNorm);
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "logBound:=" << data.logBound << ';' << std::endl;
        std::clog << "logBoundOverMinNorm:=" << data.logBoundOverMinNorm << ';' << std::endl;
        std::clog << "END DetailedHadamardBound\n" ;
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
        return DetailedHadamardBound(A).logBound;
    }

    // ----- Fast Hadamard bound



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
    inline Integer& InfinityNorm(Integer& max, const IMatrix& A, const MatrixCategories::RowColMatrixTag& tag)
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
    inline double FastHadamardBound(const IMatrix& A, const Integer& infnorm)
    {
        if (infnorm == 0) {
            return 0.0;
        }

        uint64_t n = std::max(A.rowdim(), A.coldim());
        double logBound = static_cast<double>(n) * (Givaro::logtwo(n) / 2.0 + Givaro::logtwo(infnorm));
        return logBound;
    }

    template <class IMatrix>
    inline double FastHadamardBound(const IMatrix& A, const MatrixCategories::RowColMatrixTag& tag)
    {
        Integer infnorm;
        InfinityNorm(infnorm, A, tag);
        return FastHadamardBound(A, infnorm);
    }

    template <class IMatrix>
    inline double FastHadamardBound(const IMatrix& A, const MatrixCategories::BlackboxTag& tag)
    {
        DenseMatrix<typename IMatrix::Field> ACopy(A);
        return FastHadamardBound(ACopy);
    }

    template <class IMatrix>
    inline double FastHadamardBound(const IMatrix& A)
    {
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        return FastHadamardBound(A, tag);
    }

        /**
         * Bound on the coefficients of the characteristic polynomial
         * @bib "Efficient Computation of the Characteristic Polynomial". Dumas Pernet Wan ISSAC'05.
         *
         */
    template <class IMatrix>
    inline double FastCharPolyDumasPernetWanBound(const IMatrix& A, const Integer& infnorm)
    {
		// .105815875 = 0.21163275 / 2
        return FastHadamardBound(A, infnorm) + A.coldim()*.105815875;
    }

        /**
         * A.J. Goldstein et R.L. Graham.
         * A Hadamard-type bound on the coefficients of
         * a determinant of polynomials.
         * SIAM Review, volume 15, 1973, pages 657-658.
         *
         */
    template <class IMatrix>
    inline double FastCharPolyGoldsteinGrahamBound(const IMatrix& A, const Integer& infnorm)
    {
        Integer ggb(infnorm);
        ggb *= static_cast<uint64_t>(A.coldim());
        ggb += 2;
        ggb *= infnorm;
        ++ggb;
        return Givaro::logtwo(ggb)*A.coldim()/2.0;
    }

    template <class IMatrix>
    inline double FastCharPolyHadamardBound(const IMatrix& A)
    {
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "BEG FastCharPolyHadamardBound\n" ;
#endif
        typename MatrixTraits<IMatrix>::MatrixCategory tag;
        Integer infnorm;
        InfinityNorm(infnorm, A, tag);
        const double DPWbound = FastCharPolyDumasPernetWanBound(A, infnorm);
        const double GGbound = FastCharPolyGoldsteinGrahamBound(A, infnorm);
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "DPWbound: " << DPWbound << std::endl;
        std::clog << "GGbound : " << GGbound << std::endl;
        std::clog << "END FastCharPolyHadamardBound\n" ;
#endif
        return std::min(DPWbound,GGbound);
    }


    // ----- Rational solve bound

    struct RationalSolveHadamardBoundData {
        double numLogBound;      // log2(N)
        double denLogBound;      // log2(D)
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
    typename std::enable_if<std::is_same<typename FieldTraits<typename Matrix::Field>::categoryTag, RingCategories::IntegerTag>::value,
                            RationalSolveHadamardBoundData>::type
    RationalSolveHadamardBound(const Matrix& A, const Vector& b)
    {
#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "BEG RationalSolveHadamardBound\n" ;
#endif
        RationalSolveHadamardBoundData data;

        auto hadamardBound = DetailedHadamardBound(A);
        double bLogNorm;
        vectorLogNorm(bLogNorm, b.begin(), b.end());

        data.numLogBound = hadamardBound.logBoundOverMinNorm + bLogNorm + 1.0;
        data.denLogBound = hadamardBound.logBound;
        data.solutionLogBound = data.numLogBound + data.denLogBound + 1.0;

#ifdef DEBUG_HADAMARD_BOUND
        std::clog << "numLogBound:=" << data.numLogBound << ';' << std::endl;
        std::clog << "denLogBound:=" << data.denLogBound << ';' << std::endl;
        std::clog << "END RationalSolveHadamardBound\n" ;
#endif
        return data;
    }

    /// @fixme Needed to solve-cra.h, but can't be used yet.
    template <class Matrix, class Vector>
    typename std::enable_if<std::is_same<typename FieldTraits<typename Matrix::Field>::categoryTag, RingCategories::RationalTag>::value,
                            RationalSolveHadamardBoundData>::type
    RationalSolveHadamardBound(const Matrix& A, const Vector& b)
    {
        throw NotImplementedYet("Hadamard bound on Rational matrices is not implemented yet.");
    }
}
