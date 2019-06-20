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
 *
 * Written by Clément Pernet
 */

#pragma once

#include <linbox/matrix/dense-matrix.h>
#include <linbox/matrix/sparse-matrix.h>
#include <linbox/solutions/methods.h>
namespace LinBox {
    //
    // row echelon
    //

    /**
     * \brief rowEchelon specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t rowEchelon (Matrix & E, const Matrix& A, const CategoryTag& tag, const Method::Auto& m)
    {
        return rowEchelon (E, A, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief rowEchelon specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t rowEchelon (DenseMatrix<Field>& E, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return rowEchelon (E, A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief rowEchelon specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t rowEchelon (Matrix & E, Matrix & T, const Matrix& A, const CategoryTag& tag, const Method::Auto& m)
    {
        return rowEchelon (E, A, T, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief rowEchelon specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t rowEchelon (DenseMatrix<Field>& E, DenseMatrix<Field>& T, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return rowEchelon (E, A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    //
    // row echelonize
    //

    /**
     * \brief rowEchelonize specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t rowEchelonize (Matrix & A, const CategoryTag& tag, const Method::Auto& m)
    {
        return rowEchelonize (A, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief rowEchelonize specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t rowEchelonize (DenseMatrix<Field>& A,
                                 const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return rowEchelonize (A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief rowEchelonize specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t rowEchelonize (Matrix & A, Matrix & T, const CategoryTag& tag, const Method::Auto& m)
    {
        return rowEchelonize (A, T, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief rowEchelonize specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t rowEchelonize (DenseMatrix<Field>& A, DenseMatrix<Field>& T,
                                 const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return rowEchelonize (A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    //
    // reduced row echelon
    //

    /**
     * \brief reducedRowEchelon specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t reducedRowEchelon (Matrix & E, const Matrix& A, const CategoryTag& tag, const Method::Auto& m)
    {
        return reducedRowEchelon (E, A, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief reducedRowEchelon specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t reducedRowEchelon (DenseMatrix<Field>& E, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return reducedRowEchelon (E, A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief reducedRowEchelon specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t reducedRowEchelon (Matrix & E, Matrix & T, const Matrix& A, const CategoryTag& tag, const Method::Auto& m)
    {
        return reducedRowEchelon (E, A, T, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief reducedRowEchelon specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t reducedRowEchelon (DenseMatrix<Field>& E, DenseMatrix<Field>& T, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return reducedRowEchelon (E, A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    //
    // row echelonize
    //

    /**
     * \brief reducedRowEchelonize specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t reducedRowEchelonize (Matrix & A, const CategoryTag& tag, const Method::Auto& m)
    {
        return reducedRowEchelonize (A, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief reducedRowEchelonize specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t reducedRowEchelonize (DenseMatrix<Field>& A,
                                 const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return reducedRowEchelonize (A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief reducedRowEchelonize specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t reducedRowEchelonize (Matrix & A, Matrix & T, const CategoryTag& tag, const Method::Auto& m)
    {
        return reducedRowEchelonize (A, T, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief reducedRowEchelonize specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t reducedRowEchelonize (DenseMatrix<Field>& A, DenseMatrix<Field>& T,
                                 const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return reducedRowEchelonize (A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }
    //
    // column echelon
    //

    /**
     * \brief colEchelon specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t colEchelon (Matrix & E, const Matrix& A, const CategoryTag& tag, const Method::Auto& m)
    {
        return colEchelon (E, A, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief colEchelon specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t colEchelon (DenseMatrix<Field>& E, const DenseMatrix<Field>& A,
                       const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return colEchelon (E, A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief colEchelon specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t colEchelon (Matrix & E, Matrix & T, const Matrix& A, const CategoryTag& tag, const Method::Auto& m)
    {
        return colEchelon (E, A, T, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief colEchelon specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t colEchelon (DenseMatrix<Field>& E, DenseMatrix<Field>& T, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return colEchelon (E, A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    //
    // col echelonize
    //

    /**
     * \brief colEchelonize specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t colEchelonize (Matrix & A, const CategoryTag& tag, const Method::Auto& m)
    {
        return colEchelonize (A, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief colEchelonize specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t colEchelonize (DenseMatrix<Field>& A,
                          const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return colEchelonize (A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief colEchelonize specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t colEchelonize (Matrix & A, Matrix & T, const CategoryTag& tag, const Method::Auto& m)
    {
        return colEchelonize (A, T, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief colEchelonize specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t colEchelonize (DenseMatrix<Field>& A, DenseMatrix<Field>& T,
                          const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return colEchelonize (A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    //
    // reduced col echelon
    //

    /**
     * \brief reducedColEchelon specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t reducedColEchelon (Matrix & E, const Matrix& A, const CategoryTag& tag, const Method::Auto& m)
    {
        return reducedColEchelon (E, A, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief reducedColEchelon specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t reducedColEchelon (DenseMatrix<Field>& E, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return reducedColEchelon (E, A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief reducedColEchelon specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t reducedColEchelon (Matrix & E, Matrix & T, const Matrix& A, const CategoryTag& tag, const Method::Auto& m)
    {
        return reducedColEchelon (E, A, T, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief reducedColEchelon specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t reducedColEchelon (DenseMatrix<Field>& E, DenseMatrix<Field>& T, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return reducedColEchelon (E, A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    //
    // col echelonize
    //

    /**
     * \brief reducedColEchelonize specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t reducedColEchelonize (Matrix & A, const CategoryTag& tag, const Method::Auto& m)
    {
        return reducedColEchelonize (A, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief reducedColEchelonize specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t reducedColEchelonize (DenseMatrix<Field>& A,
                                 const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return reducedColEchelonize (A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

    /**
     * \brief reducedColEchelonize specialisation for Auto.
     */
    template <class Matrix, class CategoryTag>
    size_t reducedColEchelonize (Matrix & A, Matrix & T, const CategoryTag& tag, const Method::Auto& m)
    {
        return reducedColEchelonize (A, T, tag, reinterpret_cast<const Method::Elimination&>(m));
    }

    /**
     * \brief reducedColEchelonize specialisation for Auto with DenseMatrix and ModularTag.
     */
    template <class Field>
    size_t reducedColEchelonize (DenseMatrix<Field>& A, DenseMatrix<Field>& T,
                                 const RingCategories::ModularTag& tag, const Method::Auto& m)
    {
        return reducedColEchelonize (A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

