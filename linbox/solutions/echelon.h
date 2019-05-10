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
 * Copyright (C) 2019 the LinBox group 
 * Written by Clément Pernet
 */

#pragma once

#include <iostream> // @note Needed for givaro/ring-interface to compile

#include <linbox/field/field-traits.h>
#include <linbox/solutions/methods.h>

namespace LinBox {
    //
    // row echelon form
    //

    /**
     * \brief Compute the row echelon form of a matrix, not reduced.
     *
     * Returns the row echelon form. The pivots are nonzero but not necessarily ones.
     * The form is not reduced, which means that coefficients above each pivot are not necessarily zeros.
     *
     * @param [out] E the row echelon form 
     * @param [in]  A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t rowEchelon (Matrix & E, const Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("rowEchelon<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief rowEchelon dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t rowEchelon (Matrix & E, const Matrix& A, const EchelonMethod& m)
    {
        return rowEchelon (E, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief rowEchelon dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t rowEchelon (Matrix & E, const Matrix& A)
    {
        return rowEchelon (E, A, Method::Auto());
    }


    /**
     * \brief Compute the row echelon form of a matrix, not reduced, and the related transformation matrix.
     *
     * Returns the row echelon form E and a transformation matrix T such that T . A = E. 
     * The pivots are nonzero but not necessarily ones.
     * The form is not reduced, which means that coefficients above each pivot are not necessarily zeros.
     *
     * @param [out] E the row echelon form 
     * @param [out] T the transformation matrix 
     * @param [in]  A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t rowEchelon (Matrix & E, Matrix& T, const Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("rowEchelon<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief rowEchelon dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t rowEchelon (Matrix & E, Matrix& T, const Matrix& A, const EchelonMethod& m)
    {
        return rowEchelon (E, T, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief rowEchelon dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t rowEchelon (Matrix & E, Matrix& T, const Matrix& A)
    {
        return rowEchelon (E, T, A, Method::Auto());
    }

    //
    // row echelonize
    //

    /**
     * \brief Replace the input matrix by its row echelon form, not reduced.
     *
     * The pivots are nonzero but not necessarily ones.
     * The form is not reduced, which means that coefficients above each pivot are not necessarily zeros.
     *
     * @param [in/out] A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t rowEchelonize (Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("rowEchelonize<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief rowEchelonize dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t rowEchelonize (Matrix& A, const EchelonMethod& m)
    {
        return rowEchelonize (A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief rowEchelonize dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t rowEchelonize (Matrix & A)
    {
        return rowEchelonize (A, Method::Auto());
    }


    /**
     * \brief Compute the row echelon form of a matrix, not reduced, and the related transformation matrix.
     *
     * Returns the row echelon form E and a transformation matrix T such that T . A = E. 
     * The input matrix is replaced by the row echelon form.
     * The pivots are nonzero but not necessarily ones.
     * The form is not reduced, which means that coefficients above each pivot are not necessarily zeros.
     *
     * @param [in/out] A matrix
     * @param [out] T the transformation matrix 
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t rowEchelonize (Matrix & A, Matrix& T, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("rowEchelonize<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief rowEchelonize dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t rowEchelonize (Matrix & A, Matrix& T, const EchelonMethod& m)
    {
        return rowEchelonize (A, T, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief rowEchelonize dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t rowEchelonize (Matrix & A, Matrix& T)
    {
        return rowEchelonize (A, T, Method::Auto());
    }

    //
    // reduced row echelon form
    //

    /**
     * \brief Compute the reduced row echelon form of a matrix.
     *
     * Returns the reduced row echelon form. The pivots are ones.
     * The form is reduced, which means that coefficients above each pivot are zeros.
     *
     * @param [out] E the reduced row echelon form 
     * @param [in]  A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t reducedRowEchelon (Matrix & E, const Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("reducedRowEchelon<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief reducedRowEchelon dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t reducedRowEchelon (Matrix & E, const Matrix& A, const EchelonMethod& m)
    {
        return reducedRowEchelon (E, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief reducedRowEchelon dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t reducedRowEchelon (Matrix & E, const Matrix& A)
    {
        return reducedRowEchelon (E, A, Method::Auto());
    }


    /**
     * \brief Compute the reduced row echelon form of a matrix, and the related transformation matrix.
     *
     * Returns the reduced row echelon form E and a transformation matrix T such that T . A = E. 
     * The pivots are ones. The form is reduced, which means that coefficients above each pivot are zeros.
     *
     * @param [out] E the reduced row echelon form 
     * @param [out] T the transformation matrix 
     * @param [in]  A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t reducedRowEchelon (Matrix & E, Matrix& T, const Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("reducedRowEchelon<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief reducedRowEchelon dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t reducedRowEchelon (Matrix & E, Matrix& T, const Matrix& A, const EchelonMethod& m)
    {
        return reducedRowEchelon (E, T, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief reducedRowEchelon dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t reducedRowEchelon (Matrix & E, Matrix& T, const Matrix& A)
    {
        return reducedRowEchelon (E, T, A, Method::Auto());
    }

    //
    // reduced row echelonize
    //

    /**
     * \brief Replace the input matrix by its reduced row echelon form.
     *
     * The pivots are ones. The form is reduced, which means that coefficients above each pivot are zeros.
     *
     * @param [in/out] A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t reducedRowEchelonize (Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("reducedRowEchelonize<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief reducedRowEchelonize dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t reducedRowEchelonize (Matrix& A, const EchelonMethod& m)
    {
        return reducedRowEchelonize (A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief reducedRowEchelonize dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t reducedRowEchelonize (Matrix & A)
    {
        return reducedRowEchelonize (A, Method::Auto());
    }


    /**
     * \brief Compute the reduced row echelon form of a matrix, and the related transformation matrix.
     *
     * Returns the reduced row echelon form E and a transformation matrix T such that T . A = E. 
     * The input matrix is replaced by the reduced row echelon form.
     * The pivots are ones. The form is reduced, which means that coefficients above each pivot are zeros.
     *
     * @param [in/out] A matrix
     * @param [out] T the transformation matrix 
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t reducedRowEchelonize (Matrix & A, Matrix& T, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("reducedRowEchelonize<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief reducedRowEchelonize dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t reducedRowEchelonize (Matrix & A, Matrix& T, const EchelonMethod& m)
    {
        return reducedRowEchelonize (A, T, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief reducedRowEchelonize dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t reducedRowEchelonize (Matrix & A, Matrix& T)
    {
        return reducedRowEchelonize (A, T, Method::Auto());
    }

    //
    // column echelon form
    //

    /**
     * \brief Compute the column echelon form of a matrix, not reduced.
     *
     * Returns the column echelon form. The pivots are nonzero but not necessarily ones.
     * The form is not reduced, which means that coefficients above each pivot are not necessarily zeros.
     *
     * @param [out] E the column echelon form 
     * @param [in]  A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t colEchelon (Matrix & E, const Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("colEchelon<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief colEchelon dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t colEchelon (Matrix & E, const Matrix& A, const EchelonMethod& m)
    {
        return colEchelon (E, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief colEchelon dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t colEchelon (Matrix & E, const Matrix& A)
    {
        return colEchelon (E, A, Method::Auto());
    }


    /**
     * \brief Compute the column echelon form of a matrix, not reduced, and the related transformation matrix.
     *
     * Returns the column echelon form E and a transformation matrix T such that T . A = E. 
     * The pivots are nonzero but not necessarily ones.
     * The form is not reduced, which means that coefficients above each pivot are not necessarily zeros.
     *
     * @param [out] E the column echelon form 
     * @param [out] T the transformation matrix 
     * @param [in]  A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t colEchelon (Matrix & E, Matrix& T, const Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("colEchelon<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief colEchelon dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t colEchelon (Matrix & E, Matrix& T, const Matrix& A, const EchelonMethod& m)
    {
        return colEchelon (E, T, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief colEchelon dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t colEchelon (Matrix & E, Matrix& T, const Matrix& A)
    {
        return colEchelon (E, T, A, Method::Auto());
    }

    //
    // column echelonize
    //

    /**
     * \brief Replace the input matrix by its column echelon form, not reduced.
     *
     * The pivots are nonzero but not necessarily ones.
     * The form is not reduced, which means that coefficients above each pivot are not necessarily zeros.
     *
     * @param [in/out] A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t colEchelonize (Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("colEchelonize<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief colEchelonize dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t colEchelonize (Matrix& A, const EchelonMethod& m)
    {
        return colEchelonize (A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief colEchelonize dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t colEchelonize (Matrix & A)
    {
        return colEchelonize (A, Method::Auto());
    }


    /**
     * \brief Compute the column echelon form of a matrix, not reduced, and the related transformation matrix.
     *
     * Returns the column echelon form E and a transformation matrix T such that T . A = E. 
     * The input matrix is replaced by the column echelon form.
     * The pivots are nonzero but not necessarily ones.
     * The form is not reduced, which means that coefficients above each pivot are not necessarily zeros.
     *
     * @param [in/out] A matrix
     * @param [out] T the transformation matrix 
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t colEchelonize (Matrix & A, Matrix& T, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("colEchelonize<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief colEchelonize dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t colEchelonize (Matrix & A, Matrix& T, const EchelonMethod& m)
    {
        return colEchelonize (A, T, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief colEchelonize dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t colEchelonize (Matrix & A, Matrix& T)
    {
        return colEchelonize (A, T, Method::Auto());
    }

    //
    // reduced column echelon form
    //

    /**
     * \brief Compute the reduced column echelon form of a matrix.
     *
     * Returns the reduced column echelon form. The pivots are ones.
     * The form is reduced, which means that coefficients above each pivot are zeros.
     *
     * @param [out] E the reduced column echelon form 
     * @param [in]  A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t reducedColEchelon (Matrix & E, const Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("reducedColEchelon<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief reducedColEchelon dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t reducedColEchelon (Matrix & E, const Matrix& A, const EchelonMethod& m)
    {
        return reducedColEchelon (E, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief reducedColEchelon dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t reducedColEchelon (Matrix & E, const Matrix& A)
    {
        return reducedColEchelon (E, A, Method::Auto());
    }


    /**
     * \brief Compute the reduced column echelon form of a matrix, and the related transformation matrix.
     *
     * Returns the reduced column echelon form E and a transformation matrix T such that T . A = E. 
     * The pivots are ones. The form is reduced, which means that coefficients above each pivot are zeros.
     *
     * @param [out] E the reduced column echelon form 
     * @param [out] T the transformation matrix 
     * @param [in]  A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t reducedColEchelon (Matrix & E, Matrix& T, const Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("reducedColEchelon<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief reducedColEchelon dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t reducedColEchelon (Matrix & E, Matrix& T, const Matrix& A, const EchelonMethod& m)
    {
        return reducedColEchelon (E, T, A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief reducedColEchelon dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t reducedColEchelon (Matrix & E, Matrix& T, const Matrix& A)
    {
        return reducedColEchelon (E, T, A, Method::Auto());
    }

    //
    // reduced column echelonize
    //

    /**
     * \brief Replace the input matrix by its reduced column echelon form.
     *
     * The pivots are ones. The form is reduced, which means that coefficients above each pivot are zeros.
     *
     * @param [in/out] A matrix
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t reducedColEchelonize (Matrix& A, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("reducedColEchelonize<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief reducedColEchelonize dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t reducedColEchelonize (Matrix& A, const EchelonMethod& m)
    {
        return reducedColEchelonize (A, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief reducedColEchelonize dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t reducedColEchelonize (Matrix & A)
    {
        return reducedColEchelonize (A, Method::Auto());
    }


    /**
     * \brief Compute the reduced column echelon form of a matrix, and the related transformation matrix.
     *
     * Returns the reduced column echelon form E and a transformation matrix T such that T . A = E. 
     * The input matrix is replaced by the reduced column echelon form.
     * The pivots are ones. The form is reduced, which means that coefficients above each pivot are zeros.
     *
     * @param [in/out] A matrix
     * @param [out] T the transformation matrix 
     * @param [in]  tag domain of computation
     * @param [in]  m method to use (\see solutions/method.h)
     * @return the rank of the matrix 
     */
    template <class Matrix, class CategoryTag, class EchelonMethod>
    inline size_t reducedColEchelonize (Matrix & A, Matrix& T, const CategoryTag& tag, const EchelonMethod& m)
    {
        throw LinBoxError("reducedColEchelonize<" + CategoryTag::name() + ", " + EchelonMethod::name() + "> does not exists.");
    }

    /**
     * \brief reducedColEchelonize dispatcher for automated category tag.
     */
    template <class Matrix, class EchelonMethod>
    inline size_t reducedColEchelonize (Matrix & A, Matrix& T, const EchelonMethod& m)
    {
        return reducedColEchelonize (A, T, typename FieldTraits<typename Matrix::Field>::categoryTag(), m);
    }

    /**
     * \brief reducedColEchelonize dispatcher for automated method.
     */
    template <class Matrix>
    inline size_t reducedColEchelonize (Matrix & A, Matrix& T)
    {
        return reducedColEchelonize (A, T, Method::Auto());
    }
}

//
// echelon
//

#include "./echelon/echelon-auto.h"

// Elimination
#include "./echelon/echelon-dense-elimination.h"
// #include "./echelon/echelon-elimination.h"
// #include "./echelon/echelon-sparse-elimination.h"

/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

