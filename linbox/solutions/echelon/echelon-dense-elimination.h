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
         * \brief rowEchelon specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t rowEchelon (DenseMatrix<Field>& E, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        linbox_check((A.coldim() == E.coldim()) && (A.rowdim() == E.rowdim()));

        size_t m = A.rowdim();
        size_t n = A.coldim();
        size_t* P = new size_t[m];
        size_t* Q = new size_t[n];

        size_t R = FFPACK::RowEchelonForm (F, m, n, A.getPointer(), A.getStride(), P, Q, false);

        FFPACK::getEchelonForm (F, FFLAS::FflasUpper, FFLAS::FflasUnit, m, n, R, Q,
                                A.getPointer(), A.getStride(), E.getPointer(), E.getStride());
        
        delete[] P;
        delete[] Q;
        return R;
    }


        /**
         * \brief rowEchelon with transformation specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t rowEchelon (DenseMatrix<Field>& E, DenseMatrix<Field>& T, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        linbox_check((A.coldim() == E.coldim()) && (A.rowdim() == E.rowdim()) &&
                     (A.rowdim() == T.rowdim()) && (T.rowdim() == T.coldim()) );

        size_t m = A.rowdim();
        size_t n = A.coldim();
        size_t* P = new size_t[m];
        size_t* Q = new size_t[n];
        Field& F = A.getField();

        size_t R = FFPACK::RowEchelonForm (F, m, n, A.getPointer(), A.getStride(), P, Q, false);

        FFPACK::getEchelonForm (F, FFLAS::FflasUpper, FFLAS::FflasUnit, m, n, R, Q,
                                A.getPointer(), A.getStride(), E.getPointer(), E.getStride());
        FFPACK::getEchelonTransform (F, FFLAS::FflasUpper, FFLAS::FflasUnit, m, n, R, P, Q,
                                     A.getPointer(), A.getStride(), T.getPointer(), T.getStride());
        
        delete[] P;
        delete[] Q;
        return R;
    }

        //
        // row echelonize
        //
        /**
         * \brief rowEchelonize specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t rowEchelonize (DenseMatrix<Field>& A,
                                 const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {

        size_t m = A.rowdim();
        size_t n = A.coldim();
        size_t* P = new size_t[m];
        size_t* Q = new size_t[n];

        size_t R = FFPACK::RowEchelonForm (F, m, n, A.getPointer(), A.getStride(), P, Q, false);

        FFPACK::getEchelonForm (F, FFLAS::FflasUpper, FFLAS::FflasUnit, m, n, R, Q,
                                A.getPointer(), A.getStride());
        
        delete[] P;
        delete[] Q;
        return R;
    }

        /**
         * \brief rowEchelonize with transformation specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t rowEchelonize (DenseMatrix<Field>& A, DenseMatrix<Field>& T,
                                 const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        linbox_check((A.rowdim() == T.rowdim()) && (T.rowdim() == T.coldim()));

        size_t m = A.rowdim();
        size_t n = A.coldim();
        size_t* P = new size_t[m];
        size_t* Q = new size_t[n];
        Field& F = A.getField();

        size_t R = FFPACK::RowEchelonForm (F, m, n, A.getPointer(), A.getStride(), P, Q, false);

        FFPACK::getEchelonTransform (F, FFLAS::FflasUpper, FFLAS::FflasUnit, m, n, R, P, Q,
                                     A.getPointer(), A.getStride(), T.getPointer(), T.getStride());
        FFPACK::getEchelonForm (F, FFLAS::FflasUpper, FFLAS::FflasUnit, m, n, R, Q, A.getPointer(), A.getStride());
        
        delete[] P;
        delete[] Q;
        return R;
    }

        //
        // reduced row echelon
        //

        /**
         * \brief reducedRowEchelon specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t reducedRowEchelon (DenseMatrix<Field>& E, const DenseMatrix<Field>& A,
                                     const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return reducedRowEchelon (E, A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        /**
         * \brief reducedRowEchelon with transformation specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t reducedRowEchelon (DenseMatrix<Field>& E, DenseMatrix<Field>& T, const DenseMatrix<Field>& A,
                                     const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return reducedRowEchelon (E, A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        //
        // reduced row echelonize
        //

        /**
         * \brief reducedRowEchelonize specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t reducedRowEchelonize (DenseMatrix<Field>& A,
                                        const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return reducedRowEchelonize (A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        /**
         * \brief reducedRowEchelonize with transformation specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t reducedRowEchelonize (DenseMatrix<Field>& A, DenseMatrix<Field>& T,
                                        const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return reducedRowEchelonize (A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        //
        // column echelon
        //

        /**
         * \brief colEchelon specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t colEchelon (DenseMatrix<Field>& E, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return colEchelon (E, A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        /**
         * \brief colEchelon with transformation specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t colEchelon (DenseMatrix<Field>& E, DenseMatrix<Field>& T, const DenseMatrix<Field>& A,
                              const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return colEchelon (E, A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        //
        // col echelonize
        //

        /**
         * \brief colEchelonize specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t colEchelonize (DenseMatrix<Field>& A,
                                 const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return colEchelonize (A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        /**
         * \brief colEchelonize with transformation specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t colEchelonize (DenseMatrix<Field>& A, DenseMatrix<Field>& T,
                                 const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return colEchelonize (A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        //
        // reduced column echelon
        //

        /**
         * \brief reducedColEchelon specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t reducedColEchelon (DenseMatrix<Field>& E, const DenseMatrix<Field>& A,
                                     const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return reducedColEchelon (E, A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }
        /**
         * \brief reducedColEchelon with transformation specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t reducedColEchelon (DenseMatrix<Field>& E, DenseMatrix<Field>& T, const DenseMatrix<Field>& A,
                                     const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return reducedColEchelon (E, A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        //
        // column echelonize
        //
        /**
         * \brief reducedColEchelonize specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t reducedColEchelonize (DenseMatrix<Field>& A,
                                        const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return reducedColEchelonize (A, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }

        /**
         * \brief reducedColEchelonize with transformation specialisation for DenseElimination with DenseMatrix and ModularTag.
         */
    template <class Field>
    inline size_t reducedColEchelonize (DenseMatrix<Field>& A, DenseMatrix<Field>& T,
                                        const RingCategories::ModularTag& tag, const Method::DenseElimination& m)
    {
        return reducedColEchelonize (A, T, tag, reinterpret_cast<const Method::DenseElimination&>(m));
    }
}
/* -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

