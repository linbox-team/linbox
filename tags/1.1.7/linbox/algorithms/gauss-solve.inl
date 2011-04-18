/* linbox/algorithms/gauss-solve.inl
 * Copyright (C) LinBox 2008 
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <21 Jan 10 17:01:22 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_gauss_solve_INL
#define __LINBOX_gauss_solve_INL

#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/triangular-solve.h"
#include "linbox/blackbox/permutation.h"

namespace LinBox 
{


    template <class _Field>
    template <class Matrix, class Perm, class Vector1, class Vector2> inline Vector1& 
    GaussDomain<_Field>::solve(Vector1& x, Vector1& w, unsigned long rank, const Perm& Q, const Matrix& L, const Matrix& U, const Perm& P, const Vector2& b)  const {
        
        Vector2 y(U.rowdim()), v(U.rowdim());
        
        Q.applyTranspose(y, b);
        
        lowerTriangularUnitarySolve(v, L, y);
        
        upperTriangularSolve(w, U, v);

        return P.applyTranspose(x, w);
    }

    template <class _Field>
    template <class Matrix, class Perm, class Vector1, class Vector2> inline Vector1& 
    GaussDomain<_Field>::solve(Vector1& x, unsigned long rank, const Perm& Q, const Matrix& L, const Matrix& U, const Perm& P, const Vector2& b, bool randomsol)  const {
        
        Vector1 w(U.coldim());
        if (randomsol) {
                // Random solution is in output
            typename _Field::RandIter generator(_F);
            for(typename Vector1::iterator it=w.begin()+rank;it!=w.end();++it)
                generator.random( *it );
        } else {
            for(typename Vector1::iterator it=w.begin()+rank;it!=w.end();++it)
                _F.init(*it,0);
        }
        return this->solve(x, w, rank, Q, L, U, P, b);
    }

    template <class _Field>
    template <class Matrix, class Vector1, class Vector2> inline Vector1& 
    GaussDomain<_Field>::solvein(Vector1& x, Matrix& A, const Vector2& b, bool randomsol)  const {
        
        typename Field::Element det;
        unsigned long rank;
        Matrix L(_F, A.rowdim(), A.rowdim());
        Permutation<Field> Q(A.rowdim(),_F);
        Permutation<Field> P(A.coldim(),_F);

        this->QLUPin(rank, det, Q, L, A, P, A.rowdim(), A.coldim() );

        if (! randomsol) {
                // Sets solution values to 0 for coldim()-rank columns
                // Therefore, prune unnecessary elements 
                // in those last columns of U
            for(typename Matrix::RowIterator row=A.rowBegin();
                row != A.rowEnd(); ++row) {
                if (row->size()) {
                    size_t ns=0;
                    for(typename Matrix::Row::iterator it = row->begin();
                        it != row->end(); ++it, ++ns) {
                        if (it->first >= rank) {
                            row->resize(ns);
                            break;
                        }
                    }
                }
            }
        }

        return this->solve(x, rank, Q, L, A, P, b, randomsol);
    }

} // namespace LinBox

#endif // __LINBOX_gauss_solve_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
