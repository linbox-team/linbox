/* linbox/algorithms/gauss-solve.inl
 * Copyright (C) LinBox 2008 
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <15 Sep 08 14:45:40 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 */

#ifndef __GAUSS_SOLVE_INL
#define __GAUSS_SOLVE_INL

#include "linbox/algorithms/gauss.h"
#include "linbox/algorithms/triangular-solve.h"
#include "linbox/blackbox/permutation.h"

namespace LinBox 
{


    template <class _Field>
    template <class Matrix, class Perm, class Vector1, class Vector2> inline Vector1& 
    GaussDomain<_Field>::solve(Vector1& x, unsigned long rank, const Perm& Q, const Matrix& L, const Matrix& U, const Perm& P, const Vector2& b) {
        
        Vector2 y(U.rowdim()), v(U.rowdim());
        Vector1 w(U.coldim());
        typename _Field::RandIter generator(_F);
        for(typename Vector1::iterator it=w.begin()+rank;it!=w.end();++it)
            generator.random( *it );
        
        Q.applyTranspose(y, b);
        
        lowerTriangularUnitarySolve(v, L, y);
        
        upperTriangularSolve(w, U, v);

        return P.applyTranspose(x, w);
    }

    template <class _Field>
    template <class Matrix, class Vector1, class Vector2> inline Vector1& 
    GaussDomain<_Field>::solvein(Vector1& x, Matrix& A, const Vector2& b) {
        
        typename Field::Element det;
        unsigned long rank;
        Matrix L(_F, A.rowdim(), A.rowdim());
        Permutation<Field> Q(A.rowdim(),_F);
        Permutation<Field> P(A.coldim(),_F);

        this->QLUPin(rank, det, Q, L, A, P, A.rowdim(), A.coldim() );

        return this->solve(x, rank, Q, L, A, P, b);
    }

} // namespace LinBox

#endif // __GAUSS_SOLVE_INL
