/* linbox/algorithms/gauss-solve-gf2.inl
 * Copyright (C) LinBox 2009
 *
 * Written by Jean-Guillaume Dumas <Jean-Guillaume.Dumas@imag.fr>
 * Time-stamp: <21 Jan 10 15:42:14 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 */

#ifndef __GAUSS_SOLVE_GF2_INL
#define __GAUSS_SOLVE_GF2_INL

#include "linbox/algorithms/gauss-gf2.h"
#include "linbox/algorithms/triangular-solve-gf2.h"
#include "linbox/blackbox/permutation.h"

namespace LinBox 
{


    template <class SparseSeqMatrix, class Perm, class Vector1, class Vector2>
    Vector1& GaussDomain<GF2>::solve(Vector1& x, Vector1& w, unsigned long rank, 
                   const Perm& Q, const SparseSeqMatrix& L, 
                   const SparseSeqMatrix& U, const Perm& P, 
                   const Vector2& b) const {
        
        Vector2 y(U.rowdim()), v(U.rowdim());

        Q.applyTranspose(y, b);

        lowerTriangularUnitarySolveBinary(v, L, y);
       
        upperTriangularSolveBinary(w, U, v);

        return P.applyTranspose(x, w);
    }

    template <class SparseSeqMatrix, class Perm, class Vector1, class Vector2>
    Vector1& GaussDomain<GF2>::solve(Vector1& x, unsigned long rank, 
                   const Perm& Q, const SparseSeqMatrix& L, 
                   const SparseSeqMatrix& U, const Perm& P, 
                   const Vector2& b) const {
        
        Vector1 w(U.coldim());
        const GF2 F2;
        typename GF2::RandIter generator(F2);
        for(typename Vector1::iterator it=w.begin()+rank;it!=w.end();++it)
            generator.random( *it );

        return this->solve(x, w, rank, Q, L, U, P, b);
    }

    template <class SparseSeqMatrix, class Vector1, class Vector2>
    Vector1& GaussDomain<GF2>::solvein(Vector1& x,
                     SparseSeqMatrix        &A,
                     const Vector2& b) const {
        
        typename GF2::Element det;
        unsigned long rank;
        const GF2 F2;
        SparseSeqMatrix L(F2, A.rowdim(), A.rowdim());
        Permutation<GF2> Q(A.rowdim(),F2);
        Permutation<GF2> P(A.coldim(),F2);
        
        
        
        this->QLUPin(rank, det, Q, L, A, P, A.rowdim(), A.coldim() );

        return this->solve(x, rank, Q, L, A, P, b);
    }

} // namespace LinBox

#endif // __GAUSS_SOLVE_GF2_INL
