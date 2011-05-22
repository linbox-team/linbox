/* ===================================================================
 * Copyright(C) LinBox 2008
 * Written by Jean-Guillaume Dumas
 * Triangular Solve
 * See COPYING for license information.
 * Time-stamp: <01 Oct 09 15:38:25 Jean-Guillaume.Dumas@imag.fr> 
 * ===================================================================
 */
#ifndef __LINBOX_tri_solve_gf2_INL
#define __LINBOX_tri_solve_gf2_INL

#include "linbox/vector/vector-domain.h"
#include "linbox/field/gf2.h"

namespace LinBox 
{
    template <class _Matrix, class Vector1, class Vector2> Vector1&
    upperTriangularSolveBinary (Vector1& x,
                          const _Matrix  &U,
                          const Vector2& b)
    {
//         linbox_check( x.size() == U.coldim() );
//         linbox_check( b.size() == U.rowdim() );
        typedef _Matrix Matrix;
        typedef GF2 Field;
        const GF2 F2;

        commentator.start ("Sparse Elimination Upper Triangular Solve over GF(2)", "utrsmGF2");

        typename Vector2::const_iterator vec=b.begin();
        typename Vector1::iterator res=x.begin();
        typename Matrix::const_iterator row=U.begin();

            // Find last constrained values of x, U and b
//         for( ; (res != x.end()) && (row != U.rowEnd()); ++res, ++row, ++vec) { }
        size_t last = x.size();
        if( b.size() < last ) last = b.size();
        res += last;
        row += last;
        vec += last;

        VectorCategories::DenseZeroOneVectorTag  DZOtag;
        VectorCategories::SparseZeroOneVectorTag SZOtag;
        
        bool consistant = true;
        for(typename Vector2::const_iterator bcheck=vec; bcheck != b.end(); ++bcheck) {
            if( ! F2.isZero(*bcheck) ) {
                consistant = false;
                break;
            }
        }
        if (consistant) {
            --vec; --res; --row;
            
            VectorDomain<Field> VD(F2);
            for( ; row != U.begin(); --row, --vec, --res) {
                F2.init(*res, 0UL);
                if (row->size()) {
                    typename Field::Element tmp;
                    VD.dotSpecialized(tmp, x, *row, DZOtag, SZOtag);
                    F2.addin(tmp,*vec);
                    F2.assign(*res,tmp);
                } else {
                        // Consistency check
                    if( ! F2.isZero(*vec) ) {
                        consistant = false;
                        break;
                    }
                }
            }
            F2.init(*res, 0UL);
            if (row->size()) {
                typename Field::Element tmp;
                VD.dotSpecialized(tmp, x, *row, DZOtag, SZOtag);
                F2.addin(tmp,*vec);
                F2.assign(*res,tmp);
            } else {
                    // Consistency check
                if( ! F2.isZero(*vec) ) consistant = false;
            }
        }
//         if (! consistant) throw LinboxError ("upperTriangularSolveBinary returned INCONSISTENT");
        linbox_check( consistant );
        
        commentator.stop ("done", NULL, "utrsmGF2");
        return x;
    }

    template <class _Matrix, class Vector1, class Vector2> Vector1&
    lowerTriangularUnitarySolveBinary (Vector1& x,
                                 const _Matrix  &L,
                                 const Vector2& b)
    {
        linbox_check( b.size() == L.coldim() );
        typedef _Matrix Matrix;
        const GF2 F2;

        commentator.start ("Sparse Elimination Lower Triangular Unitary Solve over GF2", "ltrsmGF2");

        typename Vector2::const_iterator vec=b.begin();
        typename Vector1::iterator res=x.begin();
        typename Matrix::const_iterator row=L.begin();

        VectorCategories::DenseZeroOneVectorTag  DZOtag;
        VectorCategories::SparseZeroOneVectorTag SZOtag;
        VectorDomain<GF2> VD(F2);
        for( ; row != L.end(); ++row, ++vec, ++res) {
            F2.init(*res, 0UL);
            GF2::Element tmp;
            VD.dotSpecialized(tmp, *row, x, SZOtag, DZOtag);
            F2.negin(tmp);
            F2.addin(tmp,*vec);
            F2.assign(*res,tmp);
        }
        
        commentator.stop ("done", NULL, "ltrsmGF2");
        return x;
    }

}
#endif //__LINBOX_tri_solve_gf2_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
