/* linbox/algorithms/gauss-det.inl
 * Copyright (C) 2009 The LinBox group
 *
// Time-stamp: <15 Jun 10 17:20:08 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 *
 * SparseElimination determinant calls
 */
#ifndef __LINBOX_gauss_det_INL
#define __LINBOX_gauss_det_INL

namespace LinBox 
{
    template <class _Field>
    template <class Matrix> inline typename GaussDomain<_Field>::Element& 
    GaussDomain<_Field>::detin(Element        &determinant,
                               Matrix        &A,
                               unsigned long  Ni,
                               unsigned long  Nj,
                               SparseEliminationTraits::PivotStrategy   reord)  const
    {
        unsigned long rank;
        if (reord == SparseEliminationTraits::PIVOT_NONE)
            NoReordering(rank, determinant, A,  Ni, Nj);
        else
            InPlaceLinearPivoting(rank, determinant, A, Ni, Nj);
        return determinant;
    }

   
    template <class _Field>
    template <class Matrix> inline typename GaussDomain<_Field>::Element& 
    GaussDomain<_Field>::detin(Element &determinant,
                               Matrix  &A,
                               SparseEliminationTraits::PivotStrategy   reord)  const
    {
        return detin(determinant, A,  A.rowdim (), A.coldim (), reord);
    }

   

    template <class _Field>
    template <class Matrix> inline typename GaussDomain<_Field>::Element& 
    GaussDomain<_Field>::det(Element        &determinant,
                             const Matrix   &A,
                             SparseEliminationTraits::PivotStrategy   reord)  const
    {
        return det(determinant, A,  A.rowdim (), A.coldim (), reord);
    }

    template <class _Field>
    template <class Matrix> inline typename GaussDomain<_Field>::Element& 
    GaussDomain<_Field>::det(Element       &determinant,
                             const Matrix  &A,
                             unsigned long  Ni,
                             unsigned long  Nj,
                             SparseEliminationTraits::PivotStrategy   reord)  const
    {
        Matrix CopyA(Ni);
        for(unsigned long i = 0; i < Ni; ++i)
            CopyA[i] = A[i];
        return detin(determinant, CopyA, Ni, Nj, reord);
    }
} // namespace LinBox

#endif // __LINBOX_gauss_det_INL
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
