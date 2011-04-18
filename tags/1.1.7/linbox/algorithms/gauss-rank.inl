/* linbox/algorithms/gauss-rank.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <15 Jun 10 17:20:20 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 *
 * SparseElimination rank calls
 */
#ifndef __LINBOX_gauss_rank_INL
#define __LINBOX_gauss_rank_INL

namespace LinBox 
{
    template <class _Field>
    template <class Matrix> unsigned long& 
    GaussDomain<_Field>::rankin(unsigned long &rank,
                                Matrix        &A,
                                unsigned long  Ni,
                                unsigned long  Nj,
                                SparseEliminationTraits::PivotStrategy   reord)  const
    {
        Element determinant;
        if (reord == SparseEliminationTraits::PIVOT_NONE)
            return NoReordering(rank, determinant, A,  Ni, Nj);
        else
            return InPlaceLinearPivoting(rank, determinant, A, Ni, Nj);
    }

   
    template <class _Field>
    template <class Matrix> unsigned long& 
    GaussDomain<_Field>::rankin(unsigned long &rank,
                                Matrix        &A,
                                SparseEliminationTraits::PivotStrategy   reord)  const
    {
        return rankin(rank, A,  A.rowdim (), A.coldim (), reord);
    }

   

    template <class _Field>
    template <class Matrix> unsigned long& 
    GaussDomain<_Field>::rank(unsigned long &rk,
                              const Matrix        &A,
                              SparseEliminationTraits::PivotStrategy   reord)  const
    {
        return rank(rk, A,  A.rowdim (), A.coldim (), reord);
    }

    template <class _Field>
    template <class Matrix> unsigned long& 
    GaussDomain<_Field>::rank(unsigned long &rank,
                              const Matrix        &A,
                              unsigned long  Ni,
                              unsigned long  Nj,
                              SparseEliminationTraits::PivotStrategy   reord)  const
    {
        Matrix CopyA(Ni);
        for(unsigned long i = 0; i < Ni; ++i)
            CopyA[i] = A[i];
        return rankin(rank, CopyA, Ni, Nj, reord);
    }
} // namespace LinBox

#endif // __LINBOX_gauss_rank_INL
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
