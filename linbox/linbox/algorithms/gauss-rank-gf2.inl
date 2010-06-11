/* linbox/algorithms/gauss-rank-gf2.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <21 Jan 10 15:08:58 Jean-Guillaume.Dumas@imag.fr> 
 *
 * See COPYING for license information.
 *
 * SparseElimination rank calls over GF2
 */

// =================================================================== //
// =================================================================== //
#ifndef __GAUSS_RANK_GF2_INL
#define __GAUSS_RANK_GF2_INL

namespace LinBox 
{
    template <class SparseSeqMatrix> unsigned long& 
    GaussDomain<GF2>::rankin(unsigned long &rank,
                                SparseSeqMatrix        &A,
                                unsigned long  Ni,
                                unsigned long  Nj,
                                SparseEliminationTraits::PivotStrategy   reord)  const
    {
        Element determinant;
        if (reord == SparseEliminationTraits::PIVOT_NONE)
            return NoReordering(rank, determinant, A,  Ni, Nj);
        else
            return InPlaceLinearPivoting(rank, determinant, A,  Ni, Nj);
    }

   
    template <class SparseSeqMatrix> unsigned long& 
    GaussDomain<GF2>::rankin(unsigned long &rank,
                                SparseSeqMatrix        &A,
                                SparseEliminationTraits::PivotStrategy   reord)  const
    {
        return rankin(rank, A,  A.rowdim (), A.coldim (), reord);
    }

   

    template <class SparseSeqMatrix> unsigned long& 
    GaussDomain<GF2>::rank(unsigned long &rk,
                              const SparseSeqMatrix        &A,
                              SparseEliminationTraits::PivotStrategy   reord)  const
    {
        return rank(rk, A,  A.rowdim (), A.coldim (), reord);
    }

    template <class SparseSeqMatrix> unsigned long& 
    GaussDomain<GF2>::rank(unsigned long &rank,
                              const SparseSeqMatrix        &A,
                              unsigned long  Ni,
                              unsigned long  Nj,
                              SparseEliminationTraits::PivotStrategy   reord)  const
    {
        SparseSeqMatrix CopyA(Ni);
        for(unsigned long i = 0; i < Ni; ++i)
            CopyA[i] = A[i];
        return rankin(rank, CopyA, Ni, Nj, reord);
    }
} // namespace LinBox

#endif // __GAUSS_GF2_INL
