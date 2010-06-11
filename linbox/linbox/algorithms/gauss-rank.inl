// =================================================================== //
// SparseElimination rank calls
// Time-stamp: <21 Jan 10 15:07:06 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef __GAUSS_RANK_INL
#define __GAUSS_RANK_INL

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
            return InPlaceLinearPivoting(rank, determinant, A,  Ni, Nj);
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

#endif // __GAUSS_INL
