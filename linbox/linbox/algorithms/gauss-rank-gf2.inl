// =================================================================== //
// SparseElimination rank calls over GF2
// Time-stamp: <01 Sep 09 15:31:30 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef __GAUSS_RANK_GF2_INL
#define __GAUSS_RANK_GF2_INL

namespace LinBox 
{
    template <class Matrix> unsigned long& 
    GaussDomain<GF2>::rankin(unsigned long &rank,
                                Matrix        &A,
                                unsigned long  Ni,
                                unsigned long  Nj,
                                SparseEliminationTraits::PivotStrategy   reord) 
    {
        Element determinant;
        if (reord == SparseEliminationTraits::PIVOT_NONE)
            return NoReordering(rank, determinant, A,  Ni, Nj);
        else
            return InPlaceLinearPivoting(rank, determinant, A,  Ni, Nj);
    }

   
    template <class Matrix> unsigned long& 
    GaussDomain<GF2>::rankin(unsigned long &rank,
                                Matrix        &A,
                                SparseEliminationTraits::PivotStrategy   reord) 
    {
        return rankin(rank, A,  A.rowdim (), A.coldim (), reord);
    }

   

    template <class Matrix> unsigned long& 
    GaussDomain<GF2>::rank(unsigned long &rank,
                              const Matrix        &A,
                              SparseEliminationTraits::PivotStrategy   reord) 
    {
        return rank(rank, A,  A.rowdim (), A.coldim (), reord);
    }

    template <class Matrix> unsigned long& 
    GaussDomain<GF2>::rank(unsigned long &rank,
                              const Matrix        &A,
                              unsigned long  Ni,
                              unsigned long  Nj,
                              SparseEliminationTraits::PivotStrategy   reord) 
    {
        Matrix CopyA(Ni);
        for(unsigned long i = 0; i < Ni; ++i)
            CopyA[i] = A[i];
        return rankin(rank, CopyA, Ni, Nj, reord);
    }
} // namespace LinBox

#endif // __GAUSS_GF2_INL
