// =================================================================== //
// SparseElimination determinant calls
// Time-stamp: <21 Jan 10 15:07:05 Jean-Guillaume.Dumas@imag.fr> 
// =================================================================== //
#ifndef __GAUSS_DET_INL
#define __GAUSS_DET_INL

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
            InPlaceLinearPivoting(rank, determinant, A,  Ni, Nj);
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

#endif // __GAUSS_INL
