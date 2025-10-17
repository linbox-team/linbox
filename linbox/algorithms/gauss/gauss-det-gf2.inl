/* linbox/algorithms/gauss-det-gf2.inl
 * Copyright (C) 2009 The LinBox group
 *
// Time-stamp: <17 Oct 25 13:47:14 Jean-Guillaume.Dumas@imag.fr>
*
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 *
 * SparseElimination determinant calls gf2 specialization
 */
#ifndef __LINBOX_gauss_det_gf2_INL
#define __LINBOX_gauss_det_gf2_INL

namespace LinBox
{
	template <class SparseSeqMatrix> inline typename GaussDomain<GF2>::Element&
	GaussDomain<GF2>::detInPlace(Element		&determinant,
                                 SparseSeqMatrix	&A,
                                 size_t	Ni,
                                 size_t	Nj,
                                 PivotStrategy   reord)  const
	{
		size_t Rank;
		if (reord == PivotStrategy::None)
			NoReordering(Rank, determinant, A, Ni, Nj);
		else {
            Permutation<GF2> P(A.field(),(int)A.coldim());
			InPlaceLinearPivoting(Rank, determinant, A, P, Ni, Nj);
        }
		return determinant;
	}


	template <class SparseSeqMatrix> inline typename GaussDomain<GF2>::Element&
	GaussDomain<GF2>::detInPlace(Element &determinant,
                                 SparseSeqMatrix  &A,
                                 PivotStrategy   reord)  const
	{
		return detInPlace(determinant, A,  A.rowdim (), A.coldim (), reord);
	}



	template <class SparseSeqMatrix> inline typename GaussDomain<GF2>::Element&
	GaussDomain<GF2>::det(Element        &determinant,
                          const SparseSeqMatrix   &A,
                          PivotStrategy   reord)  const
	{
		return det(determinant, A,  A.rowdim (), A.coldim (), reord);
	}

	template <class SparseSeqMatrix> inline typename GaussDomain<GF2>::Element&
	GaussDomain<GF2>::det(Element       &determinant,
                          const SparseSeqMatrix  &A,
                          size_t  Ni,
                          size_t  Nj,
                          PivotStrategy   reord)  const
	{
		SparseSeqMatrix CopyA(Ni);
		for(size_t i = 0; i < Ni; ++i)
			CopyA[i] = A[i];
		return detInPlace(determinant, CopyA, Ni, Nj, reord);
	}
} // namespace LinBox

#endif // __LINBOX_gauss_det_gf2_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
