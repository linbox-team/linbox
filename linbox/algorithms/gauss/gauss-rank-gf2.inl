/* linbox/algorithms/gauss-rank-gf2.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <13 Nov 17 17:00:08 Jean-Guillaume.Dumas@imag.fr>
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
 * SparseElimination rank calls over GF2
 */

// =================================================================== //
// =================================================================== //
#ifndef __LINBOX_gauss_rank_gf2_INL
#define __LINBOX_gauss_rank_gf2_INL

namespace LinBox
{
	template <class SparseSeqMatrix> size_t&
	GaussDomain<GF2>::rankInPlace(size_t &Rank,
				 SparseSeqMatrix        &A,
				 size_t  Ni,
				 size_t  Nj,
				 PivotStrategy   reord)  const
	{
		Element determinant;
		const GF2 F2;
		Permutation<GF2> P(F2,(int)A.coldim());

		if (reord == PivotStrategy::None)
			return NoReordering(Rank, determinant, A,  Ni, Nj);
		else
			return InPlaceLinearPivoting(Rank, determinant, A, P, Ni, Nj);
	}


	template <class SparseSeqMatrix> size_t&
	GaussDomain<GF2>::rankInPlace(size_t &Rank,
				 SparseSeqMatrix        &A,
				 PivotStrategy   reord)  const
	{
		return rankInPlace(Rank, A,  A.rowdim (), A.coldim (), reord);
	}



	template <class SparseSeqMatrix> size_t&
	GaussDomain<GF2>::rank(size_t &rk,
			       const SparseSeqMatrix        &A,
			       PivotStrategy   reord)  const
	{
		return rank(rk, A,  A.rowdim (), A.coldim (), reord);
	}

	template <class SparseSeqMatrix> size_t&
	GaussDomain<GF2>::rank(size_t &Rank,
			       const SparseSeqMatrix        &A,
			       size_t  Ni,
			       size_t  Nj,
			       PivotStrategy   reord)  const
	{
		SparseSeqMatrix CopyA(Ni);
		for(size_t i = 0; i < Ni; ++i)
			CopyA[i] = A[i];
		return rankInPlace(Rank, CopyA, Ni, Nj, reord);
	}
} // namespace LinBox

#endif // __LINBOX_gauss_rank_gf2_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
