/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/gauss-rank-gf2.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <14 Jun 10 15:26:01 Jean-Guillaume.Dumas@imag.fr>
 *
 * See COPYING for license information.
 *
 * SparseElimination rank calls over GF2
 */

// =================================================================== //
// =================================================================== //
#ifndef __LINBOX_gauss_rank_gf2_INL
#define __LINBOX_gauss_rank_gf2_INL

namespace LinBox
{
	template <class SparseSeqMatrix> unsigned long&
	GaussDomain<GF2>::rankin(unsigned long &Rank,
				 SparseSeqMatrix        &A,
				 unsigned long  Ni,
				 unsigned long  Nj,
				 SparseEliminationTraits::PivotStrategy   reord)  const
	{
		Element determinant;
		const GF2 F2;
		Permutation<GF2> P((int)A.coldim(),F2);

		if (reord == SparseEliminationTraits::PIVOT_NONE)
			return NoReordering(Rank, determinant, A,  Ni, Nj);
		else
			return InPlaceLinearPivoting(Rank, determinant, A, P, Ni, Nj);
	}


	template <class SparseSeqMatrix> unsigned long&
	GaussDomain<GF2>::rankin(unsigned long &Rank,
				 SparseSeqMatrix        &A,
				 SparseEliminationTraits::PivotStrategy   reord)  const
	{
		return rankin(Rank, A,  A.rowdim (), A.coldim (), reord);
	}



	template <class SparseSeqMatrix> unsigned long&
	GaussDomain<GF2>::rank(unsigned long &rk,
			       const SparseSeqMatrix        &A,
			       SparseEliminationTraits::PivotStrategy   reord)  const
	{
		return rank(rk, A,  A.rowdim (), A.coldim (), reord);
	}

	template <class SparseSeqMatrix> unsigned long&
	GaussDomain<GF2>::rank(unsigned long &Rank,
			       const SparseSeqMatrix        &A,
			       unsigned long  Ni,
			       unsigned long  Nj,
			       SparseEliminationTraits::PivotStrategy   reord)  const
	{
		SparseSeqMatrix CopyA(Ni);
		for(unsigned long i = 0; i < Ni; ++i)
			CopyA[i] = A[i];
		return rankin(Rank, CopyA, Ni, Nj, reord);
	}
} // namespace LinBox

#endif // __LINBOX_gauss_rank_gf2_INL
