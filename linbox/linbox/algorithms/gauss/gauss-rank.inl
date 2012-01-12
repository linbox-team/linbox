/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/gauss-rank.inl
 * Copyright (C) 2009 The LinBox group
 *
 * Time-stamp: <15 Jun 10 17:20:20 Jean-Guillaume.Dumas@imag.fr>
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
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 *
 * SparseElimination rank calls
 */
#ifndef __LINBOX_gauss_rank_INL
#define __LINBOX_gauss_rank_INL

namespace LinBox
{
	template <class _Field>
	template <class Matrix> unsigned long&
	GaussDomain<_Field>::rankin(unsigned long &Rank,
				    Matrix        &A,
				    unsigned long  Ni,
				    unsigned long  Nj,
				    SparseEliminationTraits::PivotStrategy   reord)  const
	{
		Element determinant;
		if (reord == SparseEliminationTraits::PIVOT_NONE)
			return NoReordering(Rank, determinant, A,  Ni, Nj);
		else
			return InPlaceLinearPivoting(Rank, determinant, A, Ni, Nj);
	}


	template <class _Field>
	template <class Matrix> unsigned long&
	GaussDomain<_Field>::rankin(unsigned long &Rank,
				    Matrix        &A,
				    SparseEliminationTraits::PivotStrategy   reord)  const
	{
		return rankin(Rank, A,  A.rowdim (), A.coldim (), reord);
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
	GaussDomain<_Field>::rank(unsigned long &Rank,
				  const Matrix        &A,
				  unsigned long  Ni,
				  unsigned long  Nj,
				  SparseEliminationTraits::PivotStrategy   reord)  const
	{
		Matrix CopyA(Ni);
		for(unsigned long i = 0; i < Ni; ++i)
			CopyA[i] = A[i];
		return rankin(Rank, CopyA, Ni, Nj, reord);
	}
} // namespace LinBox

#endif // __LINBOX_gauss_rank_INL
