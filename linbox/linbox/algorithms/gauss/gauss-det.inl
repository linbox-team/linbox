/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/algorithms/gauss-det.inl
 * Copyright (C) 2009 The LinBox group
 *
// Time-stamp: <15 Jun 10 17:20:08 Jean-Guillaume.Dumas@imag.fr>
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
		unsigned long Rank;
		if (reord == SparseEliminationTraits::PIVOT_NONE)
			NoReordering(Rank, determinant, A,  Ni, Nj);
		else
			InPlaceLinearPivoting(Rank, determinant, A, Ni, Nj);
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
