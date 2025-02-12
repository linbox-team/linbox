/* linbox/algorithms/gauss-det.inl
 * Copyright (C) 2009 The LinBox group
 *
// Time-stamp: <24 Aug 17 18:19:58 Jean-Guillaume.Dumas@imag.fr>
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
 * SparseElimination determinant calls
 */
#ifndef __LINBOX_gauss_det_INL
#define __LINBOX_gauss_det_INL

namespace LinBox
{
	template <class _Field>
	template <class _Matrix> inline typename GaussDomain<_Field>::Element&
	GaussDomain<_Field>::detInPlace(Element        &determinant,
				   _Matrix        &A,
				   size_t  Ni,
				   size_t  Nj,
				   PivotStrategy   reord)  const
	{
		size_t Rank;
		if (reord == PivotStrategy::None)
			NoReordering(Rank, determinant, A,  Ni, Nj);
		else
			InPlaceLinearPivoting(Rank, determinant, A, Ni, Nj);
		return determinant;
	}


	template <class _Field>
	template <class _Matrix> inline typename GaussDomain<_Field>::Element&
	GaussDomain<_Field>::detInPlace(Element &determinant,
				   _Matrix  &A,
				   PivotStrategy   reord)  const
	{
		return detInPlace(determinant, A,  A.rowdim (), A.coldim (), reord);
	}



	template <class _Field>
	template <class _Matrix> inline typename GaussDomain<_Field>::Element&
	GaussDomain<_Field>::det(Element        &determinant,
				 const _Matrix   &A,
				 PivotStrategy   reord)  const
	{
		return det(determinant, A,  A.rowdim (), A.coldim (), reord);
	}

	template <class _Field>
	template <class _Matrix> inline typename GaussDomain<_Field>::Element&
	GaussDomain<_Field>::det(Element       &determinant,
				 const _Matrix  &A,
				 size_t  Ni,
				 size_t  Nj,
				 PivotStrategy   reord)  const
	{
		_Matrix CopyA(Ni);
		for(size_t i = 0; i < Ni; ++i)
			CopyA[i] = A[i];
		return detInPlace(determinant, CopyA, Ni, Nj, reord);
	}
} // namespace LinBox

#endif // __LINBOX_gauss_det_INL

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
