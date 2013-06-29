/* Copyright (C) 2010 LinBox
 * Written by <brice.boyer@imag.fr>
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/** \file algorithms/dense-nullspace.h
 * @ingroup algorithm
 *
 * @brief We provide the right or left nullspace (kernel or cokernel) of a dense matrix.
 * @details Provides :
 *	- the nullspace of a matrix \p A
 *	- (soon) a random vector within the nullspace of \p A
 * @todo random nullspace vector
 */

#ifndef __LINBOX_dense_nullspace_H
#define __LINBOX_dense_nullspace_H


#include "linbox/linbox-tags.h"
#include "linbox/matrix/blas-matrix.h"

namespace LinBox
{

	/** Computes the kernel of a dense matrix using \c LQUP.
	 *
	 * Acccording to the dimensions of the input matrix, we chose different methods.
	 * @warning timings may vary and these choices were made on an experimental basis.
	 *
	 * @param F  Field
	 * @param Side  left or right from \c LinBox::SideTag
	 * @param m rows
	 * @param n cols
	 * @param A input matrix
	 * @param lda leading dimension of A
	 * @param Ker Kernel. \c NULL if \c kerdim==0
	 * @param ldk leading dimension of the kernel.
	 * @param kerdim dimension of the kernel.
	 * @return dimension of the kernel.
	 *
	 * @warning A is modified.
	 */
	template<class Field>
	size_t
	NullSpaceBasis (const Field& F, const Tag::Side Side,
			const size_t & m, const size_t & n,
			typename Field::Element * A, const size_t & lda,
			typename Field::Element *& Ker, size_t& ldk,
			size_t & kerdim) ;

	/*! Nullspace of a dense matrix on a finite field.
	 * A is modified.
	 * @param         F field
	 * @param         Side \c SideTag::Left or \c SideTag::Right nullspace.
	 * @param[in,out] A Input matrix
	 * @param[out]    Ker Nullspace of the matrix (Allocated in the routine)
	 * @param[out]    kerdim rank of the kernel
	 * @return \p kerdim
	 */
	template<class Field>
	size_t&
	NullSpaceBasis (const Field& F, const Tag::Side Side,
			BlasMatrix<Field> & A,
			BlasMatrix<Field> & Ker,
			size_t & kerdim) ;

	/*! Nullspace of a dense matrix on a finite field.
	 * A is preserved.
	 * @param      F field
	 * @param      Side \c SideTag::Left or \c SideTag::Right nullspace.
	 * @param[in]  A Input matrix
	 * @param[out] Ker Nullspace of the matrix (Allocated in the routine)
	 * @param[out] kerdim rank of the kernel
	 * @return \p kerdim
	 */
	template<class Field>
	size_t&
	NullSpaceBasis (const Field& F, const Tag::Side Side,
			const BlasMatrix<Field> & A,
			BlasMatrix<Field> & Ker,
			size_t & kerdim)
	{
		BlasMatrix<Field> B (A);
		return NullSpaceBasis<Field>(F,Side,B,Ker,kerdim);

	}


} // LinBox

#include "dense-nullspace.inl"

#endif // __LINBOX_dense_nullspace_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
