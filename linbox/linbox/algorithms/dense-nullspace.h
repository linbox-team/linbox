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
#include "linbox/matrix/dense-matrix.h"

namespace LinBox
{



	/*! Nullspace of a dense matrix on a finite field.
	 * A is modified.
	 * @param         F field
	 * @param         Side \c SideTag::Left or \c SideTag::Right nullspace.
	 * @param[in,out] A Input matrix
	 * @param[out]    Ker Nullspace of the matrix (Allocated in the routine)
	 * @param[out]    kerdim rank of the kernel
	 * @return \p kerdim
	 *
	 * @todo make it work for BlasSubmatrix too
	 */
	template<class Field>
	size_t&
	NullSpaceBasisIn (const LINBOX_enum(Tag::Side) Side,
			BlasMatrix<Field> & A,
			BlasMatrix<Field> & Ker,
			size_t & kerdim) ;

	template<class DenseMat>
	size_t&
	NullSpaceBasisIn (const LINBOX_enum(Tag::Side) Side,
			BlasSubmatrix<DenseMat> & A,
			BlasMatrix<typename DenseMat::Field> & Ker,
			size_t & kerdim) ;


	/*! Nullspace of a dense matrix on a finite field.
	 * A is preserved.
	 * @param      F field
	 * @param      Side \c SideTag::Left or \c SideTag::Right nullspace.
	 * @param[in]  A Input matrix
	 * @param[out] Ker Nullspace of the matrix (Allocated in the routine)
	 * @param[out] kerdim rank of the kernel
	 * @return \p kerdim
	 *
	 * @todo make it work for BlasSubmatrix too
	 */
	template<class Field>
	size_t&
	NullSpaceBasis (const LINBOX_enum(Tag::Side) Side,
			const BlasMatrix<Field> & A,
			BlasMatrix<Field> & Ker,
			size_t & kerdim) ;

	template<class DenseMat>
	size_t&
	NullSpaceBasis (const LINBOX_enum(Tag::Side) Side,
			const BlasSubmatrix<DenseMat> & A,
			BlasMatrix<typename DenseMat::Field> & Ker,
			size_t & kerdim);



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
