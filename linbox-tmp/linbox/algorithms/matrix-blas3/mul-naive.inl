/*  Copyright (C) 2012 the members of the LinBox group
 * Written by B. Boyer < bboyer@imag.fr >
 *
 * This file is part of the LinBox library.
 *
 * ========LICENCE========
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * LinBox is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *
 */

#ifndef __LINBOX_matrix_blas3_mul_naive_INL
#define __LINBOX_matrix_blas3_mul_naive_INL


namespace LinBox { namespace BLAS3 {
	template<class _anyMatrix, class _otherMatrix1, class _otherMatrix2>
	_anyMatrix & mul (_anyMatrix& C,
			  const _otherMatrix1& A,
			  const _otherMatrix2& B,
			  const mulMethod::naive &)
	{
		// TODO check sizes
		// TODO check fields
		// TODO check get/set Entry
		typedef typename _anyMatrix::Field Field;
		const Field &F = B.field();
		for (size_t i = 0 ; i <C.rowdim(); ++i)
			for (size_t j = 0 ; j <C.coldim(); ++j) {
				C.setEntry(i,j,F.zero);
				for (size_t k = 0 ; k <B.rowdim(); ++k)
					F.axpyin(C.refEntry(i,j),
						 A.getEntry(i, k),
						 B.getEntry(k, j));
			}
		return C;
	}
} // BLAS3
} // LinBox


#endif //  __LINBOX_matrix_blas3_mul_naive_INL

//Local Variables:
//mode: C++
//tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

