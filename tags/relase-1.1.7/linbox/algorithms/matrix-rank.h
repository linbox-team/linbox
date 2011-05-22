/* Copyright (C) 2003 LinBox
 *  Author: Zhendong Wan
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __LINBOX_matrix_rank_H
#define __LINBOX_matrix_rank_H

#include <linbox/util/debug.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/sparse.h>
#include <linbox/solutions/rank.h>

#include <linbox/algorithms/matrix-hom.h>
#include <vector>
#include <algorithm>
#include <linbox/randiter/random-prime.h>

namespace LinBox 
{    

	/** Compute the rank of an integer matrix in place over a finite field by Gaussian elimination.
	 */
	template<class _Ring, class _Field, class _RandomPrime = RandomPrimeIterator>
	class MatrixRank {

	public:

		typedef _Ring Ring;
		typedef _Field Field;

		Ring r;

		mutable _RandomPrime rp;

		MatrixRank(const Ring& _r = Ring(), const _RandomPrime& _rp = _RandomPrime() ) : r(_r), rp (_rp) {}

		~MatrixRank() {}

		//compute the integer matrix A by modulo a random prime, Monto-Carlo
		template<class IMatrix>
		long rank(const IMatrix& A) const {

			Field F (*rp);

			DenseMatrix<Field> Ap(F, A.rowdim(), A.coldim());

			MatrixHom::map(Ap, A, F);

			long result;

			result = rankIn(Ap);

			return result;
		}

		template <class Row>
		long rank(const SparseMatrix<Ring, Row>& A) const {

			Field F (*rp);
			typename SparseMatrix<Ring, Row>::template rebind<Field>::other Ap(A, F);
			long result;
			result = rankIn (Ap);
			return result;
		}


		template<class Field, class Row>
		long rankIn(SparseMatrix<Field, Row>& A) const {

			unsigned long result;

			LinBox::rank(result, A, A.field());

			return result;
		}

		// compute rank by Gauss Elimination
		long rankIn(DenseMatrix<Field>& Ap) const {

			typedef typename Field::Element Element;

			Field F = Ap.field();

			typename DenseMatrix<Field>::RowIterator cur_r, tmp_r;
			typename DenseMatrix<Field>::ColIterator cur_c, tmp_c;
			typename DenseMatrix<Field>::Row::iterator cur_rp, tmp_rp;
			typename DenseMatrix<Field>::Col::iterator tmp_cp;

			Element tmp_e;

			std::vector<Element> tmp_v(Ap.coldim());

			int offset_r = 0;

			int offset_c = 0;

			int r = 0;

			for(cur_r = Ap. rowBegin(), cur_c = Ap. colBegin(); (cur_r != Ap. rowEnd())&&(cur_c != Ap.colEnd());) {

				//try to find the pivot.
				tmp_r = cur_r;

				tmp_cp = cur_c -> begin() + offset_c;

				while ((tmp_cp != cur_c -> end()) && F.isZero(*tmp_cp)) {
					++ tmp_cp;
					++ tmp_r;
				}

				// if no pivit found
				if (tmp_cp == cur_c -> end()) {
					++ offset_r;
					++ cur_c;
					continue;
				}

				//if swicth two row if nessary. Each row in dense matrix is stored in contiguous space
				if (tmp_r != cur_r) { 

					std::copy (tmp_r -> begin(), tmp_r -> end(), tmp_v.begin());

					std::copy (cur_r -> begin(), cur_r -> end(), tmp_r -> begin());

					std::copy (tmp_v.begin(), tmp_v.end(), cur_r -> begin());
				}

				// continue gauss elimination	 
				for(tmp_r = cur_r + 1; tmp_r != Ap.rowEnd(); ++ tmp_r) {	   

					//see if need to update the row
					if (!F.isZero(*(tmp_r -> begin() + offset_r ))) {

						F.div (tmp_e, *(tmp_r -> begin() + offset_r), *(cur_r -> begin() + offset_r));

						F.negin(tmp_e);		    

						for ( cur_rp = cur_r ->begin() + offset_r,tmp_rp =  tmp_r -> begin() + offset_r; 
						      tmp_rp != tmp_r -> end(); ++ tmp_rp, ++ cur_rp )

							F.axpyin ( *tmp_rp, *cur_rp, tmp_e);

					}
				}

				++ cur_r;
				++ cur_c;
				++ offset_r;
				++ offset_c;
				++ r;

			}
			return r;
		}
	};



} // end namespace LinBox


#endif //__LINBOX_matrix_rank_H
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
