/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
 * Written by Zhendong Wan
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
 * Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */

#ifndef __LINBOX_scompose_H
#define __LINBOX_scompose_H

/*! @file scompose.h
 * Implemenatation of EGV and EGV+ algorithm.
 * Compute the perturbation A+UV and LAR
 */

#include <linbox/util/debug.h>
#include <linbox/blackbox/blackbox-interface.h>
#include <linbox/blackbox/dense.h>
#include <linbox/blackbox/compose.h>

namespace LinBox
{

	class SCompose : public  BlackboxInterface {

	public:

		//- general case, composeSmall is the linbox compose
		template<class Blackbox>

		//For EGV+ algorithm, using LAR.
		static Compose<Blackbox, Blackbox>*& compose (Compose<Blackbox, Blackbox>*& LAR,
							      const Blackbox& L,
							      const Blackbox& A,
							      const Blackbox& R)
		{

			linbox_check (L.coldim() == A.rowdim());

			linbox_check (A.coldim() == R.rowdim());

			std::vector<const Blackbox*> VB(3);

			VB[0] = &L;

			VB[1] = &A;

			VB[2] = &R;

			LAR = new Compose<Blackbox>(VB);

			return LAR;
		}

#if 0
		// specialization for dense matrix case, explicitly compute the LAR by matrix multiplication
		template <class Field>
		static Protected::DenseMatrix<Field>*& compose (Protected::DenseMatrix<Field>*& LAR,
						     const Protected::DenseMatrix<Field>& L,
						     const Protected::DenseMatrix<Field>& A,
						     const Protected::DenseMatrix<Field>& R)
		{

			linbox_check (L.coldim() == A.rowdim());

			linbox_check (A.coldim() == R.rowdim());

			LAR = new Protected::DenseMatrix<Field>(L.field(), L.rowdim(), R.coldim());

			typename Protected::DenseMatrix<Field>::ConstRowIterator crow_p;

			typename Protected::DenseMatrix<Field>::RowIterator row_p;

			std::vector<typename Field::Element> tmp(R.rowdim());

			for (row_p = LAR -> rowBegin(), crow_p = L.rowBegin();
			     row_p != LAR -> rowEnd(); ++ row_p, ++ crow_p) {

				A.applyTranspose(tmp, *crow_p);

				R.applyTranspose(*row_p, tmp);
			}

			return LAR;

		}
#endif

		template <class Field>
		static BlasBlackbox<Field>*& compose (BlasBlackbox<Field>*& LAR,
						     const BlasBlackbox<Field>& L,
						     const BlasBlackbox<Field>& A,
						     const BlasBlackbox<Field>& R)
		{

			linbox_check (L.coldim() == A.rowdim());

			linbox_check (A.coldim() == R.rowdim());

			LAR = new BlasBlackbox<Field>(L.field(), L.rowdim(), R.coldim());

			typename BlasBlackbox<Field>::ConstRowIterator crow_p;

			typename BlasBlackbox<Field>::RowIterator row_p;

			std::vector<typename Field::Element> tmp(R.rowdim());

			for (row_p = LAR -> rowBegin(), crow_p = L.rowBegin();
			     row_p != LAR -> rowEnd(); ++ row_p, ++ crow_p) {

				A.applyTranspose(tmp, *crow_p);

				R.applyTranspose(*row_p, tmp);
			}

			return LAR;

		}

#if 0
		//- Compute A + UV, for EGV algorithm, not be used any more.
		template <class Blackbox>
		static Blackbox*& composeBig (Blackbox*& AUV,
					      const Blackbox& A,
					      const Blackbox& U,
					      const Blackbox& V);



		// @brief This composeBig creates A + UV for EGV algorithm for the Protected::DenseMatrix case.
		template <class Field>
		static Protected::DenseMatrix<Field>*& composeBig (Protected::DenseMatrix<Field>*& AUV,
							const Protected::DenseMatrix<Field>& A,
							const Protected::DenseMatrix<Field>& U,
							const Protected::DenseMatrix<Field>& V) {

			linbox_check (U.rowdim() == A.rowdim());

			linbox_check (A.coldim() == V.coldim());

			AUV = new Protected::DenseMatrix<Field>(A.field(), A.rowdim(), A.coldim());

			typename Protected::DenseMatrix<Field>::ConstRowIterator crow_p;

			typename Protected::DenseMatrix<Field>::RowIterator row_p;

			for (row_p = AUV -> rowBegin(), crow_p = U.rowBegin();
			     row_p != AUV -> rowEnd(); ++ row_p, ++ crow_p) {

				V.applyTranspose(*row_p, *crow_p);

			}

			typename Protected::DenseMatrix<Field>::ConstRawIterator celt_p;
			typename Protected::DenseMatrix<Field>::RawIterator elt_p;

			for( elt_p = AUV -> rawBegin(), celt_p = A.rawBegin(); celt_p !=  A.rawEnd(); ++ elt_p, ++ celt_p)
				A.field().addin(*elt_p,*celt_p);

			return AUV;

		}
#endif



	};
}


#endif //__LINBOX_scompose_H


