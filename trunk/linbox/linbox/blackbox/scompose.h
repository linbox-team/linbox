/* Copyright (C) 2010 LinBox
 * Written by Zhendong Wan
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

#ifndef __LINBOX_scompose_H
#define __LINBOX_scompose_H

/*! @file scompose.h
 * Implemenatation of EGV and EGV+ algorithm.
 * Compute the perturbation A+UV and LAR
 */

#include "linbox/util/debug.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/blackbox/compose.h"

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


		template <class Field>
		static BlasMatrix<Field>*& compose (BlasMatrix<Field>*& LAR,
						     const BlasMatrix<Field>& L,
						     const BlasMatrix<Field>& A,
						     const BlasMatrix<Field>& R)
		{

			linbox_check (L.coldim() == A.rowdim());

			linbox_check (A.coldim() == R.rowdim());

			LAR = new BlasMatrix<Field>(L.field(), L.rowdim(), R.coldim());

			typename BlasMatrix<Field>::ConstRowIterator crow_p;

			typename BlasMatrix<Field>::RowIterator row_p;

			// std::vector<typename Field::Element> tmp(R.rowdim());
			BlasVector<Field> tmp(R.field(),R.rowdim());

			for (row_p = LAR -> rowBegin(), crow_p = L.rowBegin();
			     row_p != LAR -> rowEnd(); ++ row_p, ++ crow_p) {

				A.applyTranspose(tmp, *crow_p);

				// BlasSubvector<BlasVector<Field> > row_p_v(R.field(),*row_p); //! @warning should not be necessary (RowIterator should be a BlasSubvector...).
				// R.applyTranspose(row_p_v, tmp);
				R.applyTranspose(*row_p, tmp);
			}

			return LAR;

		}
	};
}


#endif //__LINBOX_scompose_H



// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

