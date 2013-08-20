/* Copyright (C) 2013 LinBox
 * Written by AJS <stachnik@udel.edu>
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

/*! @file   linbox/matrix/abnormal-matrix.h
 * @ingroup matrix
 * @brief
 */

#ifndef __LINBOX_ABNORMAL_MATRIX
#define __LINBOX_ABNORMAL_MATRIX

#include <stdlib.h>
#include <fstream>
#include "linbox/matrix/abnormal-helpers.h"

namespace LinBox
{

template <class Field, class Matrix>
class AbnormalMatrix {
public:
	typedef size_t Index;
	typedef typename Field::Element Element;

	AbnormalMatrix() {}

	AbnormalMatrix(const Field& F, Matrix& mat) {init(F,mat);}

	void init (const Field& F, Matrix& mat) {
		mat_=&mat;
		cols_=mat_->coldim();
		rows_=mat_->rowdim();
		field_=&F;
		helper_.init(F);
	}

	inline Element& getEntry(Element& e, Index i, Index j) {
		return mat_->getEntry(e,i,j);
	}

	inline void setEntry(Index i, Index j, const Element& e) {
		mat_->setElement(i,j,e);
	}

	template <class Mat2>
	void saxpyin(const Element& e, const Mat2& rhs,
		     Index startRow, Index startCol,
		     Index numRows, Index numCols) {
		for (Index i=0;i<numRows;++i) {
			for (Index j=0;j<numCols;++j) {
				helper_.mulacc(mat_->refEntry(i+startRow,j+startCol),
					       e,
					       rhs.getEntry(i,j));
			}
		}
	}

	void normalize() {
		for (Index i=0;i<rows_;++i) {
			for (Index j=0;j<cols_;++j) {
				helper_.normalize(mat_->refEntry(i,j));
			}
		}
	}

protected:
	Matrix* mat_;

	AbnormalHelper<Field> helper_;

	const Field* field_;

	Index cols_, rows_;
};

}

#endif // __LINBOX_ABNORMAL_MATRIX

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
