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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file   linbox/matrix/abnormal-matrix.inl
 * @ingroup matrix
 * @brief
 */

#ifndef __LINBOX_ABNORMAL_MATRIX_INL
#define __LINBOX_ABNORMAL_MATRIX_INL

#include <stdlib.h>
#include <fstream>

#include "linbox/matrix/abnormal-matrix.h"

namespace LinBox
{

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::AbnormalMatrix() :
        rowDim_(0), colDim_(0) {}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::AbnormalMatrix(const Field& F, Index r, Index c)
{
	init(F,r,c);
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
void AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::init(const Field& F, Index r, Index c)
{
	shape(r,c);
	abnormalHelper_.init(F);
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
void AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::shape(Index r, Index c)
{
	rowDim_=r;
	colDim_=c;
	if (ROW_MAJOR) {
		elts_.shape(r,c);
	} else {
		elts_.shape(c,r);
	}
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
void AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::setEntry(Index i, Index j, const Element& e)
{
	refEntry(i,j)=e;
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
typename AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::Abnormal&
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::refEntry(Index i, Index j)
{
	if (ROW_MAJOR) {
		return *elts_.refElt(i,j);
	} else {
		return *elts_.refElt(j,i);
	}
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
typename AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::Element&
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::getElement(Index i, Index j, Element& e)
{
	return e=abnormalHelper_.normalize(refEntry(i,j));
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
typename AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::AbnormalSubmatrix
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::
submat(Index startRow, Index startCol, Index endRow, Index endCol)
{
	return AbnormalSubmatrix(this,startRow,startCol,endRow,endCol);
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
typename AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::AbnormalSubmatrix
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::
wholeSlice()
{
	return AbnormalSubmatrix(this,0,0,rowDim_-1,colDim_-1);
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
typename AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::AbnormalSubmatrix
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::
rowSlice(Index i)
{
	return AbnormalSubmatrix(this,i,0,i,colDim_-1);
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
typename AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::AbnormalSubmatrix
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::
colSlice(Index j)
{
	return AbnormalSubmatrix(this,0,j,rowDim_-1,j);
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
typename AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::Iterator
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::
AbnormalSubmatrix::
begin()
{
	return Iterator(this,startRow_,startCol_);
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
typename AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::Iterator
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::
AbnormalSubmatrix::
end()
{
        return ++Iterator(this,endRow_,endCol_);
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
template <class Matrix>
void AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::
AbnormalSubmatrix::normalize(Matrix& out)
{
	for (Iterator it=begin();it!=end();++it) {
		out.setEntry(it.row_,it.col_,mat_->abnormalHelper_.normalize(*it));
	}
}

template <class Field, bool ROW_MAJOR,size_t ALIGNMENT>
template <class Matrix>
void
AbnormalMatrix<Field,ROW_MAJOR,ALIGNMENT>::
AbnormalSubmatrix::saxpyin(const Element& e, Matrix& mat)
{
        Iterator it=begin();
        Iterator itEnd=end();
        while (it != itEnd) {
                Index row=it.row_-startRow_;
                Index col=it.col_-startCol_;
		mat_->abnormalHelper_.mulacc(*it,e,mat.getEntry(row,col));
                ++it;
	}
}


}

#endif // __LINBOX_ABNORMAL_MATRIX_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
