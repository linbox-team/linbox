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

/*! @file   examples/map-sparse.inl
 * @ingroup examples
 * @brief
 */

#ifndef __LINBOX_MAP_SPARSE_INL
#define __LINBOX_MAP_SPARSE_INL

#include <stdlib.h>
#include <fstream>

namespace LinBox
{

template<class Field_>
MapSparse<Field_>::MapSparse() : numRows_(0), numCols_(0), nnz_(0) {}

template<class Field_>
MapSparse<Field_>::MapSparse(const Field& F, Index r, Index c):
	numRows_(r), numCols_(c), nnz_(0), MD_(F) {F.init(zero_,0);}

template<class Field_>
MapSparse<Field_>::MapSparse(const MapSparse& M):
	numRows_(M.numRows_), numCols_(M.numCols_),
	rowMap_(M.rowMap_), colMap_(M.colMap_),
	MD_(M.field()), nnz_(M.nnz_), zero_(M.zero_) {}

template<class Field_>
MapSparse<Field_>& MapSparse<Field_>::operator=(const MapSparse<Field_>& rhs)
{
	if (rhs==this) return;
	MD_.init(rhs.MD_);
	numCols_=rhs.numCols_;
	numRows_=rhs.numRows_;
	rowMap_=rhs.rowMap_;
	colMap_=rhs.colMap_;
	nnz_=rhs.nnz_;
	zero_=rhs.zero_;

	return *this;
}

template<class Field_>
MapSparse<Field_>::~MapSparse() {}

template<class Field_>
const Field_& MapSparse<Field_>::field() const
{
	return MD_.field();
}

template<class Field_>
void MapSparse<Field_>::setEntry(Index i, Index j, const Element& e)
{
	VectorIt it=rowMap_[i].find(j);
	if (it != rowMap_[i].end()) {
		--nnz_;
		rowMap_[i].erase(it);
	}

	if (!field().isZero(e)) {
		++nnz_;
		rowMap_[i][j]=e;
		colMap_[j][i]=e;
	}
}

template<class Field_> const typename MapSparse<Field_>::Element&
MapSparse<Field_>::getEntry(Index i, Index j) const
{
	MapConstIt row=rowMap_.find(i);
	if (row != rowMap_.end()) {
		VectorConstIt entry=(row->second).find(j);
		if (entry != (row->second.end())) {
			return (entry->second);
		}
	}
	return zero_;
}

//forall r: A_{i,r}<-A_{i,r}+k*A_{j,r}
template<class Field_>
void MapSparse<Field_>::addRow(const Element& k, Index i, Index j)
{
	MapIt rowJ=rowMap_.find(j);
	if (rowJ != rowMap_.end()) {
		for (VectorIt it=rowJ->second.begin();it!=rowJ->second.end();++it) {
			Index col=it->first;
			if (!(field().isZero(getEntry(i,col)))) {
				--nnz_;
			}
			field().axpyin(rowMap_[i][col],k,it->second);
			field().axpyin(colMap_[col][i],k,it->second);
			if (!(field().isZero(getEntry(i,col)))) {
				++nnz_;
			}
		}
	}
}

//forall r: A_{r,i}<-A_{r,i}+k*A_{r,j}
template<class Field_>
void MapSparse<Field_>::addCol(const Element& k, Index i, Index j)
{
	MapIt colJ=colMap_.find(j);
	if (colJ != colMap_.end()) {
		for (VectorIt it=colJ->second.begin();it!=colJ->second.end();++it) {
			Index row=it->first;
			if (!(field().isZero(getEntry(row,i)))) {
				--nnz_;
			} else {
				field().init(rowMap_[row][i],0);
				field().init(colMap_[i][row],0);
			}
			field().axpyin(rowMap_[row][i],k,it->second);
			field().axpyin(colMap_[i][row],k,it->second);
			if (!(field().isZero(getEntry(row,i)))) {
				++nnz_;
			}
		}
	}
}

//forall r: A_{i,r}<-k*A_{i,r}
template<class Field_>
void MapSparse<Field_>::timesRow(const Element& k, Index i)
{
	linbox_check(!(field().isZero(k)));

	MapIt row=rowMap_.find(i);
	if (row != rowMap_.end()) {
		for (VectorIt it=row->second.begin();it!=row->second.end();++it) {
			Index col=it->first;
			field().mulin(rowMap_[i][col],k);
			field().mulin(colMap_[col][i],k);
		}
	}
}

//forall r: A_{r,j}<-k*A_{r,j}
template<class Field_>
void MapSparse<Field_>::timesCol(const Element& k, Index j)
{
	linbox_check(!(field().isZero(k)));

	MapIt col=colMap_.find(j);
	if (col != colMap_.end()) {
		for (VectorIt it=col->second.begin();it!=col->second.end();++it) {
			Index row=it->first;
			field().mulin(rowMap_[row][j],k);
			field().mulin(colMap_[j][row],k);
		}
	}
}

template<class Field_>
void MapSparse<Field_>::swapRows(Index i, Index j)
{
	VectorType oldRowI=rowMap_[i];
	VectorType oldRowJ=rowMap_[j];

	rowMap_[i].clear();
	rowMap_[j].clear();

	for (VectorIt it=oldRowI.begin();it!=oldRowI.end();++it) {
		Index col=it->first;
		colMap_[col].erase(colMap_[col].find(i));
	}

	for (VectorIt it=oldRowJ.begin();it!=oldRowJ.end();++it) {
		Index col=it->first;
		colMap_[col].erase(colMap_[col].find(j));
	}

	for (VectorIt it=oldRowI.begin();it!=oldRowI.end();++it) {
		Index col=it->first;
		rowMap_[j][col]=oldRowI[col];
		colMap_[col][j]=oldRowI[col];
	}

	for (VectorIt it=oldRowJ.begin();it!=oldRowJ.end();++it) {
		Index col=it->first;
		rowMap_[i][col]=oldRowJ[col];
		colMap_[col][i]=oldRowJ[col];
	}
}

template<class Field_>
void MapSparse<Field_>::swapCols(Index i, Index j)
{
	VectorType oldColI=colMap_[i];
	VectorType oldColJ=colMap_[j];

	colMap_[i].clear();
	colMap_[j].clear();

	for (VectorIt it=oldColI.begin();it!=oldColI.end();++it) {
		Index row=it->first;
		rowMap_[row].erase(rowMap_[row].find(i));
	}

	for (VectorIt it=oldColJ.begin();it!=oldColJ.end();++it) {
		Index row=it->first;
		rowMap_[row].erase(rowMap_[row].find(j));
	}

	for (VectorIt it=oldColI.begin();it!=oldColI.end();++it) {
		Index row=it->first;
		colMap_[j][row]=oldColI[row];
		rowMap_[row][j]=oldColI[row];
	}

	for (VectorIt it=oldColJ.begin();it!=oldColJ.end();++it) {
		Index row=it->first;
		colMap_[i][row]=oldColJ[row];
		rowMap_[row][i]=oldColJ[row];
	}
}

template<class Field_>
typename MapSparse<Field_>::Index MapSparse<Field_>::nnz() const
{
	return nnz_;
}

template<class Field_>
void MapSparse<Field_>::print(std::ostream& out) const
{
	for (Index i=0;i<numRows_;++i) {
		for (Index j=0;j<numCols_;++j) {
			field().write(out,getEntry(i,j));
			if (j != (numCols_-1)) {
				out << " ";
			}
		}
		out << std::endl;
	}
}

template<class Field_>
void MapSparse<Field_>::write(std::ostream& out) const
{
	out << "%%MatrixMarket matrix coordinate integer general" << std::endl;
	out << "% written from a LinBox MapSparse" << std::endl;
        out << numRows_ << " " << numCols_ << " " << nnz_ << std::endl;
        //for (Index i = 0; i < numRows_; ++i)
        for (MapConstIt p = rowMap_.begin(); p != rowMap_.end(); ++p)
                for (VectorConstIt rp = p->second.begin(); rp != p->second.end(); ++rp)
                        field().write(out << 1+p->first << " " << 1+rp->first << " ", rp->second) << std::endl;
}

}

#endif // __LINBOX_MAP_SPARSE_INL

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
