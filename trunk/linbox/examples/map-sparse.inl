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
#include <sstream>
#include <set>
#include <utility>
#include <linbox/util/matrix-stream.h>

namespace LinBox
{

template<class Field_>
MapSparse<Field_>::MapSparse() : MD_(Field_()),numCols_(0), numRows_(0), nnz_(0) {}

template<class Field_>
MapSparse<Field_>::MapSparse(const Field& F, Index r, Index c):
	MD_(F), numCols_(c), numRows_(r), nnz_(0) {F.init(zero_,0);}

template<class Field_>
void MapSparse<Field_>::init(const Field& F, Index r, Index c) {
        MD_=F;
        F.init(zero_,0);
        shape(r,c);
}

template<class Field_>
void MapSparse<Field_>::shape(Index r, Index c) {
        rowMap_.clear();
        colMap_.clear();
        numRows_=r;
        numCols_=c;
        nnz_=0;
}

template<class Field_>
template<class Vector>
void MapSparse<Field_>::fromVector(const Vector& vec, Index r, Index c) {
        shape(r,c);
        if (numCols_==1) {
                for (Index i=0;i<numRows_;++i) {
                        setEntry(i,0,vec[i]);
                }
        } else {
                for (Index j=0;j<numCols_;++j) {
                        setEntry(0,j,vec[j]);
                }
        }
}

template<class Field_>
MapSparse<Field_>::MapSparse(const MapSparse& M):
        MD_(M.field()),
	rowMap_(M.rowMap_), colMap_(M.colMap_),
        numCols_(M.numCols_), numRows_(M.numRows_),
        nnz_(M.nnz_), zero_(M.zero_) {}

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

// A -> A' = SAS^{-1}, and A' has about nnz nonzero entries.
template<class Field_>
void MapSparse<Field_>::randomSim(Index nz, int seed)
{	typename Field::Element a;
	Index i,j;
	MersenneTwister ri;
	typename Field::RandIter r(field());
	//if (seed != 0) { ri.setSeed(seed); r.setSeed(seed); }
	if (seed != 0)
	{	ri.setSeed(seed);
		// ridiculous constructor only seeding!
		typename Field::RandIter s(field(), seed);
		r = s;
	}
	while (nnz() < nz)
	{	r.nonzerorandom(a);
		i = ri.randomIntRange(0, rowdim()); j = ri.randomIntRange(0, coldim());
		addCol(a, i, j);
		//std::cout << nnz() << std::endl;
		field().negin(a);
		addRow(a, i, j);
	}
}

// A -> A' = UAV, with U and V nonsingular, and A' has about nnz nonzero entries.
template<class Field_>
void MapSparse<Field_>::randomEquiv(Index nz, int seed)
{	typename Field::Element a;
	Index i,j;
	MersenneTwister ri;
	typename Field::RandIter r(field());
	if (seed != 0)
	{	ri.setSeed(seed);
		// ridiculous seeding!
		typename Field::RandIter s(field(), seed);
		r = s;
	}
	bool flip = true;
	while (nnz() < nz)
	{	r.nonzerorandom(a);
		i = ri.randomIntRange(0, rowdim()); j = ri.randomIntRange(0, coldim());
		if (flip) addCol(a, i, j);
		else addRow(a, i, j);
		flip = not flip;
	}
}

template<class Field_>
std::ostream& MapSparse<Field_>::print(std::ostream& out) const
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
        return out;
}

template<class Field_>
std::ostream& MapSparse<Field_>::write(std::ostream& out) const
{
	out << "%%MatrixMarket matrix coordinate integer general" << std::endl;
	out << "% written from a LinBox MapSparse" << std::endl;
        out << numRows_ << " " << numCols_ << " " << nnz_ << std::endl;
        //for (Index i = 0; i < numRows_; ++i)
        for (MapConstIt p = rowMap_.begin(); p != rowMap_.end(); ++p)
                for (VectorConstIt rp = p->second.begin(); rp != p->second.end(); ++rp)
                        field().write(out << 1+p->first << " " << 1+rp->first << " ", rp->second) << std::endl;
        out << std::endl;
        return out;
}

template<class Field_>
std::istream& MapSparse<Field_>::read(std::istream& in) {
        Index r,c;
        Element d;
        field().init(d);
        MatrixStream<Field> ms(field(),in);
        ms.getDimensions(r,c);
        shape(r,c);
        while (ms.nextTriple(r,c,d)) setEntry(r,c,d);
        return in;
}

template<class Field_>
typename MapSparse<Field_>::Index MapSparse<Field_>::rowdim() const
{
        return numRows_;
}

template<class Field_>
typename MapSparse<Field_>::Index MapSparse<Field_>::coldim() const
{
        return numCols_;
}

template<class Field_>
void MapSparse<Field_>::transpose()
{
        rowMap_.swap(colMap_);
        int temp=numCols_;numCols_=numRows_;numRows_=temp;
}

template<class Field_>
template<class Matrix>
void MapSparse<Field_>::copy(Matrix& mat) const
{
        std::stringstream ss;
        write(ss);
        mat.read(ss);
        mat.finalize();
}

template<class Field_>
template<class Matrix>
void MapSparse<Field_>::copyFrom(Matrix& mat)
{
        std::stringstream ss;
        mat.write(ss);
        read(ss);
}

template<class Field_>
template<class Vector>
void MapSparse<Field_>::toVector(Vector& vec) const
{
        if (numCols_ == 1) {
                for (Index i=0;i<numRows_;++i) {
                        vec[i]=getEntry(i,0);
                }
        } else {
                for (Index j=0;j<numCols_;++j) {
                        vec[j]=getEntry(0,j);
                }
        }
}

template<class Field_>
void MapSparse<Field_>::scaleMat(const Element& k)
{
        for (size_t i=0;i<numRows_;++i) {
                timesRow(k,i);
        }
}

template<class Field_>
bool MapSparse<Field_>::areEqual(MapSparse<Field_> rhs) const
{
        if (rhs.numCols_ != numCols_) {
                return false;
        }
        if (rhs.numRows_ != numRows_) {
                return false;
        }
        if (nnz_ != rhs.nnz_) {
                return false;
        }

        for (MapConstIt rowIt=rowMap_.begin();rowIt!=rowMap_.end();++rowIt) {
                for (VectorConstIt eltIt=rowIt->second.begin();eltIt!=rowIt->second.end();++eltIt) {
                        if (eltIt->second != rhs.getEntry(rowIt->first,eltIt->first)) {
                                return false;
                        }
                }
        }
        return true;
}

template<class Field>
void MapSparse<Field>::generateDenseRandMat(MapSparse<Field>& mat, int q)
{
        typedef typename Field::Element Element;

        size_t m=mat.rowdim(),n=mat.coldim();
	Element d;

        for (size_t i=0;i<m;++i) {
                for (size_t j=0;j<n;++j) {
                        mat.field().init(d,randRange(0,q));
                        mat.setEntry(i,j,d);
                }
        }
}

template<class Field>
void MapSparse<Field>::generateRandMat(MapSparse<Field>& mat, int nnz, int q)
{
        typedef typename Field::Element Element;

        size_t m=mat.rowdim(),n=mat.coldim();
	Element d;

        typedef std::pair<size_t,size_t> CoordPair;
        typedef std::set<CoordPair> PairSet;
        PairSet pairs;

	for(int i = 0; i < (int)nnz; ++i) {
                size_t row,col;
                do {
                        row = randRange(0,m);
                        col = randRange(0,n);
                } while (pairs.count(CoordPair(row,col))!=0);

                mat.field().init(d, randRange(1,q));
                mat.setEntry(row,col,d);
                pairs.insert(CoordPair(row,col));
        }
}

template<class Field>
void MapSparse<Field>::generateScaledIdent(MapSparse<Field>& mat, int alpha)
{
        typedef typename Field::Element Element;
        size_t m=mat.rowdim(),n=mat.coldim();
        size_t minDim=(m<n)?m:n;
	Element d;
        mat.field().init(d,alpha);

        for (size_t i=0;i<minDim;++i) {
                mat.setEntry(i,i,d);
        }
}

template<class Field>
int MapSparse<Field>::randRange(int start, int end)
{
        double rval = rand();
        static const double NORMALIZING_CONSTANT = 1.0/(1.0+RAND_MAX);
        double normedRVal = rval*NORMALIZING_CONSTANT;
        double rangeSize = end-start;
        int offset = rangeSize*normedRVal;
        return start+offset;
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
