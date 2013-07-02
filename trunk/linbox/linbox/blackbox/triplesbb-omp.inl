/* linbox/blackbox/triplesbb-omp.inl
 * Copyright (c) Linbox
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========

 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * with mods by bds
 */

/** @file blackbox/triplesbb-omp.h
 * @ingroup blackbox
 * @brief NO DOC
 */

#ifndef __LINBOX_triplesbb_omp_INL
#define __LINBOX_triplesbb_omp_INL

#include <algorithm>
#include <iostream>
using std::istream;
using std::ostream;
#include <omp.h>
#include "linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/field/hom.h"
#include "linbox/util/matrix-stream.h"

#include <vector>

namespace LinBox
{

template<class Field> TriplesBBOMP<Field>::TriplesBBOMP() : sortType_(TRIPLES_UNSORTED) {}
template<class Field> TriplesBBOMP<Field>::~TriplesBBOMP() {}

template<class Field> TriplesBBOMP<Field>::TriplesBBOMP(const Field& F, istream& in) : field_(&F) {
	read(in);
}

template<class Field>
istream& TriplesBBOMP<Field>::read(istream& in){
	Index r, c;
	typename Field::Element v; field().init(v);
	MatrixStream<Field> ms(field(), in);
	ms.getDimensions(r, c);
	shape(field(), r, c);
	while (ms.nextTriple(r, c, v)) setEntry(r, c, v);
	return in;
}

template<class Field>
ostream& TriplesBBOMP<Field>::write(ostream& out){
	Index r, c;
	out << "%%MatrixMarket matrix coordinate integer general" << std::endl;
	out << "% written from a LinBox TriplesBBOMP" << std::endl;
	out << rowdim() <<" " << coldim() << " " << size() << std::endl;
	for (Index k = 0; k < size(); ++k) {
		Triple t = data_[k];
		field().write(out << t.row()+1 << " " << t.col()+1 << " ", t.elt) << std::endl;
	}
	return out;
}

template<class Field>
TriplesBBOMP<Field>& TriplesBBOMP<Field>::shape(const Field& F, Index r, Index c)
{ field_ = &F; data_.clear(); rows_ = r; cols_ = c; sortType_ = TRIPLES_UNSORTED; return *this; }

template<class Field> TriplesBBOMP<Field>::TriplesBBOMP(const Field& F, Index r, Index c)
: field_(&F), data_(), rows_(r), cols_(c), sortType_(TRIPLES_UNSORTED) {}

template<class Field>
TriplesBBOMP<Field>::TriplesBBOMP(const TriplesBBOMP<Field> & B)
: field_ ( B.field_ ), data_ ( B.data_ ), rows_ ( B.rows_ ), cols_ ( B.cols_ ),
  sortType_ ( B.sortType_ ),
  colBlocks_ (B.colBlocks__), rowBlocks_ (B.rowBlocks_)
{}

template<class Field>
TriplesBBOMP<Field> & TriplesBBOMP<Field>::operator=(const TriplesBBOMP<Field> & rhs)
{
	if (rhs == this)
		return ;
	field_ = rhs.field_;
	data_ = rhs.data_;
	rows_ = rhs.rows_;
	cols_ = rhs.cols_;
	sortType_ = rhs.sortType_;
        colBlocks_ = rhs.colBlocks_;
        rowBlocks_ = rhs.rowBlocks_;
	return *this;
}

template<class Field>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field>::apply(OutVector & y, const InVector & x) const
{

	linbox_check( rowdim() == y.size() );
	linbox_check( coldim() == x.size() );

#pragma omp parallel
        {
#pragma omp for schedule (static,100)
                for (Index i = 0; i < y.size(); ++i) {
                        field().init(y[i], field().zero);
                }
#pragma omp for schedule (static,1)
                for (Index i=0;i<rowBlocks_.size();++i) {
                        for (Index j=0;j<rowBlocks_[i].size();++j) {
                                Index start=rowBlocks_[i][j].start_;
                                Index end=rowBlocks_[i][j].end_;
                                for (Index k=start;k<end;++k) {
                                        const Triple& t = data_[k];
                                        field().axpyin(y[t.getRow()], t.elt, x[t.getCol()]);
                                }
                        }
                }
        }
        return y;
}

template<class Field>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field>::seqApply(OutVector & y, const InVector & x) const
{

	linbox_check( rowdim() == y.size() );
	linbox_check( coldim() == x.size() );
	for (Index i = 0; i < y.size(); ++i) field().init(y[i], field().zero);
	for (Index k = 0; k < data_.size(); ++k) {
		const Triple& t = data_[k];
		field().axpyin(y[t.getRow()], t.elt, x[t.getCol()]);
	}
	return y;
}

template<class Field>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field>::applyTranspose(OutVector & y, const InVector & x) const
{
	linbox_check( coldim() == y.size() );
	linbox_check( rowdim() == x.size() );
	for (Index i = 0; i < y.size(); ++i) field().init(y[i], field().zero);
	for (Index k = 0; k < data_.size(); ++k) {
		const Triple& t = data_[k];
		field().axpyin(y[t.col()], t.elt, x[t.row()]);
	}
	return y;
}

template<class Field>
size_t TriplesBBOMP<Field>::rowdim() const { return rows_; }

template<class Field>
size_t TriplesBBOMP<Field>::coldim() const { return cols_; }

template<class Field>
const Field & TriplesBBOMP<Field>::field() const { return *field_; }

template<class Field>
size_t TriplesBBOMP<Field>::size() const { return data_.size(); }

template<class Field>
void TriplesBBOMP<Field>::sortRow()
{
        if ((sortType_ & TRIPLES_ROW_MAJOR_SORT) != 0) {
                return;
        }
        sortBlock();

        typedef std::map<Index,std::vector<Block> > BlockMap;
        typedef typename BlockMap::iterator BlockMapIt;

        BlockMap blocks;

        Index blockSize=1024;
        Index blockMask=~(blockSize-1);
        Index startRow=0;
        Index startIx=0;
        for (Index k=0;k<data_.size();++k) {
                if (data_[k].getRow()>startRow+blockSize) {
                        BlockMapIt vec=blocks.find(startRow);
                        if (vec==blocks.end()) {
                                blocks.insert(make_pair(startRow,std::vector<Block>(1,Block(startIx,k))));
                        } else {
                                vec->second.push_back(Block(startIx,k));
                        }
                        startRow=data_[k].getRow()&blockMask;
                        startIx=k;
                }
        }
        BlockMapIt vec=blocks.find(startRow);
        if (vec==blocks.end()) {
                blocks.insert(make_pair(startRow,std::vector<Block>(1,Block(startIx,data_.size()))));
        } else {
                vec->second.push_back(Block(startIx,data_.size()));
        }

        rowBlocks_.clear();
        for (BlockMapIt it=blocks.begin();it!=blocks.end();++it) {
                rowBlocks_.push_back(it->second);
        }

        sortType_ |= TRIPLES_ROW_MAJOR_SORT;
}

template<class Field>
void TriplesBBOMP<Field>::sortColumn()
{
        if (sortType_ & TRIPLES_COL_MAJOR_SORT != 0) {
                return;
        }

        sortBlock();

        typedef std::map<Index,std::vector<Block> > BlockMap;
        typedef typename BlockMap::iterator BlockMapIt;

        BlockMap blocks;

        Index blockSize=1024;
        Index blockMask=~(blockSize-1);
        Index startCol=0;
        Index startIx=0;
        for (Index k=0;k<data_.size();++k) {
                if (data_[k].getCol()>startCol+blockSize) {
                        BlockMapIt vec=blocks.find(startCol);
                        if (vec==blocks.end()) {
                                blocks.insert(make_pair(startCol,std::vector<Block>(1,Block(startIx,k))));
                        } else {
                                vec->second.push_back(Block(startIx,k));
                        }
                        startCol=data_[k].getCol()&blockMask;
                        startIx=k;
                }
        }
        BlockMapIt vec=blocks.find(startCol);
        if (vec==blocks.end()) {
                blocks.insert(make_pair(startCol,std::vector<Block>(1,Block(startIx,data_.size()))));
        } else {
                vec->second.push_back(Block(startIx,data_.size()));
        }

        colBlocks_.clear();
        for (BlockMapIt it=blocks.begin();it!=blocks.end();++it) {
                colBlocks_.push_back(it->second);
        }

        sortType_ |= TRIPLES_COL_MAJOR_SORT;
}

template<class Field>
void TriplesBBOMP<Field>::sortBlock()
{
        if (sortType_ & TRIPLES_BLOCK_SORT != 0) { return; }

        for (Index k=0; k<data_.size();++k) {
                data_[k].toBlock();
        }

        std::sort(data_.begin(),data_.end(),Triple::compareBlockTriples);

        for (Index k=0;k<data_.size();++k) {
                data_[k].fromBlock();
        }

        sortType_=TRIPLES_BLOCK_SORT;
}

template<class Field>
void TriplesBBOMP<Field>::setEntry(Index i, Index j, const typename Field::Element & e)
{
	sortType_ = TRIPLES_UNSORTED;
	data_.push_back(Triple(i, j, e));
}

template<class Field>
typename Field::Element& TriplesBBOMP<Field>::getEntry(typename Field::Element& e, Index i, Index j) const
{
	for (Index k = 0; k < data_.size(); ++k)
		if (data_[k].row == i and data_[k].col == j)
			return e = data_[k].elt;
	return e = field().zero;
}

} // namespace LinBox

#endif // __LINBOX_triplesbb_INL


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
