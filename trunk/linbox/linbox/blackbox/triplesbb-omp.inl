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

template<class MatDom> TriplesBBOMP<MatDom>::TriplesBBOMP() : sortType_(TRIPLES_UNSORTED) {}
template<class MatDom> TriplesBBOMP<MatDom>::~TriplesBBOMP() {}

template<class MatDom> TriplesBBOMP<MatDom>::
TriplesBBOMP(const Field& F, istream& in) : field_(&F), MD_(F)
{
	read(in);
}

template<class MatDom>
istream& TriplesBBOMP<MatDom>::read(istream& in){
	Index r, c;
	typename Field::Element v; field().init(v);
	MatrixStream<Field> ms(field(), in);
	ms.getDimensions(r, c);
	shape(field(), r, c);
	while (ms.nextTriple(r, c, v)) setEntry(r, c, v);
	return in;
}

template<class MatDom>
ostream& TriplesBBOMP<MatDom>::write(ostream& out){
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

template<class MatDom>
TriplesBBOMP<MatDom>& TriplesBBOMP<MatDom>::shape(const Field& F, Index r, Index c)
{ field_ = &F; MD_=F; data_.clear(); rows_ = r; cols_ = c; sortType_ = TRIPLES_UNSORTED; return *this; }

template<class MatDom> TriplesBBOMP<MatDom>::
TriplesBBOMP(const Field& F, Index r, Index c)
        : field_(&F),MD_(F), data_(), rows_(r), cols_(c),
          sortType_(TRIPLES_UNSORTED) {}

template<class MatDom>
TriplesBBOMP<MatDom>::TriplesBBOMP(const TriplesBBOMP<MatDom> & B)
        : field_ ( B.field_ ), MD_(B.MD_), data_ ( B.data_ ),
          rows_ ( B.rows_ ), cols_ ( B.cols_ ),
          sortType_ ( B.sortType_ ),
          colBlocks_ (B.colBlocks__), rowBlocks_ (B.rowBlocks_)
{}

template<class MatDom>
TriplesBBOMP<MatDom> & TriplesBBOMP<MatDom>::operator=(const TriplesBBOMP<MatDom> & rhs)
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

template<class MatDom>
typename TriplesBBOMP<MatDom>::Matrix& TriplesBBOMP<MatDom>::
apply_right(typename TriplesBBOMP<MatDom>::Matrix &Y,
           const typename TriplesBBOMP<MatDom>::Matrix &X) const
{
        Y.zero();
#pragma omp parallel
        {
                typename MatrixDomain::Submatrix Yr, Xr; // row submatrices
#pragma omp for schedule (static,1)
                for (Index i=0;i<rowBlocks_.size();++i) {
                        for (Index j=0;j<rowBlocks_[i].size();++j) {
                                Index start=rowBlocks_[i][j].start_;
                                Index end=rowBlocks_[i][j].end_;
                                for (Index k=start;k<end;++k) {
                                        const Triple& t = data_[k];
                                        Yr.submatrix(Y,0,t.getCol(),Y.rowdim(),1);
                                        Xr.submatrix(X,0,t.getRow(),X.rowdim(),1);
                                        MD_.axpyin(Yr, t.elt, Xr);
                                }
                        }
                }
        }
	return Y;
}

template<class MatDom>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<MatDom>::apply(OutVector & y, const InVector & x) const
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

template<class MatDom>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<MatDom>::seqApply(OutVector & y, const InVector & x) const
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

template<class MatDom>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<MatDom>::applyTranspose(OutVector & y, const InVector & x) const
{
	linbox_check( coldim() == y.size() );
	linbox_check( rowdim() == x.size() );
	for (Index i = 0; i < y.size(); ++i) field().init(y[i], field().zero);
	for (Index k = 0; k < data_.size(); ++k) {
		const Triple& t = data_[k];
		field().axpyin(y[t.getCol()], t.elt, x[t.getRow()]);
	}
	return y;
}

template<class MatDom>
size_t TriplesBBOMP<MatDom>::rowdim() const { return rows_; }

template<class MatDom>
size_t TriplesBBOMP<MatDom>::coldim() const { return cols_; }

template<class MatDom>
const typename MatDom::Field& TriplesBBOMP<MatDom>::
field() const { return *field_; }

template<class MatDom> const MatDom& TriplesBBOMP<MatDom>::
domain() const { return MD_; }


template<class MatDom>
size_t TriplesBBOMP<MatDom>::size() const { return data_.size(); }

template<class MatDom>
void TriplesBBOMP<MatDom>::sortRow()
{
        if ((sortType_ & TRIPLES_ROW_MAJOR_SORT) != 0) {
                return;
        }
        sortBlock();

        typedef std::map<Index,std::vector<Block> > BlockMap;
        typedef typename BlockMap::iterator BlockMapIt;

        BlockMap blocks;

        Index blockMask=~(65535);
        Index startRowBlock=0;
        Index startIx=0;
        for (Index k=0;k<data_.size();++k) {
                Index rowBlock=(data_[k].getRow())&blockMask;
                if (rowBlock!=startRowBlock) {
                        BlockMapIt vec=blocks.find(startRowBlock);
                        if (vec==blocks.end()) {
                                blocks.insert(make_pair(startRowBlock,std::vector<Block>(1,Block(startIx,k))));
                        } else {
                                vec->second.push_back(Block(startIx,k));
                        }
                        startRowBlock=rowBlock;
                        startIx=k;
                }
        }
        BlockMapIt vec=blocks.find(startRowBlock);
        if (vec==blocks.end()) {
                blocks.insert(make_pair(startRowBlock,std::vector<Block>(1,Block(startIx,data_.size()))));
        } else {
                vec->second.push_back(Block(startIx,data_.size()));
        }

        rowBlocks_.clear();
        for (BlockMapIt it=blocks.begin();it!=blocks.end();++it) {
                rowBlocks_.push_back(it->second);
        }

        sortType_ |= TRIPLES_ROW_MAJOR_SORT;
}

template<class MatDom>
void TriplesBBOMP<MatDom>::sortColumn()
{
        if (sortType_ & TRIPLES_COL_MAJOR_SORT != 0) {
                return;
        }

        sortBlock();

        typedef std::map<Index,std::vector<Block> > BlockMap;
        typedef typename BlockMap::iterator BlockMapIt;

        BlockMap blocks;

        Index blockMask=~(1023);
        Index startColBlock=0;
        Index startIx=0;
        for (Index k=0;k<data_.size();++k) {
                Index colBlock=data_[k].getCol&blockMask;
                if (colBlock!=startColBlock) {
                        BlockMapIt vec=blocks.find(startColBlock);
                        if (vec==blocks.end()) {
                                blocks.insert(make_pair(startColBlock,std::vector<Block>(1,Block(startIx,k))));
                        } else {
                                vec->second.push_back(Block(startIx,k));
                        }
                        startColBlock=data_[k].getCol()&blockMask;
                        startIx=k;
                }
        }
        BlockMapIt vec=blocks.find(startColBlock);
        if (vec==blocks.end()) {
                blocks.insert(make_pair(startColBlock,std::vector<Block>(1,Block(startIx,data_.size()))));
        } else {
                vec->second.push_back(Block(startIx,data_.size()));
        }

        colBlocks_.clear();
        for (BlockMapIt it=blocks.begin();it!=blocks.end();++it) {
                colBlocks_.push_back(it->second);
        }

        sortType_ |= TRIPLES_COL_MAJOR_SORT;
}

template<class MatDom>
void TriplesBBOMP<MatDom>::sortBlock()
{
        if ((sortType_ & TRIPLES_BLOCK_SORT) != 0) {return; }
        for (Index k=0; k<data_.size();++k) {
                data_[k].toBlock();
        }

        std::sort(data_.begin(),data_.end(),Triple::compareBlockTriples);

        for (Index k=0;k<data_.size();++k) {
                data_[k].fromBlock();
        }

        sortType_=TRIPLES_BLOCK_SORT;
}

template<class MatDom>
void TriplesBBOMP<MatDom>::setEntry(Index i, Index j, const typename Field::Element & e)
{
	sortType_ = TRIPLES_UNSORTED;
	data_.push_back(Triple(i, j, e));
}

template<class MatDom>
typename MatDom::Field::Element& TriplesBBOMP<MatDom>::
getEntry(typename Field::Element& e, Index i, Index j) const
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
