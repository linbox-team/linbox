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

template<class Field_>
void TriplesBBOMP<Field_>::computeRowBlocks()
{
	typedef std::map<Index,Index> IntervalSet;
	typedef typename IntervalSet::iterator IntervalIterator;
	typedef std::pair<Index,Index> Interval;
	IntervalSet intervals;

	for (Index k=0;k<superBlocks_.size();++k) {
		Index start=superBlocks_[k].getStartRow();
		Index end=superBlocks_[k].getEndRow();

		//invariant: no existing interval is dominated by any other

		//if the new interval is not dominated by any other, insert it
		IntervalIterator it=intervals.lower_bound(start);
		bool insert=false;
		if (it == intervals.end()) {
			insert=true;
		} else if (it->first == start) {
			if (it->second < end) {
				insert=true;
			}
		} else {
			if (it != intervals.begin()) {
				--it;
			}
			if (it->first > start) {
				insert=true;
			} else if (it->second < start) {
				insert=true;
			}
		}
		if (insert) {
			intervals[start]=end;
			it=intervals.upper_bound(start);

			//remove all intervals dominated by the new one
			while ((it!=intervals.end()) && (it->first <= end)) {
				++it;
			}
			intervals.erase(intervals.upper_bound(start),it);
		}
	}

	typedef std::map<Index,std::vector<Block> > RowBlockMap;

	RowBlockMap blocksByRow;
	for (IntervalIterator it=intervals.begin();it!=intervals.end();++it) {
		blocksByRow[it->first];
	}

	for (Index k=0;k<superBlocks_.size();++k) {
		Index start=superBlocks_[k].getStartRow();

		IntervalIterator it=intervals.lower_bound(start);
		if (it != intervals.begin()) {
			--it;
		}
		blocksByRow[it->first].push_back(superBlocks_[k]);
	}

	rowBlocks_.clear();
	for (typename RowBlockMap::iterator it=blocksByRow.begin();
	     it!=blocksByRow.end();
	     ++it) {
		rowBlocks_.push_back(it->second);
	}
}

template<class Field_> TriplesBBOMP<Field_>::TriplesBBOMP() : sortType_(TRIPLES_UNSORTED) {}
template<class Field_> TriplesBBOMP<Field_>::~TriplesBBOMP() {}

template<class Field_> TriplesBBOMP<Field_>::
TriplesBBOMP(const Field_& F, istream& in) : MD_(F), data_(),rows_(0),cols_(0)
{
	read(in);
}

template<class Field_>
istream& TriplesBBOMP<Field_>::read(istream& in){
	Index r, c;
	typename Field::Element v; field().init(v);
	MatrixStream<Field> ms(field(), in);
	ms.getDimensions(r, c);
	shape(field(), r, c);
	while (ms.nextTriple(r, c, v)) setEntry(r, c, v);
	return in;
}

template<class Field_>
ostream& TriplesBBOMP<Field_>::write(ostream& out){
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

template<class Field_>
TriplesBBOMP<Field_>& TriplesBBOMP<Field_>::shape(const Field& F, Index r, Index c)
{ MD_=F; data_.clear(); rows_ = r; cols_ = c; sortType_ = TRIPLES_UNSORTED; return *this; }

template<class Field_> TriplesBBOMP<Field_>::
TriplesBBOMP(const Field& F, Index r, Index c)
        : MD_(F), data_(), rows_(r), cols_(c),
          sortType_(TRIPLES_UNSORTED) {}

template<class Field_>
TriplesBBOMP<Field_>::TriplesBBOMP(const TriplesBBOMP<Field_> & B)
        : MD_(B.MD_), data_ ( B.data_ ),
          rows_ ( B.rows_ ), cols_ ( B.cols_ ),
          sortType_ ( B.sortType_ ),
          rowBlocks_(B.rowBlocks_),colBlocks_(B.colBlocks_),
          superBlocks_(B.superBlocks)
{}

template<class Field_>
TriplesBBOMP<Field_> & TriplesBBOMP<Field_>::operator=(const TriplesBBOMP<Field_> & rhs)
{
	if (rhs == this)
		return ;
        MD_.init(rhs.field);
	data_ = rhs.data_;
	rows_ = rhs.rows_;
	cols_ = rhs.cols_;
	sortType_ = rhs.sortType_;
        rowBlocks_=rhs.rowBlocks_;
        colBlocks_=rhs.colBlocks_;
        superBlocks_=rhs.superBlocks_;
	return *this;
}

template<class Field_>
void TriplesBBOMP<Field_>::fromBlock(TriplesBBOMP<Field_>::Coord& coord)
{
	Llong temp,final;
	Index localRow=0,localCol=0;
#ifdef USING_64_BIT_SIZE_T
	if (coord.blockIxArr[1]!=0) {
		final=coord.blockIxArr[1];
		temp=(final^(final>>1))&0x2222222222222222;
		final^=temp^(temp<<1);
		temp=(final^(final>>2))&0x0C0C0C0C0C0C0C0C;
		final^=temp^(temp<<2);
		temp=(final^(final>>4))&0x00F000F000F000F0;
		final^=temp^(temp<<4);
		temp=(final^(final>>8))&0x0000FF000000FF00;
		final^=temp^(temp<<8);
		temp=(final^(final>>16))&0x00000000FFFF0000;
		final^=temp^(temp<<16);
		localRow=(Index)(final&(~0xFFFFFFFF));
		localCol=(Index)((final<<32)&(~0xFFFFFFFF));
	}
#endif
	final=coord.blockIx;
	temp=(final^(final>>1))&0x2222222222222222;
	final^=temp^(temp<<1);
	temp=(final^(final>>2))&0x0C0C0C0C0C0C0C0C;
	final^=temp^(temp<<2);
	temp=(final^(final>>4))&0x00F000F000F000F0;
	final^=temp^(temp<<4);
	temp=(final^(final>>8))&0x0000FF000000FF00;
	final^=temp^(temp<<8);
	temp=(final^(final>>16))&0x00000000FFFF0000;
	final^=temp^(temp<<16);
	localRow|=(Index)((final>>32)&0xFFFFFFFF);
	localCol|=(Index)(final&0xFFFFFFFF);
	coord.rowCol[0]=localRow;
	coord.rowCol[1]=localCol;
}


template<class Field_>
void TriplesBBOMP<Field_>::toBlock(TriplesBBOMP<Field_>::Coord& coord)
{
	Llong temp,final;
	Index localRow=coord.rowCol[0],localCol=coord.rowCol[1];
#ifdef USING_64_BIT_SIZE_T
	Index highOrderLocalRow=localRow&(~0xFFFFFFFF);
	Index highOrderLocalCol=localCol&(~0xFFFFFFFF);
	Llong highOrderFinal;
	coord.blockIxArr[1]=0;
	if ((highOrderLocalRow != 0) || (highOrderLocalCol != 0)) {
		highOrderFinal=((Llong)localRow)|(((Llong)localCol)>>32);
		temp=(highOrderFinal^(highOrderFinal>>16))&0x00000000FFFF0000;
		highOrderFinal^=temp^(temp<<16);
		temp=(highOrderFinal^(highOrderFinal>>8))&0x0000FF000000FF00;
		highOrderFinal^=temp^(temp<<8);
		temp=(highOrderFinal^(highOrderFinal>>4))&0x00F000F000F000F0;
		highOrderFinal^=temp^(temp<<4);
		temp=(highOrderFinal^(highOrderFinal>>2))&0x0C0C0C0C0C0C0C0C;
		highOrderFinal^=temp^(temp<<2);
		temp=(highOrderFinal^(highOrderFinal>>1))&0x2222222222222222;
		highOrderFinal^=temp^(temp<<1);
		coord.blockIxArr[1]=highOrderFinal;
	}
	localCol&=0xFFFFFFFF;
	localRow&=0xFFFFFFFF;
#endif
	final=(((Llong)localRow)<<32)|((Llong)localCol);
	temp=(final^(final>>16))&0x00000000FFFF0000;
	final^=temp^(temp<<16);
	temp=(final^(final>>8))&0x0000FF000000FF00;
	final^=temp^(temp<<8);
	temp=(final^(final>>4))&0x00F000F000F000F0;
	final^=temp^(temp<<4);
	temp=(final^(final>>2))&0x0C0C0C0C0C0C0C0C;
	final^=temp^(temp<<2);
	temp=(final^(final>>1))&0x2222222222222222;
	final^=temp^(temp<<1);
	coord.blockIx=final;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field_>::apply(OutVector & y, const InVector & x) const
{
	linbox_check( coldim() == x.size() );
	linbox_check( rowdim() == y.size() );


        std::vector<FieldAXPY<Field> > yTemp(y.size(),FieldAXPY<Field>(field()));

#pragma omp parallel
        {
#pragma omp for schedule (static,1)
                for (Index i=0;i<rowBlocks_.size();++i) {
                        for (Index j=0;j<rowBlocks_[i].size();++j) {
                                Index start=rowBlocks_[i][j].start_;
                                Index end=rowBlocks_[i][j].end_;
                                for (Index k=start;k<end;++k) {
                                        const Triple& t = data_[k];
                                        yTemp[t.getRow()].mulacc(t.getElt(),x[t.getCol()]);
                                }
                        }
                }
#pragma omp for schedule (static,100)
                for (Index i = 0; i < y.size(); ++i) {
                        yTemp[i].get(y[i]);
                }
        }
        return y;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field_>::seqApply(OutVector & y, const InVector & x) const
{

	linbox_check( coldim() == x.size() );
	linbox_check( rowdim() == y.size() );

	for (Index i = 0; i < y.size(); ++i) field().init(y[i], field().zero);
	for (Index k = 0; k < data_.size(); ++k) {
		const Triple& t = data_[k];
		field().axpyin(y[t.getRow()], t.getElt(), x[t.getCol()]);
	}
	return y;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field_>::applyTranspose(OutVector & y, const InVector & x) const
{
	linbox_check( coldim() == y.size() );
	linbox_check( rowdim() == x.size() );
	for (Index i = 0; i < y.size(); ++i) field().init(y[i], field().zero);
	for (Index k = 0; k < data_.size(); ++k) {
		const Triple& t = data_[k];
		field().axpyin(y[t.getCol()], t.getElt(), x[t.getRow()]);
	}
	return y;
}

template<class Field_>
size_t TriplesBBOMP<Field_>::rowdim() const { return rows_; }

template<class Field_>
size_t TriplesBBOMP<Field_>::coldim() const { return cols_; }

template<class Field_>
const Field_& TriplesBBOMP<Field_>::
field() const { return MD_.field();}

template<class Field_>
size_t TriplesBBOMP<Field_>::size() const { return data_.size(); }

template<class Field_>
void TriplesBBOMP<Field_>::finalize()
{
        sortBlock();
}


template<class Field_>
void TriplesBBOMP<Field_>::mergeBlocks(Index maxSize,
				       std::vector<int>& mergeCounts,
				       std::vector<bool>& mergeFlags)
{
        int i=0;
        while (mergeCounts[i]==4) {
                mergeCounts[i]=0;
                ++mergeCounts[i+1];

		bool shouldMerge=!(mergeFlags[i]);

                //compute nnz of tentative new block
                int numSuperBlocks=superBlocks_.size();
                int newBlockNNZ=0;
                newBlockNNZ+=superBlocks_[numSuperBlocks-1].nnz();
                newBlockNNZ+=superBlocks_[numSuperBlocks-2].nnz();
                newBlockNNZ+=superBlocks_[numSuperBlocks-3].nnz();
                newBlockNNZ+=superBlocks_[numSuperBlocks-4].nnz();

		shouldMerge = shouldMerge && newBlockNNZ<MAX_BLOCK_NNZ;

		superBlocks_[numSuperBlocks-4].fromBlock();
		int newBlockSize = 4*(superBlocks_[numSuperBlocks-4].blockSize());
		superBlocks_[numSuperBlocks-4].toBlock();

		shouldMerge = shouldMerge && (newBlockSize <= maxSize);

		if (shouldMerge) {
			int dataEnd=superBlocks_[numSuperBlocks-1].end_;
			int blockEnd=superBlocks_[numSuperBlocks-1].blockEnd_.blockIx;
                        superBlocks_[numSuperBlocks-4].expand(dataEnd,blockEnd);
                        superBlocks_.pop_back();
                        superBlocks_.pop_back();
                        superBlocks_.pop_back();
                } else {
                        mergeFlags[i+1]=true;
                }
                mergeFlags[i]=false;
		++i;
        }
}

template<class Field_>
void TriplesBBOMP<Field_>::sortBlock()
{
        if ((sortType_ & TRIPLES_SORTED) != 0) {return; }
        for (Index k=0; k<data_.size();++k) {
                data_[k].toBlock();
        }

        std::sort(data_.begin(),data_.end(),Triple::compareBlockTriples);


        superBlocks_.clear();
        Index startBlockIx=0,startDataIx=0;

        std::vector<int> mergeCounts(64,0);
        std::vector<bool> mergeFlags(64,false);

	int maxBlockSize = ((rowdim()*coldim())>>BLOCK_SIZE_BITS)/(EXPECTED_THREADS*EXPECTED_THREADS);
	maxBlockSize = (maxBlockSize==0)?1:maxBlockSize;

        for (Index k=0;k<data_.size();++k) {
                Index blockIx=(data_[k].blockIx()&BLOCK_MASK)>>BLOCK_SIZE_BITS;
                //make new blocks:
                while (startBlockIx < blockIx) {
                        superBlocks_.push_back(Block(startDataIx,k,startBlockIx));
                        startDataIx=k;
                        ++startBlockIx;
                        ++(mergeCounts[0]);
                        mergeBlocks(maxBlockSize,mergeCounts,mergeFlags);
                }
        }

        superBlocks_.push_back(Block(startDataIx,data_.size(),startBlockIx));
        ++(mergeCounts[0]);
        mergeBlocks(maxBlockSize,mergeCounts,mergeFlags);

        for (Index k=0;k<data_.size();++k) {
                data_[k].fromBlock();
        }

	for (Index i=0;i<superBlocks_.size();++i) {
		superBlocks_[i].fromBlock();
	}
	computeRowBlocks();

        sortType_=TRIPLES_SORTED;
}

template<class Field_>
void TriplesBBOMP<Field_>::setEntry(Index i, Index j, const typename Field::Element & e)
{
	sortType_ = TRIPLES_UNSORTED;
	data_.push_back(Triple(i, j, e));
}

template<class Field_>
typename Field_::Element& TriplesBBOMP<Field_>::
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
