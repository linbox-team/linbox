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

template<class Field_> void
TriplesBBOMP<Field_>::nonOverlappingIntervals(BlockListIt startIt,
                                              BlockListIt endIt,
                                              IntervalSet& intervals,
                                              const int rowOrCol)
{
        for (;startIt!=endIt;++startIt) {
                Index start,end;
                if (rowOrCol==CHUNK_BY_ROW) {
                        start=(*startIt).getStartRow();
                        end=(*startIt).getEndRow();
                } else {
                        start=(*startIt).getStartCol();
                        end=(*startIt).getEndCol();
                }

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
}

template<class Field_>
void TriplesBBOMP<Field_>::combineIntervals(BlockListIt startIt,
                                            BlockListIt endIt,
                                            IntervalSet& intervals,
                                            VectorChunks& chunksOut,
                                            const int rowOrCol)
{
	typedef std::map<Index,std::vector<Block> > VectorBlockMap;

	VectorBlockMap vectorMap;
	for (IntervalIterator it=intervals.begin();it!=intervals.end();++it) {
		vectorMap[it->first];
	}

        for (;startIt!=endIt;++startIt) {
                Index start;
                if (rowOrCol==CHUNK_BY_ROW) {
                        start=startIt->getStartRow();
                } else {
                        start=startIt->getStartCol();
                }

		IntervalIterator it=intervals.lower_bound(start);
		if (it != intervals.begin()) {
			--it;
		}
		vectorMap[it->first].push_back(*startIt);
	}

        chunksOut.clear();
        for (typename VectorBlockMap::iterator it=vectorMap.begin();
             it!=vectorMap.end();
             ++it) {
                chunksOut.push_back(it->second);
        }
}

template<class Field_>
void TriplesBBOMP<Field_>::computeVectors(const int rowOrCol)
{
	IntervalSet intervals;
        VectorChunks chunks;

        if (rowOrCol==CHUNK_BY_ROW) {
                rowBlocks_.clear();
        } else {
                colBlocks_.clear();
        }

        std::sort(superBlocks_.begin(),superBlocks_.end(),Block::compareBlockSizes);

        Index k=0;
        Index numSuperBlocks=superBlocks_.size();
        while (k<numSuperBlocks) {
                Coord blockSize=superBlocks_[k].blockSize();
                BlockListIt start=superBlocks_.begin()+k;
                while (superBlocks_[k].blockSize()==blockSize) {
                        ++k;
                        if (k>=numSuperBlocks) {
                                break;
                        }
                }
                BlockListIt end=superBlocks_.begin()+k;

                intervals.clear();
                nonOverlappingIntervals(start,end,intervals,rowOrCol);
                chunks.clear();
                combineIntervals(start,end,intervals,chunks,rowOrCol);
                if (rowOrCol==CHUNK_BY_ROW) {
                        rowBlocks_.push_back(chunks);
                } else {
                        colBlocks_.push_back(chunks);
                }
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
		field().write(out << t.getRow()+1 << " " << t.getCol()+1 << " ", t.getElt()) << std::endl;
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
template<class Mat1, class Mat2> Mat1& TriplesBBOMP<Field_>::
applyLeft(Mat1 &Y, const Mat2 &X) const
{
        Y.zero();

#pragma omp parallel
        {
                Matrix Yc,Xc;

                for (Index chunkSizeIx=0;chunkSizeIx<rowBlocks_.size();++chunkSizeIx) {
#pragma omp for schedule (static,1)
                        for (Index rowChunk=0;rowChunk<rowBlocks_[chunkSizeIx].size();++rowChunk) {
                                for (Index subChunk=0;subChunk<rowBlocks_[chunkSizeIx][rowChunk].size();++subChunk) {
                                        Index start=rowBlocks_[chunkSizeIx][rowChunk][subChunk].start_;
                                        Index end=rowBlocks_[chunkSizeIx][rowChunk][subChunk].end_;
                                        for (Index k=start;k<end;++k) {
                                                const Triple& t = data_[k];
                                                Yc.submatrix(Y,t.getRow(),0,1,Y.coldim());
                                                Xc.submatrix(X,t.getCol(),0,1,X.coldim());
                                                MD_.saxpyin(Yc,t.getElt(),Xc);
                                        }
                                }

                        }
                }

        }
        return Y;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field_>::apply(OutVector & y, const InVector & x) const
{
	linbox_check( coldim() == x.size() );
	linbox_check( rowdim() == y.size() );

        std::vector<FieldAXPY<Field> > yTemp(y.size(),FieldAXPY<Field>(field()));

        for (Index chunkSizeIx=0;chunkSizeIx<rowBlocks_.size();++chunkSizeIx) {
#pragma omp parallel for schedule (static,1)
                for (Index rowChunk=0;rowChunk<rowBlocks_[chunkSizeIx].size();++rowChunk) {
                        for (Index subChunk=0;subChunk<rowBlocks_[chunkSizeIx][rowChunk].size();++subChunk) {
                                Index start=rowBlocks_[chunkSizeIx][rowChunk][subChunk].start_;
                                Index end=rowBlocks_[chunkSizeIx][rowChunk][subChunk].end_;
                                for (Index k=start;k<end;++k) {
                                        const Triple& t = data_[k];
                                        yTemp[t.getRow()].mulacc(t.getElt(),x[t.getCol()]);
                                }
                        }
                        
                }
        }
#pragma omp parallel for schedule (static,100)
        for (Index i = 0; i < y.size(); ++i) {
                yTemp[i].get(y[i]);
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
void TriplesBBOMP<Field_>::splitBlock(TriplesBBOMP<Field_>::Block block)
{
        if (block.nnz()<MAX_BLOCK_NNZ) {
                --block.blockEnd_;
                superBlocks_.push_back(block);
                return;
        }

        Triple blockStart(block.blockStart_);
        Triple blockEnd(block.blockEnd_);
        Coord blockSize=block.blockSize();
        Coord quarterSize=blockSize>>2;
        Triple midOne(blockStart.coord_+quarterSize);
        Triple midTwo(midOne.coord_+quarterSize);
        Triple midThree(midTwo.coord_+quarterSize);

	typedef typename std::vector<Triple>::iterator BlockVecIt;

        BlockVecIt dataBegin=data_.begin()+block.start_;
        BlockVecIt dataEnd=data_.begin()+block.end_;
        Index midOneData,midTwoData,midThreeData;

        midOneData=std::upper_bound(dataBegin,dataEnd,midOne,Triple::compareBlockTriples)-data_.begin();
        midTwoData=std::upper_bound(dataBegin,dataEnd,midTwo,Triple::compareBlockTriples)-data_.begin();
        midThreeData=std::upper_bound(dataBegin,dataEnd,midThree,Triple::compareBlockTriples)-data_.begin();

        Block blockOne(block.start_,midOneData,block.blockStart_,midOne.coord_);
        Block blockTwo(midOneData,midTwoData,midOne.coord_,midTwo.coord_);
        Block blockThree(midTwoData,midThreeData,midTwo.coord_,midThree.coord_);
        Block blockFour(midThreeData,block.end_,midThree.coord_,block.blockEnd_);

        splitBlock(blockOne);
        splitBlock(blockTwo);
        splitBlock(blockThree);
        splitBlock(blockFour);
}

template<class Field_>
void TriplesBBOMP<Field_>::sortBlock()
{
        double start=omp_get_wtime();
        if ((sortType_ & TRIPLES_SORTED) != 0) {return; }
        for (Index k=0; k<data_.size();++k) {
                data_[k].toBlock();
        }

        std::sort(data_.begin(),data_.end(),Triple::compareBlockTriples);

        superBlocks_.clear();

        Index rowBound=roundUpIndex(rowdim());
        Index colBound=roundUpIndex(coldim());
        Index rowColBound=(rowBound<colBound)?colBound:rowBound;
        Coord startBlockIx(0,0),endBlockIx(0,rowColBound);
        coordToBlock(endBlockIx);
        Block initialBlock(0,data_.size(),startBlockIx,endBlockIx);

        splitBlock(initialBlock);

        for (Index k=0;k<data_.size();++k) {
                data_[k].fromBlock();
        }

	for (Index i=0;i<superBlocks_.size();++i) {
		superBlocks_[i].fromBlock();
	}
        computeVectors(CHUNK_BY_ROW);
        computeVectors(CHUNK_BY_COL);

        sortType_=TRIPLES_SORTED;
        //std::cerr << omp_get_wtime()-start << std::endl;
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
