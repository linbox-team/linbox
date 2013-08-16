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
#include "linbox/matrix/abnormal-matrix.h"

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
	typedef std::map<Index,BlockList> VectorBlockMap;

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
void TriplesBBOMP<Field_>::computeVectors(SizedChunks& sizedChunks,
					  BlockList& superBlocks,
					  const int rowOrCol)
{
	IntervalSet intervals;
        VectorChunks chunks;

        Index k=0;
        Index numSuperBlocks=superBlocks.size();
        while (k<numSuperBlocks) {
                TriplesCoord blockSize=superBlocks[k].blockSize();
                BlockListIt start=superBlocks.begin()+k;
                while (superBlocks[k].blockSize()==blockSize) {
                        ++k;
                        if (k>=numSuperBlocks) {
                                break;
                        }
                }
                BlockListIt end=superBlocks.begin()+k;

                intervals.clear();
                nonOverlappingIntervals(start,end,intervals,rowOrCol);
                chunks.clear();
                combineIntervals(start,end,intervals,chunks,rowOrCol);
		sizedChunks.push_back(chunks);
        }
}

template<class Field_> TriplesBBOMP<Field_>::TriplesBBOMP() {}
template<class Field_> TriplesBBOMP<Field_>::~TriplesBBOMP() {}

template<class Field_> TriplesBBOMP<Field_>::
TriplesBBOMP(const Field_& F, istream& in) : MD_(F)
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
        : MD_(F), rows_(r), cols_(c),
          sortType_(TRIPLES_UNSORTED) {}

template<class Field_>
TriplesBBOMP<Field_>::TriplesBBOMP(const TriplesBBOMP<Field_> & B)
        : MD_(B.MD_), data_ ( B.data_ ),
          rows_ ( B.rows_ ), cols_ ( B.cols_ ),
          sortType_ ( B.sortType_ ),
          rowBlocks_(B.rowBlocks_),colBlocks_(B.colBlocks_)
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
	return *this;
}

template<class Field_>
template<class Mat1, class Mat2> Mat1& TriplesBBOMP<Field_>::
applyLeft(Mat1 &Y, const Mat2 &X) const
{
        typedef AbnormalMatrix<Field_,true,64> AbnormalMat;
        typedef typename AbnormalMat::AbnormalSubmatrix AbnormalSubmat;
        AbnormalMat YTemp(field(),Y.rowdim(),Y.coldim());

#pragma omp parallel
	{
                Matrix Xc;

		Index numBlockSizes=rowBlocks_.size();
		for (Index chunkSizeIx=0;chunkSizeIx<numBlockSizes;++chunkSizeIx) {
			const VectorChunks *rowChunks=&(rowBlocks_[chunkSizeIx]);
			Index numChunks=rowChunks->size();
#pragma omp for schedule (static,1)
			for (Index rowChunk=0;rowChunk<numChunks;++rowChunk) {
				const BlockList *blocks=&((*rowChunks)[rowChunk]);
				Index numBlocks=blocks->size();
				for (Index block=0;block<numBlocks;++block) {
                                        const DataBlock *dataBlock=&((*blocks)[block]);
					for (Index k=0;k<dataBlock->elts_.size();++k) {
                                                const Index row=dataBlock->getRow(k);
                                                const Index col=dataBlock->getCol(k);
                                                AbnormalSubmat Yc=YTemp.rowSlice(row);
                                                Xc.submatrix(X,col,0,1,X.coldim());
                                                Yc.saxpyin(dataBlock->elts_[k],Xc);
                                        }
                                }
                        }
                }
        }
        YTemp.wholeSlice().normalize(Y);
        return Y;
}

template<class Field_>
template<class Mat1, class Mat2> Mat1& TriplesBBOMP<Field_>::
applyRight(Mat1 &Y, const Mat2 &X) const
{
        typedef AbnormalMatrix<Field_,true,64> AbnormalMat;
        typedef typename AbnormalMat::AbnormalSubmatrix AbnormalSubmat;
        AbnormalMat YTemp(field(),Y.rowdim(),Y.coldim());

#pragma omp parallel
	{
                Matrix Xr;

		Index numBlockSizes=colBlocks_.size();
		for (Index chunkSizeIx=0;chunkSizeIx<numBlockSizes;++chunkSizeIx) {
			const VectorChunks *colChunks=&(colBlocks_[chunkSizeIx]);
			Index numChunks=colChunks->size();
#pragma omp for schedule (static,1)
			for (Index colChunk=0;colChunk<numChunks;++colChunk) {
				const BlockList *blocks=&((*colChunks)[colChunk]);
				Index numBlocks=blocks->size();
				for (Index block=0;block<numBlocks;++block) {
                                        const DataBlock *dataBlock=&((*blocks)[block]);
					for (Index k=0;k<dataBlock->elts_.size();++k) {
                                                const Index row=dataBlock->getRow(k);
                                                const Index col=dataBlock->getCol(k);
                                                AbnormalSubmat Yr=YTemp.colSlice(col);
						Xr.submatrix(X,0,row,X.rowdim(),1);
                                                Yr.saxpyin(dataBlock->elts_[k],Xr);
                                        }
                                }
                        }
                }
        }
        YTemp.wholeSlice().normalize(Y);
        return Y;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field_>::apply(OutVector & y, const InVector & x) const
{
	linbox_check( coldim() == x.size() );
	linbox_check( rowdim() == y.size() );

	uint8_t* yTempSpace=new uint8_t[sizeof(Field_)*y.size()+CACHE_ALIGNMENT];
	size_t spacePtr=(size_t)yTempSpace;
	FieldAXPY<Field_>* yTemp=(FieldAXPY<Field_>*)(spacePtr+CACHE_ALIGNMENT-(spacePtr%CACHE_ALIGNMENT));


#pragma omp parallel
	{
                const Field_& fieldRef=field();
#pragma omp for schedule (static,1024)
                for (size_t i=0;i<y.size();++i) {
                        new ((void*)(yTemp+i)) FieldAXPY<Field_>(fieldRef);
                }

		Index numBlockSizes=rowBlocks_.size();
		for (Index chunkSizeIx=0;chunkSizeIx<numBlockSizes;++chunkSizeIx) {
			const VectorChunks *rowChunks=&(rowBlocks_[chunkSizeIx]);
			Index numChunks=rowChunks->size();
#pragma omp for schedule (static,1)
			for (Index rowChunk=0;rowChunk<numChunks;++rowChunk) {
				const BlockList *blocks=&((*rowChunks)[rowChunk]);
				Index numBlocks=blocks->size();
				for (Index block=0;block<numBlocks;++block) {
                                        const DataBlock *dataBlock=&((*blocks)[block]);
					for (Index k=0;k<dataBlock->elts_.size();++k) {
                                                const Index row=dataBlock->getRow(k);
                                                const Index col=dataBlock->getCol(k);
                                                yTemp[row].mulacc(dataBlock->elts_[k],x[col]);
					}
				}
			}
			//implicit barrier
		}
#pragma omp for schedule (static,1024)
		for (Index i = 0; i < y.size(); ++i) {
			yTemp[i].get(y[i]);
			yTemp[i].~FieldAXPY<Field_>();
		}
	}

	delete[] yTempSpace;
        return y;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector & TriplesBBOMP<Field_>::applyTranspose(OutVector & y, const InVector & x) const
{
	linbox_check( coldim() == y.size() );
	linbox_check( rowdim() == x.size() );

	uint8_t* yTempSpace=new uint8_t[sizeof(Field_)*y.size()+CACHE_ALIGNMENT];
	size_t spacePtr=(size_t)yTempSpace;
	FieldAXPY<Field_>* yTemp=(FieldAXPY<Field_>*)(spacePtr+CACHE_ALIGNMENT-(spacePtr%CACHE_ALIGNMENT));


#pragma omp parallel
	{
                const Field_& fieldRef=field();
#pragma omp for schedule (static,1024)
                for (size_t i=0;i<y.size();++i) {
                        new ((void*)(yTemp+i)) FieldAXPY<Field_>(fieldRef);
                }

		Index numBlockSizes=colBlocks_.size();
		for (Index chunkSizeIx=0;chunkSizeIx<numBlockSizes;++chunkSizeIx) {
			const VectorChunks *colChunks=&(colBlocks_[chunkSizeIx]);
			Index numChunks=colChunks->size();
#pragma omp for schedule (static,1)
			for (Index colChunk=0;colChunk<numChunks;++colChunk) {
				const BlockList *blocks=&((*colChunks)[colChunk]);
				Index numBlocks=blocks->size();
				for (Index block=0;block<numBlocks;++block) {
                                        const DataBlock *dataBlock=&((*blocks)[block]);
					for (Index k=0;k<dataBlock->elts_.size();++k) {
                                                const Index row=dataBlock->getRow(k);
                                                const Index col=dataBlock->getCol(k);
                                                yTemp[col].mulacc(dataBlock->elts_[k],x[row]);
					}
				}
			}
			//implicit barrier
		}
#pragma omp for schedule (static,1024)
		for (Index i = 0; i < y.size(); ++i) {
			yTemp[i].get(y[i]);
			yTemp[i].~FieldAXPY<Field_>();
		}
	}

	delete[] yTempSpace;
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
void TriplesBBOMP<Field_>::splitBlock(RefBlockList& superBlocks,TriplesBlock startBlock)
{
	typedef typename std::vector<Triple>::iterator BlockVecIt;

	std::vector<TriplesBlock> blocks;
	if (startBlock.nnz()<MAX_BLOCK_NNZ) {
		--startBlock.blockEnd_;
		superBlocks.push_back(startBlock);
	} else {
		blocks.push_back(startBlock);
	}

	while(!(blocks.empty())) {
		TriplesBlock block=blocks.back();
		blocks.pop_back();

		Triple blockStart(block.blockStart_);
		Triple blockEnd(block.blockEnd_);
		TriplesCoord blockSize=block.blockSize();
		TriplesCoord quarterSize=blockSize>>2;
		Triple midOne(blockStart.coord_+quarterSize);
		Triple midTwo(midOne.coord_+quarterSize);
		Triple midThree(midTwo.coord_+quarterSize);
		
		BlockVecIt dataBegin=data_.begin()+block.start_;
		BlockVecIt dataEnd=data_.begin()+block.end_;
		Index midOneData,midTwoData,midThreeData;

                Triple midOneMM(midOne), midTwoMM(midTwo), midThreeMM(midThree);
                --midOneMM.coord_; --midTwoMM.coord_; --midThreeMM.coord_;
		midOneData=std::upper_bound(dataBegin,dataEnd,midOneMM,Triple::compareBlockTriples)-data_.begin();
		midTwoData=std::upper_bound(dataBegin,dataEnd,midTwoMM,Triple::compareBlockTriples)-data_.begin();
		midThreeData=std::upper_bound(dataBegin,dataEnd,midThreeMM,Triple::compareBlockTriples)-data_.begin();

		TriplesBlock blockOne(block.start_,midOneData,block.blockStart_,midOne.coord_);
		TriplesBlock blockTwo(midOneData,midTwoData,midOne.coord_,midTwo.coord_);
		TriplesBlock blockThree(midTwoData,midThreeData,midTwo.coord_,midThree.coord_);
		TriplesBlock blockFour(midThreeData,block.end_,midThree.coord_,block.blockEnd_);

		if (blockOne.nnz()<MAX_BLOCK_NNZ) {
			--blockOne.blockEnd_;
			superBlocks.push_back(blockOne);
		} else {
			blocks.push_back(blockOne);
		}
		if (blockTwo.nnz()<MAX_BLOCK_NNZ) {
			--blockTwo.blockEnd_;
			superBlocks.push_back(blockTwo);
		} else {
			blocks.push_back(blockTwo);
		}
		if (blockThree.nnz()<MAX_BLOCK_NNZ) {
			--blockThree.blockEnd_;
			superBlocks.push_back(blockThree);
		} else {
			blocks.push_back(blockThree);
		}
		if (blockFour.nnz()<MAX_BLOCK_NNZ) {
			--blockFour.blockEnd_;
			superBlocks.push_back(blockFour);
		} else {
			blocks.push_back(blockFour);
		}
	}
}

template<class Field_>
void TriplesBBOMP<Field_>::toDataBlock(const RefBlockList& superBlocks,
				       BlockList& dataBlocks)
{
	dataBlocks.clear();
	dataBlocks.resize(superBlocks.size());
	for (size_t i=0;i<superBlocks.size();++i) {
		TriplesBlock block=superBlocks[i];
		int ccf;
		
		Index rowDiff=block.getStartRow()^block.getEndRow();
		Index colDiff=block.getStartCol()^block.getEndCol();
		Index combinedDiff=rowDiff|colDiff;
		if ((combinedDiff>>8)==0) {
			ccf=3;
		} else if ((combinedDiff>>16) == 0) {
			ccf=2;
		} else if ((combinedDiff>>32) == 0) {
			ccf=1;
		} else {
			ccf=0;
		}

		typename DataBlock::TripleListIt startIt=data_.begin()+block.start_;
		typename DataBlock::TripleListIt endIt=data_.begin()+block.end_;
		block.toBlock();
		TriplesCoord blockSize=block.blockSize();
                ++blockSize;
		block.fromBlock();
		dataBlocks[i].init(block.blockStart_,block.blockEnd_,
				   blockSize,startIt,endIt,ccf);
				   
	}
}

template<class Field_>
void TriplesBBOMP<Field_>::finalize()
{
        if ((sortType_ & TRIPLES_SORTED) != 0) {return; }
        for (Index k=0; k<data_.size();++k) {
                data_[k].toBlock();
        }

        std::stable_sort(data_.begin(),data_.end(),Triple::compareBlockTriples);
	std::vector<Triple> tempData;
	tempData.reserve(data_.size());
	if (!(data_.empty())) {
		tempData.push_back(data_[0]);
	}
	for (size_t k=1;k<data_.size();++k) {
		if (!(data_[k].coord_==tempData.back().coord_)) {
			tempData.push_back(data_[k]);
		}
	}
	data_.swap(tempData);

	RefBlockList superBlocks;
        Index rowBound=roundUpIndex(rowdim());
        Index colBound=roundUpIndex(coldim());
        Index rowColBound=(rowBound<colBound)?colBound:rowBound;
        TriplesCoord startBlockIx(0,0),endBlockIx(0,rowColBound);
        coordToBlock(endBlockIx);
        TriplesBlock initialBlock(0,data_.size(),startBlockIx,endBlockIx);

        splitBlock(superBlocks,initialBlock);

        for (Index k=0;k<data_.size();++k) {
                data_[k].fromBlock();
        }

        std::sort(superBlocks.begin(),superBlocks.end(),TriplesBlock::compareBlockSizes);

	for (Index i=0;i<superBlocks.size();++i) {
		superBlocks[i].fromBlock();
	}

	std::vector<DataBlock> dataBlocks;
	toDataBlock(superBlocks,dataBlocks);
	
	rowBlocks_.clear();
        computeVectors(rowBlocks_,dataBlocks,CHUNK_BY_ROW);
	colBlocks_.clear();
        computeVectors(colBlocks_,dataBlocks,CHUNK_BY_COL);

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
	for (Index k = data_.size(); k > 0; --k)
		if (data_[k-1].getRow() == i and data_[k-1].getCol() == j)
			return e = data_[k-1].getElt();
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
