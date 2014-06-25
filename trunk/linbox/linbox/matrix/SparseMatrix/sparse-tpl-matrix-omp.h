/* linbox/blackbox/triplesbb-omp.h
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

 * Written by Alex Stachnik <stachnik@udel.edu>
 */

/*! @file matrix/SparseMatrix/sparse-tpl-matrix-omp.h
 * @ingroup sparsematrix
 * @ingroup omp
 * @brief NO DOC
 */

#ifndef __LINBOX_matrix_sparsematrix_sparse_tpl_matrix_omp_H
#define __LINBOX_matrix_sparsematrix_sparse_tpl_matrix_omp_H


#include <algorithm>
#include <iostream>
#include <map>
#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/field/hom.h"
#include "linbox/matrix/SparseMatrix/triples-coord.h"

#include <vector>


namespace LinBox
{

union TriplesSmallCoord {
	uint8_t c[8];
	uint16_t s[4];
	uint32_t i[2];
	uint64_t l;
};

template <class Element>
struct TriplesBBTriple {
	TriplesBBTriple(TriplesCoord c):coord_(c) {}
	TriplesBBTriple(Index r, Index c, const Element& e){init(r, c, e);}
	void init(Index r, Index c, const Element& e)
	{
		row()=r;
		col()=c;
		elt_ = e;
	}
	inline Index& row() {return coord_.rowCol[0];}
	inline Index& col() {return coord_.rowCol[1];}
	inline Element& elt() {return elt_;}
	inline const Index& getRow() const {return coord_.rowCol[0];}
	inline const Index& getCol() const {return coord_.rowCol[1];}
	inline const Element& getElt() const {return elt_;}
	void toBlock() {coordToBlock(coord_);}
	void fromBlock() {coordFromBlock(coord_);}
	static bool compareBlockTriples(const TriplesBBTriple<Element>& lhs,
					const TriplesBBTriple<Element>& rhs) {
		return lhs.coord_<rhs.coord_;
	}
	TriplesCoord coord_;
	Element elt_;
};

struct TriplesBlock {
	Index start_,end_;
	TriplesCoord blockStart_, blockEnd_;

	TriplesBlock(Index start,Index end,TriplesCoord blockStart,TriplesCoord blockEnd) :
		start_(start),end_(end),
		blockStart_(blockStart),blockEnd_(blockEnd) {
	}
	Index nnz() const {return end_-start_;}
	TriplesCoord blockSize() const {return blockEnd_-blockStart_;}
	void toBlock() {
		coordToBlock(blockStart_);
		coordToBlock(blockEnd_);
	}
	void fromBlock() {
		coordFromBlock(blockStart_);
		coordFromBlock(blockEnd_);

	}
	inline const Index& getStartRow() const {
		return blockStart_.rowCol[0];
	}
	inline const Index& getStartCol() const {
		return blockStart_.rowCol[1];
	}
	inline const Index& getEndRow() const {
		return blockEnd_.rowCol[0];
	}
	inline const Index& getEndCol() const {
		return blockEnd_.rowCol[1];
	}
	static bool compareBlockSizes(const TriplesBlock& lhs,
				      const TriplesBlock& rhs) {
		return lhs.blockSize()<rhs.blockSize();
	}
};


template <class Element>
struct TriplesDataBlock {
	TriplesCoord blockStart_, blockEnd_, blockSize_;
	std::vector<TriplesSmallCoord> rowIxs_,colIxs_;
	std::vector<Element> elts_;
	int ccf_;

	typedef TriplesBBTriple<Element> Triple;
	typedef typename std::vector<Triple> TripleList;
	typedef typename TripleList::const_iterator TripleListIt;

        inline Index getRow(int ix) const {
                Index retVal=getStartRow();
                const TriplesSmallCoord *rowIxPtr=&(rowIxs_[ix>>ccf_]);
                switch (ccf_) {
                case 0:
                        return rowIxPtr->l;
                case 1:
                        return (retVal|(rowIxPtr->i[ix&1]));
                case 2:
                        return (retVal|(rowIxPtr->s[ix&3]));
                case 3:
                default:
                        return (retVal|(rowIxPtr->c[ix&7]));
                }
        }
        inline Index getCol(int ix) const {
                Index retVal=getStartCol();
                const TriplesSmallCoord *colIxPtr=&(colIxs_[ix>>ccf_]);
                switch (ccf_) {
                case 0:
                        return colIxPtr->l;
                case 1:
                        return (retVal|(colIxPtr->i[ix&1]));
                case 2:
                        return (retVal|(colIxPtr->s[ix&3]));
                case 3:
                default:
                        return (retVal|(colIxPtr->c[ix&7]));
                }
        }

	void init(TriplesCoord blockStart,
		  TriplesCoord blockEnd,
		  TriplesCoord myBlockSize, // blockSize shadows blockSize()
		  TripleListIt triplesStart,
		  TripleListIt triplesEnd,
		  int ccf) {
		blockStart_=blockStart;
		blockEnd_=blockEnd;
		blockSize_=myBlockSize;
		ccf_=ccf;
                int numElts=(int)(triplesEnd-triplesStart);
		rowIxs_.resize(1+(numElts>>ccf));
		colIxs_.resize(1+(numElts>>ccf));
		elts_.reserve(numElts);

		int coordOffset=0, coordIx=0;
		while (triplesStart != triplesEnd) {
			elts_.push_back(triplesStart->getElt());

			switch (ccf) {
			case 0:
				rowIxs_[coordIx].l=(uint64_t)(triplesStart->getRow());
				colIxs_[coordIx].l=(uint64_t)(triplesStart->getCol());
				break;
			case 1:
				rowIxs_[coordIx].i[coordOffset]=
					(uint32_t)(triplesStart->getRow());
				colIxs_[coordIx].i[coordOffset]=
					(uint32_t)(triplesStart->getCol());
				break;
			case 2:
				rowIxs_[coordIx].s[coordOffset]=
					(uint16_t)(triplesStart->getRow());
				colIxs_[coordIx].s[coordOffset]=
					(uint16_t)(triplesStart->getCol());
				break;
			case 3:
				rowIxs_[coordIx].c[coordOffset]=
					(uint8_t)(triplesStart->getRow());
				colIxs_[coordIx].c[coordOffset]=
					(uint8_t)(triplesStart->getCol());
				break;
			}

			++coordOffset;
			if (coordOffset==(1<<ccf)) {
				coordOffset=0;
				++coordIx;
			}
			++triplesStart;
		}
	}

	inline TriplesCoord blockSize() const {
		return blockSize_;
	}
	inline Index getStartRow() const {
		return blockStart_.rowCol[0];
	}
	inline Index getStartCol() const {
		return blockStart_.rowCol[1];
	}
	inline Index getEndRow() const {
		return blockEnd_.rowCol[0];
	}
	inline Index getEndCol() const {
		return blockEnd_.rowCol[1];
	}
};


/** \brief wrapper for NAG Sparse Matrix format.
 *
 \ingroup blackbox
 * Sparse matrix representation which stores nonzero entries by i,j,value triples.
 */
template<class Field_>
class SparseMatrix<Field_, SparseMatrixFormat::TPL_omp> : public BlackboxInterface {

	public:
        typedef Field_ Field;
	typedef typename MatrixDomain<Field>::Matrix Matrix;
	typedef typename Field::Element Element;
	typedef SparseMatrix<Field, SparseMatrixFormat::TPL_omp> Self_t;

	// Default constructor.
	SparseMatrix();

	SparseMatrix(const SparseMatrix & B);

	SparseMatrix & operator=(const SparseMatrix & B);

	SparseMatrix(const Field& F, std::istream& in);

	std::istream& read(std::istream& in);

	std::ostream& write(std::ostream& out);

	~SparseMatrix();

	SparseMatrix(const Field& F, Index r = 0, Index c = 0);

        void finalize();

	// (re)shape the matrix.  Any prior entries are abandoned.
	SparseMatrix& shape(const Field& F, Index r = 0, Index c = 0);

	// need cstor from matrix stream, read, write

	// Element e is added in the i,j position.
	void setEntry(Index i, Index j, const Element & e);

	// Element e is set to the i,j entry.
	Element& getEntry(Element& e, Index i, Index j) const;

	/** y <- A x.
	 */
	template<class OutVector, class InVector>
	OutVector & apply(OutVector &, const InVector &) const;

        // Mul with this on left: Y <- AX. Requires conformal shapes.
	// Requires Y != X
	template<class Mat1, class Mat2>
	Mat1 & applyLeft(Mat1 &Y, const Mat2 &X) const;

        // Mul with this on right: Y <- XA. Requires conformal shapes.
	// Requires Y != X
	template<class Mat1, class Mat2>
	Mat1 & applyRight(Mat1 &Y, const Mat2 &X) const;

	/** y <- A^T x.
	 */
	template<class OutVector, class InVector>
	OutVector & applyTranspose(OutVector &, const InVector &) const;

	Index rowdim() const;

	Index coldim() const;

	const Field & field() const;

	/* Returns number of non-zero entries */
	size_t size() const;

	template<typename Tp1_>
	struct rebind {
		typedef SparseMatrix<Tp1_> other;
		// BB : is this working at all ?
		void operator() (other & Ap, const Self_t& A)
		{
			Hom <typename Self_t::Field, Tp1_> hom( A.field(), Ap.field());

			typedef typename Tp1_::Element otherElt;
			typedef typename std::vector<otherElt> othervec;
			typedef typename std::vector<Element> selfvec;
			typedef typename othervec::iterator otheriter;
			typedef typename selfvec::const_iterator selfiter;
			otheriter vp_p; selfiter v_p;

			Ap.data_.resize(A.data.size());
			for (v_p = A.data_.begin(), vp_p = Ap.data_.begin();
			     v_p != A.data.end(); ++ v_p, ++ vp_p)
				hom.image (vp_p->elt, v_p->elt);
		}
	};

        //For debugging:
        typename std::vector<std::vector<std::vector<TriplesDataBlock<Element> > > >& getRowBlocks() {
                return rowBlocks_;
        }

protected:

	const static Index MAX_BLOCK_NNZ=1024;
        const static size_t CACHE_ALIGNMENT=64;

        inline static int roundUpIndex(Index val) {
                int exponent=0;
                --val;
                while (val > 0) {
                        ++exponent;
                        val=val>>1;
                }
                return 1<<exponent;
        }

	typedef TriplesBBTriple<Element> Triple;
	typedef TriplesDataBlock<Element> DataBlock;

        typedef std::vector<TriplesBlock> RefBlockList;

        typedef std::vector<DataBlock> BlockList;
        typedef typename BlockList::iterator BlockListIt;
        typedef std::vector<BlockList> VectorChunks;
        typedef typename VectorChunks::iterator VectorChunkIt;
        typedef std::vector<VectorChunks> SizedChunks;

	typedef std::map<Index,Index> IntervalSet;
	typedef typename IntervalSet::iterator IntervalIterator;
	typedef std::pair<Index,Index> Interval;

	void splitBlock(RefBlockList &superBlocks, TriplesBlock block);

	void toDataBlock(const RefBlockList& superBlocks,
			 BlockList& dataBlocks);

	void computeVectors(SizedChunks &sizedChunks,
			    BlockList &superBlocks,
			    const int rowOrCol);

        void combineIntervals(BlockListIt startIt,
                              BlockListIt endIt,
                              IntervalSet& intervals,
                              VectorChunks& chunksOut,
                              const int rowOrCol);

        void nonOverlappingIntervals(BlockListIt startIt,
                                     BlockListIt endIt,
                                     IntervalSet& intervals,
                                     const int rowOrCol);

	MatrixDomain<Field> MD_;

	std::vector<Triple> data_;

	Index rows_, cols_;

#define TRIPLES_UNSORTED 1
#define TRIPLES_SORTED 2

	//Either TRIPLES_SORTED or TRIPLES_UNSORTED
	int sortType_;

#define CHUNK_BY_ROW 1
#define CHUNK_BY_COL 2

        SizedChunks rowBlocks_;

        SizedChunks colBlocks_;
  }; // SparseMatrix

} // namespace LinBox

#include "sparse-tpl-matrix-omp.inl"

#endif // __LINBOX_matrix_sparsematrix_sparse_tpl_matrix_omp_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
