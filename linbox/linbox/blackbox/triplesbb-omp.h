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

 * Written by Rich Seagraves <seagrave@cis.udel.edu>
 * with mods by bds
 */

/** @file blackbox/triplesbb-omp.h
 * @ingroup blackbox
 * @brief NO DOC
 */

#ifndef __LINBOX_triplesbb_omp_H
#define __LINBOX_triplesbb_omp_H

#include <limits.h>
#include <algorithm>
#include <iostream>
#include <map>
using std::istream;
using std::ostream;
using std::max;
#include "linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/field/hom.h"

#include <vector>

#if ULONG_MAX > 4294967295UL
#define USING_64_BIT_SIZE_T 1
#endif

namespace LinBox
{

/** \brief wrapper for NAG Sparse Matrix format.
 *
 \ingroup blackbox
 * Sparse matrix representation which stores nonzero entries by i,j,value triples.
 */
template<class Field_>
class TriplesBBOMP : public BlackboxInterface {

	public:
        typedef Field_ Field;
	typedef typename MatrixDomain<Field>::Matrix Matrix;
	typedef typename Field::Element Element;
	typedef TriplesBBOMP<Field> Self_t;
	typedef size_t Index; // would prefer a signed type
	typedef long long Llong;

	// Default constructor.
	TriplesBBOMP();

	TriplesBBOMP(const TriplesBBOMP & B);

	TriplesBBOMP & operator=(const TriplesBBOMP & B);

	TriplesBBOMP(const Field& F, istream& in);

	istream& read(istream& in);

	ostream& write(ostream& out);

	~TriplesBBOMP();

	TriplesBBOMP(const Field& F, Index r = 0, Index c = 0);

        void finalize();

	// (re)shape the matrix.  Any prior entries are abandoned.
	TriplesBBOMP& shape(const Field& F, Index r = 0, Index c = 0);

	// need cstor from matrix stream, read, write

	// Element e is added in the i,j position.
	void setEntry(Index i, Index j, const Element & e);

	// Element e is set to the i,j entry.
	Element& getEntry(Element& e, Index i, Index j) const;

	/** y <- A x.
	 *
	 *  Performance will be better if A is in rowMajor or colMajor order.
	 *
	 *  If this were to be used extensively for sparse black box ops,
	 *  optimizations would be desirable.
	 */
	template<class OutVector, class InVector>
	OutVector & apply(OutVector &, const InVector &) const;

	template<class OutVector, class InVector>
	OutVector & seqApply(OutVector &, const InVector &) const;

	/// Y <- AX, requires conformal shapes.
	//Matrix & apply_right(Matrix &Y, const Matrix &X) const;

	/** y <- A^T x.
	 *
	 *  Performance will be better if A is in rowMajor or colMajor order.
	 *
	 *  If this were to be used extensively for sparse black box ops,
	 *  optimizations would be desirable.
	 */
	template<class OutVector, class InVector>
	OutVector & applyTranspose(OutVector &, const InVector &) const;

	size_t rowdim() const;

	size_t coldim() const;

	const Field & field() const;

	/* Returns number of non-zero entries */
	size_t size() const;

	template<typename Tp1_>
	struct rebind {
		typedef TriplesBBOMP<Tp1_> other;
		void operator() (other & Ap, const Self_t& A, const Tp1_& F)
		{
			Hom <typename Self_t::Field, Tp1_> hom( A.field(), F);

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

protected:

	const static Index BLOCK_SIZE_BITS=10;
	const static Index BLOCK_MASK=~((1<<BLOCK_SIZE_BITS)-1);
	const static Index MAX_BLOCK_NNZ=4096;
	const static Index EXPECTED_THREADS=64;
	void sortBlock();

	void mergeBlocks(Index maxSize,std::vector<int>& mergeCounts, std::vector<bool>& mergeFlags);

	MatrixDomain<Field> MD_;

	union Coord {
		Index rowCol[2];
#ifdef USING_64_BIT_SIZE_T
		Llong blockIxArr[2]; // little endian
#endif
		Llong blockIx;
	};

	static void toBlock(Coord& coord);

	static void fromBlock(Coord& coord);

        struct Triple {
                Triple(Index& r, Index& c, const Element& e){init(r, c, e);}
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
		void toBlock() {TriplesBBOMP<Field_>::toBlock(coord_);}
		void fromBlock() {TriplesBBOMP<Field_>::fromBlock(coord_);}
                inline Llong& blockIx() {return coord_.blockIx;}
		static bool compareBlockTriples(const Triple& lhs,
						const Triple& rhs) {
#ifdef USING_64_BIT_SIZE_T
			return (lhs.coord_.blockIxArr[1]<rhs.coord_.blockIxArr[1])||
				((!(lhs.coord_.blockIxArr[1]>rhs.coord_.blockIxArr[1])) &&
				 (lhs.coord_.blockIxArr[0]<rhs.coord_.blockIxArr[0]));
#else
			return lhs.coord_.blockIx<rhs.coord_.blockIx;
#endif
		}
	protected:
                union Coord coord_;
                Element elt_;
        };

	std::vector<Triple> data_;

	Index rows_, cols_;

#define TRIPLES_UNSORTED 1
#define TRIPLES_SORTED 2

	//Either TRIPLES_SORTED or TRIPLES_UNSORTED, maybe extend this later
	int sortType_;

	struct Block {
		Index start_,end_;
		Coord blockStart_, blockEnd_;
		Index blockSize_;

		Block(Index start,Index end,Index blockIx) :
			start_(start),end_(end),blockSize_(1) {
			blockStart_.blockIx=blockIx;
			blockEnd_.blockIx=blockIx;
		}

                Index nnz() const {return end_-start_;}
                Index blockSize() const {return blockSize_;}
		void expand(Index newDataEnd, Index newBlockEnd) {
			end_=newDataEnd;
			blockEnd_.blockIx=newBlockEnd;
			blockSize_=1+blockEnd_.blockIx-blockStart_.blockIx;
		}
		void toBlock() {
			TriplesBBOMP<Field_>::toBlock(blockStart_);
			TriplesBBOMP<Field_>::toBlock(blockEnd_);
		}
		void fromBlock() {
			TriplesBBOMP<Field_>::fromBlock(blockStart_);
			TriplesBBOMP<Field_>::fromBlock(blockEnd_);
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
	};

	void computeRowBlocks();

        typedef std::vector<std::vector<Block> > BlockList;

        std::vector<Block> superBlocks_;

        BlockList rowBlocks_;

        BlockList colBlocks_;
  }; // TriplesBBOMP

} // namespace LinBox

#include "linbox/blackbox/triplesbb-omp.inl"

#endif // __LINBOX_triplesbb_omp_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
