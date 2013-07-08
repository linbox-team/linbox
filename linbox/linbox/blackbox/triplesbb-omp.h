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

// TODO: Grab this from configure, or implement it in a smarter way
#define USING_64_BIT_SIZE_T 1

namespace LinBox
{

/** \brief wrapper for NAG Sparse Matrix format.
 *
 \ingroup blackbox
 * Sparse matrix representation which stores nonzero entries by i,j,value triples.
 */
template<class MatDom>
class TriplesBBOMP : public BlackboxInterface {

	public:
        typedef MatDom MatrixDomain;
	typedef typename MatrixDomain::Field Field;
	typedef typename MatrixDomain::Submatrix Matrix;
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
	Matrix & apply_right(Matrix &Y, const Matrix &X) const;

	/** y <- A^T x.
	 *
	 *  Performance will be better if A is in rowMajor or colMajor order.
	 *
	 *  If this were to be used extensively for sparse black box ops,
	 *  optimizations would be desirable.
	 */
	template<class OutVector, class InVector>
	OutVector & applyTranspose(OutVector &, const InVector &) const;

	void sortBlock();

	void sortRow();

	void sortColumn();

	size_t rowdim() const;

	size_t coldim() const;

	const Field & field() const;

        const MatrixDomain & domain() const;

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
	const Field *field_; // The field used by this class
        MatrixDomain MD_;

        union Coord {
                Index rowCol[2];
#ifdef USING_64_BIT_SIZE_T
                Llong blockIxArr[2]; // little endian
#endif
                Llong blockIx;
        };

        struct Triple {
                union Coord coord;
                Element elt;
                Triple(Index& r, Index& c, const Element& e)
                { init(r, c, e); }
                void init(Index r, Index c, const Element& e)
                {
                        row()=r;
                        col()=c;
                        elt = e;
                }
                inline Index& row() {return coord.rowCol[0];}
                inline Index& col() {return coord.rowCol[1];}
                inline const Index& getRow() const {return coord.rowCol[0];}
                inline const Index& getCol() const {return coord.rowCol[1];}
                inline Llong& blockIx() {return coord.blockIx;}
                void toBlock()
                {
                        Llong temp,final;
                        Index localRow=row(),localCol=col();
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

                void fromBlock()
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
                        row()=localRow;
                        col()=localCol;
                }
                static bool compareBlockTriples(const Triple& left, const Triple& right) {
#ifdef USING_64_BIT_SIZE_T
                        return (left.coord.blockIxArr[1]<right.coord.blockIxArr[1])||
                                ((!(left.coord.blockIxArr[1]>right.coord.blockIxArr[1])) &&
                                 (left.coord.blockIxArr[0]<right.coord.blockIxArr[0]));
#else
                        return left.coord.blockIx<right.coord.blockIx;
#endif
                }
        };

	// the data
	std::vector<Triple> data_;

	/// The number of rows, columns
	Index rows_, cols_;

        #define TRIPLES_UNSORTED 1
        #define TRIPLES_ROW_MAJOR_SORT 2
        #define TRIPLES_COL_MAJOR_SORT 4
        #define TRIPLES_BLOCK_SORT 8

	int sortType_;

	struct Block {
		Index start_,end_;
		Block(Index start,Index end) : start_(start),end_(end) {}
	};

        typedef std::vector<std::vector<Block> > BlockList;

        BlockList colBlocks_;

        BlockList rowBlocks_;

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
