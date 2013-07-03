/* linbox/blackbox/triplesbb.h
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

/** @file blackbox/triplesbb.h
 * @ingroup blackbox
 * @brief NO DOC
 */

#ifndef __LINBOX_triplesbb_H
#define __LINBOX_triplesbb_H

#include <algorithm>
#include <iostream>
using std::istream;
using std::ostream;
using std::max;
#include "linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/field/hom.h"

#include <vector>

namespace LinBox
{

	/** \brief wrapper for NAG Sparse Matrix format.
	 *
	 \ingroup blackbox
	 * Sparse matrix representation which stores nonzero entries by i,j,value triples.
	 */
  template<class MatDom>
  class TriplesBB : public BlackboxInterface {

	public:
	typedef MatDom MatrixDomain;
	typedef typename MatrixDomain::Field Field;
	typedef typename MatrixDomain::Submatrix Matrix;
	typedef typename Field::Element Element;
	typedef TriplesBB<MatrixDomain> Self_t;
	typedef size_t Index; // would prefer a signed type

	// Default constructor.
	TriplesBB();

	TriplesBB(const TriplesBB & B);

	TriplesBB & operator=(const TriplesBB & B);

	TriplesBB(const Field& F, istream& in);

	istream& read(istream& in);

	ostream& write(ostream& out);

	~TriplesBB();

	TriplesBB(const Field& F, Index r = 0, Index c = 0);

	// (re)initialize the matrix.  Any prior entries are abandoned.
	TriplesBB& init(const Field& F, Index r = 0, Index c = 0);

	// need cstor from matrix stream, read, write

	// Element e is added in the i,j position.
	void setEntry(Index i, Index j, const Element & e);

	// Element e is set to the i,j entry.
	Element& getEntry(Element& e, Index i, Index j) const;

	/// Y <- AX, requires conformal shapes.
	Matrix & apply_right(Matrix &Y, const Matrix &X) const;

	/// Y <- XA, requires conformal shapes.
	Matrix & apply_left(Matrix &Y, const Matrix &X) const;

	/** y <- A x.
	 *
	 *  Performance will be better if A is in rowMajor or colMajor order.
	 *
	 *  If this were to be used extensively for sparse black box ops,
	 *  optimizations would be desirable.
	 */
	template<class OutVector, class InVector>
	OutVector & apply(OutVector &y, const InVector &x) const;

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

	const MatrixDomain & domain() const;

	/* Returns number of non-zero entries */
	size_t size() const;

	template<typename Tp1_>
	struct rebind {
		typedef TriplesBB<Tp1_> other;
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
	const Field *field_; // The field of the entries
	MatrixDomain MD_; // The matrix domain for apply_left and apply_right

	struct Triple { Index row; Index col; Element elt;
		Triple(Index& r, Index& c, const Element& e)
		{ init(r, c, e); }
		void init(Index& r, Index& c, const Element& e)
		{ row = r; col = c; elt = e; }
	};

	// the data
	std::vector<Triple> data_;

	/// The number of rows, columns
	Index rows_, cols_;

	int rowMajor_; // 1 = rowMajor, -1 = colMajor, 0 = unsorted.

  }; // TriplesBB

} // namespace LinBox

#include "linbox/blackbox/triplesbb.inl"

#endif // __LINBOX_triplesbb_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
