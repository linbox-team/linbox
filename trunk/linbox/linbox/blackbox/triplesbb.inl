/* linbox/blackbox/triplesbb.inl
 * Copyright (c) LinBox
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
 * with mods by bds and Alex Stachnik
 */

/** @file blackbox/triplesbb.h
 * @ingroup blackbox
 * @brief A COO (vector of i,j,v triples) sparse matrix rep.
 */

#ifndef __LINBOX_triplesbb_INL
#define __LINBOX_triplesbb_INL

#include <algorithm>
#include <iostream>
//#include <vector>

#include "linbox/util/matrix-stream.h"

namespace LinBox
{

// template<class Field_> TriplesBB<Field_>::
// TriplesBB() : MD_(), data_(), rows_(0), cols_(0), sort_(unsorted) {}

template<class Field_> TriplesBB<Field_>::
~TriplesBB() {}

template<class Field_> TriplesBB<Field_>::
TriplesBB(const Field& F, std::istream& in)
: MD_(F), data_(), rows_(0), cols_(0), sort_(unsorted)
{ read(in); }

template<class Field_> std::istream& TriplesBB<Field_>::
read(std::istream& in){
	Index r, c;
	typename Field::Element v; field().init(v);
	MatrixStream<Field> ms(field(), in);
	ms.getDimensions(r, c);
	init(field(), r, c);
	while (ms.nextTriple(r, c, v)) setEntry(r, c, v);
	return in;
}

template<class Field_> std::ostream& TriplesBB<Field_>::
write(std::ostream& out){
	out << "%%MatrixMarket matrix coordinate integer general" << std::endl;
	out << "% written from a LinBox TriplesBB" << std::endl;
	out << rowdim() <<" " << coldim() << " " << size() << std::endl;
	for (Index k = 0; k < size(); ++k) {
		Triple t = data_[k];
		field().write(out << t.row+1 << " " << t.col+1 << " ", t.elt) << std::endl;
	}
	return out;
}

template<class Field_> TriplesBB<Field_>& TriplesBB<Field_>::
init(const Field& F, Index r, Index c)
{ MD_.init(F), data_.clear(); rows_ = r; cols_ = c; sort_ = unsorted; return *this; }

template<class Field_> TriplesBB<Field_>::
TriplesBB(const Field& F, Index r, Index c)
: MD_(F), data_(), rows_(r), cols_(c), sort_(unsorted) {}

template<class Field_> TriplesBB<Field_>::
TriplesBB(const TriplesBB<Field_> & B)
: MD_(B.field()), data_ ( B.data_ ), rows_ ( B.rows_ ), cols_ ( B.cols_ ),
   sort_ ( B.sort_ )
{}

template<class Field_> TriplesBB<Field_>& TriplesBB<Field_>::
operator=(const TriplesBB<Field_> & rhs)
{	if (rhs == this) return ;
	MD_.init(rhs.field_);
	data_ = rhs.data_;
	rows_ = rhs.rows_;
	cols_ = rhs.cols_;
	sort_ = rhs.sort_;
	return *this;
}

template<class Field_>
template<class Mat1, class Mat2> Mat1& TriplesBB<Field_>::
// I do not like this need to templatize -bds
//typename TriplesBB<Field_>::Matrix& TriplesBB<Field_>::
applyLeft // Y = AX
	( /*typename TriplesBB<Field_>::Matrix*/Mat1 &Y,
	  const /*typename TriplesBB<Field_>::Matrix*/Mat2 &X
	) const
{	Y.zero();
	Matrix Yc, Xc;// row submatrices
	for (Index k = 0; k < data_.size(); ++k) {
		Triple t = data_[k];
		Yc.submatrix(Y,t.row,0,1,Y.coldim());
		Xc.submatrix(X,t.col,0,1,X.coldim());
		MD_.saxpyin(Yc, t.elt, Xc);
	}
	return Y;
}

//template<class Field_> typename TriplesBB<Field_>::Matrix& TriplesBB<Field_>::
template<class Field_>
template<class Mat1, class Mat2> Mat1& TriplesBB<Field_>::
applyRight // Y = XA
	( /*TriplesBB<Field_>::Matrix*/Mat1 &Y,
	  const /*typename TriplesBB<Field_>::Matrix*/Mat2 &X
	) const
{	Y.zero();
	Matrix Yr, Xr; // row submatrices
	for (Index k = 0; k < data_.size(); ++k) {
		Triple t = data_[k];
		Yr.submatrix(Y,0,t.col,Y.rowdim(),1);
		Xr.submatrix(X,0,t.row,X.rowdim(),1);
		MD_.saxpyin(Yr, t.elt, Xr);
	}
	return Y;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector & TriplesBB<Field_>::apply(OutVector & y, const InVector & x) const
{
	linbox_check( rowdim() == y.size() );
	linbox_check( coldim() == x.size() );
	for (Index i = 0; i < y.size(); ++i) field().init(y[i], field().zero);
	for (Index k = 0; k < data_.size(); ++k) {
		Triple t = data_[k];
		field().axpyin(y[t.row], t.elt, x[t.col]);
	}
	return y;
}

template<class Field_>
template<class OutVector, class InVector>
OutVector & TriplesBB<Field_>::applyTranspose(OutVector & y, const InVector & x) const
{
	linbox_check( coldim() == y.size() );
	linbox_check( rowdim() == x.size() );
	for (Index i = 0; i < y.size(); ++i) field().init(y[i], field().zero);
	for (Index k = 0; k < data_.size(); ++k) {
		const Triple& t = data_[k];
		field().axpyin(y[t.col], t.elt, x[t.row]);
	}
	return y;
}

template<class Field_> size_t TriplesBB<Field_>::
rowdim() const { return rows_; }

template<class Field_> size_t TriplesBB<Field_>::
coldim() const { return cols_; }

template<class Field_> const Field_& TriplesBB<Field_>::
field() const { return MD_.field(); }

template<class Field_> size_t TriplesBB<Field_>::
size() const { return data_.size(); }

template<class Field_> void TriplesBB<Field_>::
setEntry(Index i, Index j, const typename Field::Element & e)
{
	sort_ = unsorted;
	data_.push_back(Triple(i, j, e));
}

template<class Field_> void TriplesBB<Field_>::
finalize(sortPolicy s) { /* sort according to policy */ sort_ = s; }

template<class Field_> typename Field_::Element& TriplesBB<Field_>::
getEntry(typename Field_::Element& e, Index i, Index j) const
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
