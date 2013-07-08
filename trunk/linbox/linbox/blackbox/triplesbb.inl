/* linbox/blackbox/triplesbb.inl
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

#ifndef __LINBOX_triplesbb_INL
#define __LINBOX_triplesbb_INL

#include <algorithm>
#include <iostream>
using std::istream;
using std::ostream;
#include "linbox-config.h"
#include "linbox/util/debug.h"
#include "linbox/util/field-axpy.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/field/hom.h"
#include "linbox/util/matrix-stream.h"

#include <omp.h>

#include <vector>

namespace LinBox
{

template<class MatDom> TriplesBB<MatDom>::
TriplesBB() {}

template<class MatDom> TriplesBB<MatDom>::
~TriplesBB() {}

template<class MatDom> TriplesBB<MatDom>::
TriplesBB(const Field& F, istream& in) : field_(&F), MD_(F) { read(in); }

template<class MatDom> istream& TriplesBB<MatDom>::
read(istream& in){
	Index r, c;
	typename Field::Element v; field().init(v);
	MatrixStream<Field> ms(field(), in);
	ms.getDimensions(r, c);
	init(field(), r, c);
	while (ms.nextTriple(r, c, v)) setEntry(r, c, v);
	return in;
}

template<class MatDom> ostream& TriplesBB<MatDom>::
write(ostream& out){
	Index r, c;
	out << "%%MatrixMarket matrix coordinate integer general" << std::endl;
	out << "% written from a LinBox TriplesBB" << std::endl;
	out << rowdim() <<" " << coldim() << " " << size() << std::endl;
	for (Index k = 0; k < size(); ++k) {
		Triple t = data_[k];
		field().write(out << t.row+1 << " " << t.col+1 << " ", t.elt) << std::endl;
	}
	return out;
}

template<class MatDom> TriplesBB<MatDom>& TriplesBB<MatDom>::
init(const Field& F, Index r, Index c)
{ field_ = &F; MD_.init(F), data_.clear(); rows_ = r; cols_ = c; rowMajor_ = 0; return *this; }

template<class MatDom> TriplesBB<MatDom>::
TriplesBB(const Field& F, Index r, Index c)
: field_(&F), data_(), rows_(r), cols_(c), rowMajor_(0) {}

template<class MatDom> TriplesBB<MatDom>::
TriplesBB(const TriplesBB<MatDom> & B)
: field_ ( B.field_ ), MD_(B.field()), data_ ( B.data_ ), rows_ ( B.rows_ ), cols_ ( B.cols_ ),
   rowMajor_ ( B.rowMajor_ )
{}

template<class MatDom> TriplesBB<MatDom>& TriplesBB<MatDom>::
operator=(const TriplesBB<MatDom> & rhs)
{	if (rhs == this) return ;
	field_ = rhs.field_;
	MD_.init(rhs.field_);
	data_ = rhs.data_;
	rows_ = rhs.rows_;
	cols_ = rhs.cols_;
	rowMajor_ = rhs.rowMajor_;
	return *this;
}

template<class MatDom> typename TriplesBB<MatDom>::Matrix& TriplesBB<MatDom>::
apply_left // Y = AX
	( typename TriplesBB<MatDom>::Matrix &Y,
	  const typename TriplesBB<MatDom>::Matrix &X
	) const 
{	Y.zero();
	typename MatrixDomain::Submatrix Yc, Xc;// col submatrices
	for (Index k = 0; k < data_.size(); ++k) {
		Triple t = data_[k];
		Yc.submatrix(Y,t.row,0,1,Y.coldim());
		Xc.submatrix(X,t.col,0,1,X.coldim());
		MD_.axpyin(Yc, t.elt, Xc);
	}
	return Y;
}

template<class MatDom> typename TriplesBB<MatDom>::Matrix& TriplesBB<MatDom>::
apply_right // Y = XA
	( typename TriplesBB<MatDom>::Matrix &Y,
	  const typename TriplesBB<MatDom>::Matrix &X
	) const 
{	Y.zero();
	typename MatrixDomain::Submatrix Yr, Xr; // row submatrices
	for (Index k = 0; k < data_.size(); ++k) {
		Triple t = data_[k];
		Yr.submatrix(Y,0,t.col,Y.rowdim(),1);
		Xr.submatrix(X,0,t.row,X.rowdim(),1);
		MD_.axpyin(Yr, t.elt, Xr);
	}
	return Y;
}

template<class MatDom> 
template<class OutVector, class InVector>
OutVector & TriplesBB<MatDom>::apply(OutVector & y, const InVector & x) const
{
	linbox_check( rowdim() == y.size() );
	linbox_check( coldim() == x.size() );
        double start = omp_get_wtime();
	for (Index i = 0; i < y.size(); ++i) field().init(y[i], field().zero);
	for (Index k = 0; k < data_.size(); ++k) {
		Triple t = data_[k];
		field().axpyin(y[t.row], t.elt, x[t.col]);
	}
	return y;
}

template<class MatDom> 
template<class OutVector, class InVector>
OutVector & TriplesBB<MatDom>::applyTranspose(OutVector & y, const InVector & x) const
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

template<class MatDom> size_t TriplesBB<MatDom>::
rowdim() const { return rows_; }

template<class MatDom> size_t TriplesBB<MatDom>::
coldim() const { return cols_; }

template<class MatDom> const typename MatDom::Field& TriplesBB<MatDom>::
field() const { return *field_; }

template<class MatDom> const MatDom& TriplesBB<MatDom>::
domain() const { return MD_; }

template<class MatDom> size_t TriplesBB<MatDom>::
size() const { return data_.size(); }

template<class MatDom> void TriplesBB<MatDom>::
setEntry(Index i, Index j, const typename Field::Element & e)
{
	rowMajor_ = false;
	data_.push_back(Triple(i, j, e));
}

template<class MatDom> typename MatDom::Field::Element& TriplesBB<MatDom>::
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
