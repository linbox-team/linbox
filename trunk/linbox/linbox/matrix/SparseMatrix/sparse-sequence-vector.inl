/* linbox/matrix/sparse.inl
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               1999-2001 William J Turner,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Based on sparse-base.h by William J Turner <wjturner@math.ncsu.edu>
 *
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/sparse-base.inl to matrix/sparse.inl
 * ------------------------------------
 * 2002-11-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 *   - Renamed ColOfRowsIterator to RowIterator
 *   - Named template argument _Row rather than Row; add a typedef to Row
 * ------------------------------------
 *
 *
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
 *.
 */

#ifndef __LINBOX_matrix_sparse_sequence_INL
#define __LINBOX_matrix_sparse_sequence_INL

namespace LinBox { namespace Protected {

	template <class Field, class Row>
	SparseMatrixGeneric<Field,Row,VectorCategories::SparseSequenceVectorTag> ::SparseMatrixGeneric( MatrixStream<Field>& ms ) :
		_field(ms.field())
		,_MD(ms.field()),_AT(*this)
		,_matA(0), _m(0), _n(0)

	{
		Element val;
		size_t i, j;
		while( ms.nextTriple(i,j,val) ) {
			if( i >= _m ) {
				_m = i + 1;
				_matA.resize( _m );
			}
			if( j >= _n ) _n = j + 1;
			setEntry(i,j,val);
		}
		if( ms.getError() > END_OF_MATRIX )
			throw ms.reportError(__func__,__LINE__);
		if( !ms.getDimensions( i, _n ) )
			throw ms.reportError(__func__,__LINE__);
		if( i > _m ) {
			_m = i;
			_matA.resize(_m);
		}
	}

	template <class Field, class Row>
	void SparseMatrixGeneric<Field, Row, VectorCategories::SparseSequenceVectorTag > ::setEntry (size_t i, size_t j, const typename Field::Element &value)
	{
		typedef typename Row::value_type value_type;
		Row &v = _matA[i];
		typename Row::iterator iter;

		if (v.size () == 0) {
			v.push_back ( value_type((unsigned)j, value));
		}
		else {
			iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

			if (iter == v.end () || iter->first != j)
				iter = v.insert (iter, value_type((unsigned)j, value));
			else
				iter->second = value;
		}
	}

	template <class Field, class Row>
	typename Field::Element &SparseMatrixGeneric<Field, Row, VectorCategories::SparseSequenceVectorTag > ::refEntry (size_t i, size_t j)
	{

		Row &v = _matA[i];
		typename Row::iterator iter;

		if (v.size () == 0) {
			v.push_back (std::pair <size_t, Element> (j, field().zero));
			return v.front ().second;
		}
		else {
			iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

			if (iter == v.end () || iter->first != j)
				iter = v.insert (iter, std::pair <size_t, Element> (j, field().zero));

			return iter->second;
		}
	}

	template <class Field, class Row>
	const typename Field::Element &SparseMatrixGeneric<Field, Row, VectorCategories::SparseSequenceVectorTag > ::getEntry (size_t i, size_t j) const
	{

		const Row &v = _matA[i];
		typename Row::const_iterator iter;

		if (v.size () == 0)
			return field().zero;
		else {
			iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

			if (iter == v.end () || iter->first != j)
				return field().zero;
			else
				return iter->second;
		}
	}


	template <class Field, class Row>
	template <class Vector>
	Vector &SparseMatrixGeneric<Field, Row, VectorCategories::SparseSequenceVectorTag >::columnDensity (Vector &v) const
	{
		unsigned int row = 0;

		for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
			typename Row::const_iterator j = i.begin ();

			for (; j != i->begin (); ++j)
				++v[j->first];
		}

		return v;
	}


	template <class Field, class Row>
	SparseMatrixGeneric<Field, Row, VectorCategories::SparseSequenceVectorTag >
	&SparseMatrixGeneric<Field, Row, VectorCategories::SparseSequenceVectorTag >::transpose (SparseMatrixGeneric &AT) const
	{
		unsigned int row = 0;

		for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
			typename Row::const_iterator j = i.begin ();

			for (; j != i->begin (); ++j)
				AT._matA[j->first].push_back (std::pair<size_t, Element> (row, j->second));
		}

		return AT;
	}


} // namespace LinBox
} // namespace Protected

#endif // __LINBOX_matrix_sparse_sequence_INL


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
