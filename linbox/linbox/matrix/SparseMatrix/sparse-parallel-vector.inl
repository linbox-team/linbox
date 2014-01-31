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

#ifndef __LINBOX_matrix_sparse_parallel_INL
#define __LINBOX_matrix_sparse_parallel_INL

namespace LinBox {

	template <class Field, class Row>
	std::ostream &SparseMatrixWriteHelper<Protected::SparseMatrixGeneric<Field, Row,
	VectorCategories::SparseParallelVectorTag > >::writeTriple (const Protected::SparseMatrixGeneric<Field, Row> &A
								    , std::ostream &os
								    , bool oneBased)
	{
		typename Matrix::Rep::const_iterator i;
		typename Row::first_type::const_iterator j_idx;
		typename Row::second_type::const_iterator j_elt;
		const Field & F = A.field();

		size_t i_idx;
		for (i = A.getRep().begin (), i_idx = 0; i != A.getRep().end (); i++, i_idx++) {
			for (j_idx = i->first.begin (), j_elt = i->second.begin ();
			     j_idx != i->first.end ();
			     ++j_idx, ++j_elt)
			{
				if (oneBased)
					os << i_idx +1 << ' ' << *j_idx +1 << ' ';
				else
					os << i_idx << ' ' << *j_idx << ' ';
				F.write (os, *j_elt);
				os << std::endl;
			}
		}
		return os ;


	}

	template <class Field, class Row>
	std::ostream &SparseMatrixWriteHelper<Protected::SparseMatrixGeneric<Field, Row,
	VectorCategories::SparseParallelVectorTag > >::writePretty (const Protected::SparseMatrixGeneric<Field, Row> &A
								    , std::ostream &os
								    , std::string begmat
								    , std::string endmat
								    , std::string begrow
								    , std::string endrow
								    , std::string sepelt
								    , std::string seprow
								   )
	{
		typename Matrix::Rep::const_iterator i;
		typename Row::first_type::const_iterator j_idx;
		typename Row::second_type::const_iterator j_elt;
		const Field & F = A.field();

		size_t i_idx, j_idx_1;
		bool firstrow=true;

		os << begmat;

		for (i = A.getRep().begin (), i_idx = 0; i != A.getRep().end (); i++, i_idx++) {
			if (firstrow) {
				os << begrow;
				firstrow =false;
			}
			else
				os << seprow << begrow ;

			j_idx = i->first.begin ();
			j_elt = i->second.begin ();

			for (j_idx_1 = 0; j_idx_1 < A.coldim(); j_idx_1++) {
				if (j_idx == i->first.end () || j_idx_1 != *j_idx)
					F.write (os, F.zero);
				else {
					F.write (os, *j_elt);
					++j_idx;
					++j_elt;
				}

				if (j_idx_1 < A.coldim() - 1)
					os << sepelt << ' ';
			}

			os << endrow;
		}

		os << endmat;

		return os ;
	}


	//! @bug this is reall the "generic" one
	template <class Field, class Row>
	std::ostream &SparseMatrixWriteHelper<Protected::SparseMatrixGeneric<Field, Row,
	VectorCategories::SparseParallelVectorTag > >::write (const Protected::SparseMatrixGeneric<Field, Row> &A
							      , std::ostream &os
							      , LINBOX_enum(Tag::FileFormat) format)
	{

		// Avoid massive unneeded overhead in the case that this
		// printing is disabled
		if (not os)
			return os;

		switch (format) {
		case Tag::FileFormat::Detect:
			throw PreconditionFailed (__func__, __LINE__, "format != Tag::FileFormat::Detect");
			//break//BB: unreachable;

		case Tag::FileFormat::Turner:
			return writeTriple(A,os,false);


		case Tag::FileFormat::OneBased:
			return writeTriple(A,os,true);

		case Tag::FileFormat::Guillaume:
			os << A.rowdim() << ' ' << A.coldim() << " M" << std::endl;

			writeTriple(A,os,true);

			os << "0 0 0" << std::endl;

			return os;

		case Tag::FileFormat::Matlab:
			return writePretty(A,os,"[","]","","; ",",","");
			// std::string begmat = "[";
			// std::string endmat = "]";
			// std::string begrow = "";
			// std::string endrow = "; ";
			// std::string sepelt  = ",";
			// std::string seprow  = "";

		case Tag::FileFormat::Maple:
			return writePretty(A,os,"[","]","["," ]",",",", ");
			// std::string begmat = "[";
			// std::string endmat = "]";
			// std::string begrow = "[";
			// std::string endrow = " ]";
			// std::string sepelt = ",";
			// std::string seprow = ", "

		case Tag::FileFormat::Pretty:
			return writePretty(A,os,"",""," [ ","]\n"," ","");
			// std::string begmat = "";
			// std::string endmat = "";
			// std::string begrow = " [ ";
			// std::string endrow = "]\n";
			// std::string sepelt  = " ";
			// std::string seprow  = "";

		case Tag::FileFormat::MatrixMarket:
			writeMMCoordHeader(os, A, A.size(), "SparseMatrixGeneric");
			return writeTriple(A,os,true);


		case Tag::FileFormat::MagmaCpt:
			os << "sparse matrix written in MagmaCpt form is not implemented" << std::endl;
			break;
		default:
			os << "sparse matrix written in Tag::FileFormat::" << (int)format << " is not implemented" << std::endl;
		}

		return os;
	}

} // LinBox

namespace LinBox { namespace Protected {

	template <class Field, class Row>
	SparseMatrixGeneric<Field,Row,VectorCategories::SparseParallelVectorTag> ::SparseMatrixGeneric( MatrixStream<Field>& ms ) :
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
	void SparseMatrixGeneric<Field, Row, VectorCategories::SparseParallelVectorTag > ::setEntry (size_t i, size_t j, const typename Field::Element &value)
	{
		while (_matA.size() < i + 1) _matA.push_back(Row());
		_m = _matA.size();

		Row &v = _matA[i];
		typename Row::first_type::iterator iter;

		if (v.first.size () == 0) {
			v.first.push_back (j);
			v.second.push_back (value);
		}
		else {
			iter = std::lower_bound (v.first.begin (), v.first.end (), j);

			if (iter == v.first.end () || *iter != j) {
				iter = v.first.insert (iter, j);
				v.second.insert (v.second.begin () + (iter - v.first.begin ()), value);
			}
			else
				*(v.second.begin () + (iter - v.first.begin ())) = value;
		}
	}

	template <class Field, class Row>
	typename Field::Element &SparseMatrixGeneric<Field, Row, VectorCategories::SparseParallelVectorTag > ::refEntry (size_t i, size_t j)
	{

		Row &v = _matA[i];
		typename Row::first_type::iterator iter;
		typename Row::second_type::iterator iter_elt;

		if (v.first.size () == 0) {
			v.first.push_back (j);
			v.second.push_back (field().zero);
			return v.second.front ();
		}
		else {
			iter = std::lower_bound (v.first.begin (), v.first.end (), j);

			if (iter == v.first.end () || *iter != j) {
				iter = v.first.insert (iter, j);
				iter_elt = v.second.insert (v.second.begin () + (iter - v.first.begin ()), field().zero);
			}
			else
				iter_elt = v.second.begin () + (iter - v.first.begin ());

			return *iter_elt;
		}
	}

	template <class Field, class Row>
	const typename Field::Element &SparseMatrixGeneric<Field, Row, VectorCategories::SparseParallelVectorTag > ::getEntry (size_t i, size_t j) const
	{

		const Row &v = _matA[i];
		typename Row::first_type::const_iterator iter;

		if (v.first.size () == 0)
			return field().zero;
		else {
			iter = std::lower_bound (v.first.begin (), v.first.end (), j);

			if (iter == v.first.end () || *iter != j)
				return field().zero;
			else
				return *(v.second.begin () + (iter - v.first.begin ()));
		}
	}


	template <class Field, class Row>
	template <class Vector>
	Vector &SparseMatrixGeneric<Field, Row, VectorCategories::SparseParallelVectorTag >::columnDensity (Vector &v) const
	{
		unsigned int row = 0;

		for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
			typename Row::first_type::const_iterator j_idx = i->first.begin ();

			for (; j_idx != i->first.end (); ++j_idx)
				++v[*j_idx];
		}

		return v;
	}



	template <class Field, class Row>
	SparseMatrixGeneric<Field, Row, VectorCategories::SparseParallelVectorTag >
	&SparseMatrixGeneric<Field, Row, VectorCategories::SparseParallelVectorTag >::transpose (SparseMatrixGeneric &AT) const
	{
		unsigned int row = 0;

		for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
			typename Row::first_type::const_iterator j_idx = i->first.begin ();
			typename Row::second_type::const_iterator j_elt = i->second.begin ();

			for (; j_idx != i->first.end (); ++j_idx, ++j_elt) {
				AT._matA[*j_idx].first.push_back (row);
				AT._matA[*j_idx].second.push_back (*j_elt);
			}
		}

		return AT;
	}


} // namespace LinBox
} // namespace Protected

#endif // __LINBOX_matrix_sparse_parallel_INL


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
