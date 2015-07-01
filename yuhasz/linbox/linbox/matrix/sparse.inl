/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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
 * See COPYING for license information.
 */

#ifndef __MATRIX_SPARSE_INL
#define __MATRIX_SPARSE_INL

#include "linbox-config.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>

#include "linbox/matrix/sparse.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/util/debug.h"
#include <linbox/util/commentator.h>

namespace LinBox
{

#ifndef __LINBOX_XMLENABLED

template <class Element, class Row, class Trait>
template <class Field>
std::istream &SparseMatrixReadWriteHelper<Element, Row, Trait>
	::readTurner (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf)
{
	size_t i, j;

	A._A.clear ();
	A._A.resize (A._m);

	do {
		std::istringstream str (buf);

		str >> i;

		if (i == (size_t) -1) break; // return also if row index is -1
		str >> j;
		F.read (str, A.refEntry (i, j));

		is.getline (buf, 80);
	} while (is);

	return is;
}

template <class Element, class Row, class Trait>
template <class Field>
std::istream &SparseMatrixReadWriteHelper<Element, Row, Trait>
	::readGuillaume (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf)
{
	size_t i, j;

	std::istringstream str (buf);

	str >> A._m >> A._n;

	A._A.clear ();
	A._A.resize (A._m);

	while (is >> i) {
		if (i == 0 || i == (size_t) -1) break;
		is >> j;
		if (i > A._m || j > A._n)
			throw InvalidMatrixInput ();
		F.read (is, A.refEntry (i - 1, j - 1));
	}

	return is;

}

template <class Element, class Row, class Trait>
template <class Field>
std::istream &SparseMatrixReadWriteHelper<Element, Row, Trait>
	::readMatlab (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf)
{
	size_t i = 0, j = 0;
	char c;
	Element a_ij;

	while (1) {
		do is >> c; while (is && !isdigit (c));
		if (!is) break;

		is.putback (c);

		F.read (is, a_ij);
		A.setEntry (i, j++, a_ij);

		do is >> c; while (is && c != ',' && c != ';' && c != ']');
		if (!is) break;;

		if (c == ';') {
			++i;
			j = 0;
		}
		else if (c == ']') break;
	}

	return is;
}

template <class Element, class Row, class Trait>
template <class Field>
std::istream &SparseMatrixReadWriteHelper<Element, Row, Trait>
	::readPretty (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf)
{
	size_t i, j;
	Element a_ij;
	char c;

	A._m = 0;
	A._A.clear ();

	i = 0;

	do {
		A._m++;
		A._A.push_back (Row ());

		std::istringstream str (buf);

		do str >> c; while (isspace (c));
		if (c != '[')
			throw InvalidMatrixInput ();

		j = 0;

		while (str) {
			do str >> c; while (isspace (c));
			if (!str || c == ']') break;
			F.read (str, a_ij);

			j++;
			if (j > A._n)
				A._n++;

			if (a_ij != 0)
				A.setEntry (i, j, a_ij);
		}

		is.getline (buf, 80);

		i++;
	} while (is);

	return is;

}

template <class Element, class Row, class Trait>
template <class Field>
std::istream &SparseMatrixReadWriteHelper<Element, Row, Trait>
	::read (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F,
		typename SparseMatrixWriteHelper<Element, Row, Trait>::Format format)
{
	char buf[80];
	char c;

	is.getline (buf, 80);
	std::istringstream str (buf);

	switch (format) {
	    case FORMAT_DETECT:
		do str >> c; while (isspace (c));

		if (c == '[') {
			if (strchr (buf, ';') != NULL)
				readMatlab (A, is, F, buf);
			else
				readPretty (A, is, F, buf);
		} else if (isdigit (c)) {
			do str >> c; while (str && (isspace (c) || isdigit (c)));

			if (c == 'M')
				return readGuillaume (A, is, F, buf);
			else
				return readTurner (A, is, F, buf);
		} else
			throw InvalidMatrixInput ();

	    case FORMAT_TURNER:
		return readTurner (A, is, F, buf);
	    case FORMAT_GUILLAUME:
		return readGuillaume (A, is, F, buf);
	    case FORMAT_MATLAB:
		return readMatlab (A, is, F, buf);
	    case FORMAT_PRETTY:
		return readPretty (A, is, F, buf);
	}

	return is;
}

template <class Element, class Row, class Trait>
template <class Field>
std::ostream &SparseMatrixWriteHelper<Element, Row, Trait>
	::write (const SparseMatrixBase<Element, Row> &A, std::ostream &os, const Field &F, Format format)
{
	typename SparseMatrixBase<Element, Row>::Rep::const_iterator i;
	typename Row::const_iterator j;
	typename Field::Element zero;
	size_t i_idx, j_idx;
	//int col_width;
	integer c;

	// Avoid massive unneeded overhead in the case that this
	// printing is disabled
	if (commentator.isNullStream (os))
		return os;

	switch (format) {
	    case FORMAT_DETECT:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "format != FORMAT_DETECT");
		break;

	    case FORMAT_TURNER:
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j = i->begin (), j_idx = 0; j != i->end (); j++, j_idx++) {
				os << i_idx << ' ' << j->first << ' ';
				F.write (os, j->second);
				os << std::endl;
			}
		}

		break;

	    case FORMAT_GUILLAUME:
		os << A._m << ' ' << A._n << " M" << std::endl;

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j = i->begin (), j_idx = 0; j != i->end (); j++, j_idx++) {
				os << i_idx + 1 << ' ' << j->first + 1 << ' ';
				F.write (os, j->second);
				os << std::endl;
			}
		}

		os << "0 0 0" << std::endl;

		break;

	    case FORMAT_MATLAB:
		F.init (zero, 0);

		os << "[";

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			j = i->begin ();

			for (j_idx = 0; j_idx < A._n; j_idx++) {
				if (j == i->end () || j_idx != j->first)
					F.write (os, zero);
				else {
					F.write (os, j->second);
					j++;
				}

				if (j_idx < A._n - 1)
					os << ", ";
			}

			os << "; ";
		}

		os << "]" << std::endl;

		break;

	    case FORMAT_PRETTY:
		F.characteristic (c);
		//col_width = (int) ceil (log ((double) c) / M_LN10);
		F.init (zero, 0);

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			os << "  [ ";

			j = i->begin ();

			for (j_idx = 0; j_idx < A._n; j_idx++) {
				//os.width (col_width);

				if (j == i->end () || j_idx != j->first)
					F.write (os, zero);
				else {
					F.write (os, j->second);
					j++;
				}

				os << ' ';
			}

			os << ']' << std::endl;
		}

		break;
	}

	return os;
}

template <class Element, class Row, class RowTrait>
template <class Field>
std::ostream &SparseMatrixWriteHelper<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >
	::write (const SparseMatrixBase<Element, Row> &A, std::ostream &os, const Field &F, Format format)
{
	typename SparseMatrixBase<Element, Row>::Rep::const_iterator i;
	typename Row::first_type::const_iterator j_idx;
	typename Row::second_type::const_iterator j_elt;
	typename Field::Element zero;
	size_t i_idx, j_idx_1, col_idx;
	//int col_width;
	integer c;

	// Avoid massive unneeded overhead in the case that this
	// printing is disabled
	if (commentator.isNullStream (os))
		return os;

	switch (format) {
	    case FORMAT_DETECT:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "format != FORMAT_DETECT");
		break;

	    case FORMAT_TURNER:
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j_idx = i->first.begin (), j_elt = i->second.begin ();
			     j_idx != i->first.end ();
			     ++j_idx, ++j_elt)
			{
				os << i_idx << ' ' << *j_idx << ' ';
				F.write (os, *j_elt);
				os << std::endl;
			}
		}

		break;

	    case FORMAT_GUILLAUME:
		os << A._m << ' ' << A._n << " M" << std::endl;

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j_idx = i->first.begin (), j_elt = i->second.begin ();
			     j_idx != i->first.end ();
			     ++j_idx, ++j_elt)
			{
				os << i_idx + 1 << ' ' << *j_idx + 1 << ' ';
				F.write (os, *j_elt);
				os << std::endl;
			}
		}

		os << "0 0 0" << std::endl;

		break;

	    case FORMAT_MATLAB:
		F.init (zero, 0);

		os << "[";

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			j_idx = i->first.begin ();
			j_elt = i->second.begin ();

			for (j_idx_1 = 0; j_idx_1 < A._n; j_idx_1++) {
				if (j_idx == i->first.end () || j_idx_1 != *j_idx)
					F.write (os, zero);
				else {
					F.write (os, *j_elt);
					++j_idx;
					++j_elt;
				}

				if (j_idx_1 < A._n - 1)
					os << ", ";
			}

			os << "; ";
		}

		os << "]" << std::endl;

		break;

	    case FORMAT_PRETTY:
		F.characteristic (c);
		//col_width = (int) ceil (log ((double) c) / M_LN10);
		F.init (zero, 0);

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			os << "  [ ";

			j_idx = i->first.begin ();
			j_elt = i->second.begin ();

			for (col_idx = 0; col_idx < A._n; col_idx++) {
				//nos.width (col_width);

				if (j_idx == i->first.end () || col_idx != *j_idx)
					F.write (os, zero);
				else {
					F.write (os, *j_elt);
					++j_idx; ++j_elt;
				}

				os << ' ';
			}

			os << ']' << std::endl;
		}

		break;
	}

	return os;
}

#else

template<class Element, class Row, class RowTrait>
ostream &SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >::write(ostream &out) const {
	Writer W;
	if( toTag(W)) 
		W.write(out);
		
	return out;
}

template<class Element, class Row, class RowTrait>
ostream &SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >::write(ostream &out) const 
{
	Writer W;
	if( toTag(W) ) 
		W.write(out);
	return out;
}

template<class Element, class Row, class RowTrait>
ostream &SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag<RowTrait> >::write(ostream &out) const
{
	Writer W;
	if( toTag(W) ) 
		W.write(out);
	return out;
}

template<class Element, class Row, class RowTrait>
bool SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >::toTag(Writer &W) const
{

	vector<Element> elem;
	vector<size_t> row, col;
	size_t i;
	string holder;
	typename Row::const_iterator iter;

	W.setTagName("MatrixOver");
	W.setAttribute("rows", Writer::numToString(holder, _m));
	W.setAttribute("cols", Writer::numToString(holder, _n));
	W.setAttribute("implDetail", "sparse-sequence");

	W.addTagChild();
	W.setTagName("sparseMatrix");
	
	for(i = 0; i < _m; ++i) {
		
		for(iter = _A[i].begin(); iter != _A[i].end(); ++iter) {
			row.push_back(i);
			col.push_back(iter->first);
			elem.push_back(iter->second);
		}
	}

	W.addTagChild();
	W.setTagName("index");
	W.addNumericalList(row);
	W.upToParent();
	
	W.addTagChild();
	W.setTagName("index");
	W.addNumericalList(col);
	W.upToParent();

	W.addTagChild();
	W.setTagName("entry");
	W.addNumericalList(elem);
	W.upToParent();

	W.upToParent();

	return true;
}


template<class Element, class Row, class RowTrait>
SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >::SparseMatrixBase(Reader &R)
{

	Element e;
	vector<Element> elem;
	typename vector<Element>::const_iterator iter;
	vector<size_t> row, col;
	vector<size_t>::const_iterator i1, i2;
	typename Row::iterator ii;
	size_t i;

	if(!R.expectTagName("MatrixOver")) return;
	if(!R.expectAttributeNum("rows", _m) || !R.expectAttributeNum("cols", _n)) return;

	if(!R.expectChildTag()) return;

	R.traverseChild();
	if(R.checkTagName("field")) { // skip the field if there is one
		R.upToParent();
		if(!R.getNextChild()) {
			R.setErrorString("Got a matrix that just had a field!");
			R.setErrorCode(Reader::OTHER);
			return;
		}

		if(!R.expectChildTag()) return;
		R.traverseChild();
	}

	if(R.checkTagName("diag")) {
		if(!R.expectChildTag()) return;
		R.traverseChild();

		if(!R.expectTagName("entry") || !R.expectTagNumVector(elem)) return;

		_A.resize(_m);

		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		for(i = 0, iter = elem.begin(); i < _m && i < _n; ++i, ++iter) {
			_A[i].push_back(std::pair<size_t, Element>(i, *iter));
		}
	}
	else if(R.checkTagName("scalar")) {

		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagNum(e)) return;
		R.upToParent();

		R.upToParent();
		R.getPrevChild();

		_A.resize(_m);
		for(i = 0; i < _m && i < _n; ++i) {
			_A[i].push_back(std::pair<size_t, Element>(i, e));
		}
	}
	else if(R.checkTagName("zero-one")) {
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(row)) return;
		R.upToParent();

		if(!R.getNextChild()) {
			R.setErrorString("Could not find column indices for zero-one matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}

		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") ||  !R.expectTagNumVector(col)) return;
		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		e = integer(1);
		_A.resize(_m);
		for(i1 = row.begin(), i2 = col.begin(); i1 != row.end(); ++i1, ++i2) {
			
			ii = std::lower_bound(_A[*i1].begin(), _A[*i1].end(), *i2, VectorWrapper::CompareSparseEntries<Element>());
			if(ii == _A[*i1].end() || ii->first != *i2) {
				_A[*i1].insert(ii, std::pair<size_t, Element>(*i2, e));
			}
		}
		
	}
	else if(!R.expectTagName("sparseMatrix")) 
		return;
	else {

		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(row)) return;
		R.upToParent();

		if(!R.getNextChild()) {
			R.setErrorString("Could not find columnar indices for sparse matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}

		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(col)) return;
		R.upToParent();

		if(!R.getNextChild()) {
			R.setErrorString("Could not find entries for sparse matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("entry") || !R.expectTagNumVector(elem)) return;
		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		_A.resize(_m);
		for(i1 = row.begin(), i2 = col.begin(), iter = elem.begin(); i1 != row.end(); ++i1, ++i2, ++iter) {
			
			ii = std::lower_bound(_A[*i1].begin(), _A[*i1].end(), *i2, VectorWrapper::CompareSparseEntries<Element>());
			if(ii == _A[*i1].end() || ii->first != *i2) {
				_A[*i1].insert(ii, std::pair<size_t, Element>(*i2, *iter));
			}
		}
	}

	return;
}


template<class Element, class Row, class RowTrait>
bool SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag<RowTrait> >::toTag(Writer &W) const
{

	vector<Element> elem;
	vector<size_t> row, col;
	size_t i;
	string holder;
	typename Row::const_iterator iter;

	W.setTagName("MatrixOver");
	W.setAttribute("rows", Writer::numToString(holder, _m));
	W.setAttribute("cols", Writer::numToString(holder, _n));
	W.setAttribute("implDetail", "sparse-associative");

	W.addTagChild();
	W.setTagName("sparseMatrix");
	
	for(i = 0; i < _m; ++i) {
		
		for(iter = _A[i].begin(); iter != _A[i].end(); ++iter) {
			row.push_back(i);
			col.push_back(iter->first);
			elem.push_back(iter->second);
		}
	}

	W.addTagChild();
	W.setTagName("index");
	W.addNumericalList(row);
	W.upToParent();
	
	W.addTagChild();
	W.setTagName("index");
	W.addNumericalList(col);
	W.upToParent();

	W.addTagChild();
	W.setTagName("entry");
	W.addNumericalList(elem);
	W.upToParent();

	W.upToParent();

	return true;
}


template<class Element, class Row, class RowTrait>
SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag<RowTrait> >::SparseMatrixBase(Reader &R)
{

	Element e;
	vector<Element> elem;
	typename vector<Element>::const_iterator iter;
	vector<size_t> row, col;
	vector<size_t>::const_iterator i1, i2;
	size_t i;

	if(!R.expectTagName("MatrixOver")) 
		return;

	

	if(!R.expectAttributeNum("rows", _m) || !R.expectAttributeNum("cols", _n)) 
		return;


	if(!R.expectChildTag()) return;

	R.traverseChild();
	if(R.checkTagName("field")) { // skip the field if there is one
		R.upToParent();
		if(!R.getNextChild()) {
			R.setErrorString("Got a matrix with a field and no data.");
			R.setErrorCode(Reader::OTHER);
			return;
		}

		if(!R.expectChildTag()) return;
		R.traverseChild();
	}

	if(R.checkTagName("diag")) {
		if(!R.expectChildTag()) return;
		R.traverseChild();

		if(!R.expectTagName("entry") || !R.expectTagNumVector(elem)) return;

		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		_A.resize(_m);

		for(i = 0, iter = elem.begin(); i < _m && i < _n; ++i, ++iter) {
			_A[i].insert(std::pair<size_t, Element>(i, *iter));
		}
	}
	else if(R.checkTagName("scalar")) {

		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagNum(e)) return;
		R.upToParent();

		R.upToParent();
		R.getPrevChild();

		_A.resize(_m);
		for(i = 0; i < _m && i < _n; ++i) {
			_A[i].insert(std::pair<size_t, Element>(i, e));
		}
	}
	else if(R.checkTagName("zero-one")) {
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(row)) return;
		R.upToParent();

		if(!R.getNextChild()) {
			R.setErrorString("Could not find columnar indices for zero-one matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		if( !R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(col)) return;

		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		e = integer(1);
		_A.resize(_m);
		for(i1 = row.begin(), i2 = col.begin(); i1 != row.end(); ++i1, ++i2) {
			_A[*i1].insert(std::pair<size_t, Element>(*i2, e));
		}
	}
	else if(!R.expectTagName("sparseMatrix")) 
		return;
	else {
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(row)) return;
		R.upToParent();


		if(!R.getNextChild()) {
			R.setErrorString("Could not find columnar indices for sparse matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(col)) return;
		R.upToParent();

		if(!R.getNextChild()) {
			R.setErrorString("Could not find entries for sparse matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("entry") || !R.expectTagNumVector(elem)) return;
		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		_A.resize(_m);
		for(i1 = row.begin(), i2 = col.begin(), iter = elem.begin(); i1 != row.end(); ++i1, ++i2, ++iter) {
			_A[*i1].insert(std::pair<size_t, Element>(*i2, *iter));
		}

	}

	return;
}


template<class Element, class Row, class RowTrait>
bool SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >::toTag(Writer &W) const
{

	vector<Element> elem;
	vector<size_t> row, col;
	size_t i;
	string holder;
	typename Row::first_type::const_iterator iter1;
	typename Row::second_type::const_iterator iter2;

	W.setTagName("MatrixOver");
	W.setAttribute("rows", Writer::numToString(holder, _m));
	W.setAttribute("cols", Writer::numToString(holder, _n));
	W.setAttribute("implDetail", "sparse-parallel");

	W.addTagChild();
	W.setTagName("sparseMatrix");
	
	for(i = 0; i < _m; ++i) {
		
		for(iter1 = _A[i].first.begin(), iter2 = _A[i].second.begin(); iter1 != _A[i].first.end(); ++iter1, ++iter2) {
			row.push_back(i);
			col.push_back(*iter1);
			elem.push_back(*iter2);
		}
	}

	W.addTagChild();
	W.setTagName("index");
	W.addNumericalList(row);
	W.upToParent();
	
	W.addTagChild();
	W.setTagName("index");
	W.addNumericalList(col);
	W.upToParent();

	W.addTagChild();
	W.setTagName("entry");
	W.addNumericalList(elem);
	W.upToParent();

	W.upToParent();

	return true;
}


template<class Element, class Row, class RowTrait>
SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >::SparseMatrixBase(Reader &R)
{

	Element e;
	vector<Element> elem;
	typename vector<Element>::const_iterator iter;
	vector<size_t> row, col;
	vector<size_t>::const_iterator i1, i2;
	typename Row::first_type::iterator fi;
	size_t i;

	if(!R.expectTagName("MatrixOver")) return;
	if(!R.expectAttributeNum("rows", _m) || !R.expectAttributeNum("cols", _n)) return;

	if(!R.expectChildTag()) return;

	R.traverseChild();
	if(R.checkTagName("field")) { // skip the field if there is one
		R.upToParent();
		if(!R.getNextChild()) {
			R.setErrorString("Got a matrix with a field and no data.");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		if(!R.expectChildTag()) return;
		R.traverseChild();
	}

	if(R.checkTagName("diag")) {
		if(!R.expectChildTag()) return;
		R.traverseChild();

		if(!R.expectTagName("entry") || !R.expectTagNumVector(elem)) return;
		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		_A.resize(_m);

		for(i = 0, iter = elem.begin(); i < _m && i < _n; ++i, ++iter) {
			_A[i].first.push_back(i);
			_A[i].second.push_back(*iter);
		}
	}
	else if(R.checkTagName("scalar")) {

		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagNum(e)) return;
		R.upToParent();

		R.upToParent();
		R.getPrevChild();

		_A.resize(_m);
		for(i = 0; i < _m && i < _n; ++i) {
			_A[i].first.push_back(i);
			_A[i].second.push_back(e);
		}
	}
	else if(R.checkTagName("zero-one")) {
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(row)) return;
		R.upToParent();

		if(!R.getNextChild()) {
			R.setErrorString("Could not find columnar indices for zero one matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(col)) return;
		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		e = integer(1);
		_A.resize(_m);
		for(i1 = row.begin(), i2 = col.begin(); i1 != row.end(); ++i1, ++i2) {
			fi = std::lower_bound (_A[*i1].first.begin (), _A[*i1].first.end (), *i2);

			if (fi == _A[*i1].first.end () || *fi != *i2) {
				fi = _A[*i1].first.insert (fi, *i2);
				_A[*i1].second.insert (_A[*i1].second.begin () + (fi - _A[*i1].first.begin ()), e);
			}
		}
	}

	else if(!R.expectTagName("sparseMatrix")) 
		return;
	else {
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(row)) return;
		R.upToParent();

		if(!R.getNextChild()) {
			R.setErrorString("Could not find columnar indices for sparse matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("index") || !R.expectTagNumVector(col)) return;
		R.upToParent();

		if(!R.getNextChild()) {
			R.setErrorString("Could not find entries of sparse matrix");
			R.setErrorCode(Reader::OTHER);
			return;
		}
		if(!R.expectChildTag()) return;
		R.traverseChild();
		if(!R.expectTagName("entry") || !R.expectTagNumVector(elem)) return;
		R.upToParent();
		R.upToParent();
		R.getPrevChild();

		_A.resize(_m);


		
		for(i1 = row.begin(), i2 = col.begin(), iter = elem.begin(); i1 != row.end(); ++i1, ++i2, ++iter) {
			fi = std::lower_bound (_A[*i1].first.begin (), _A[*i1].first.end (), *i2);

			if (fi == _A[*i1].first.end () || *fi != *i2) {
				fi = _A[*i1].first.insert (fi, *i2);
				_A[*i1].second.insert (_A[*i1].second.begin () + (fi - _A[*i1].first.begin ()), *iter);
			}
		}
		
	}
	return;
}

#endif



template <class Element, class Row, class RowTrait>
void SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	Row &v = _A[i];
	typename Row::iterator iter;

	if (v.size () == 0) {
		v.push_back (std::pair <size_t, Element> (j, value));
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

		if (iter == v.end () || iter->first != j)
			iter = v.insert (iter, std::pair <size_t, Element> (j, value));
	}
}

template <class Element, class Row, class RowTrait>
Element &SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >
	::refEntry (size_t i, size_t j) 
{
	static Element zero;

	Row &v = _A[i];
	typename Row::iterator iter;

	if (v.size () == 0) {
		v.push_back (std::pair <size_t, Element> (j, zero));
		return v.front ().second;
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

		if (iter == v.end () || iter->first != j)
			iter = v.insert (iter, std::pair <size_t, Element> (j, zero));

		return iter->second;
	}
}

template <class Element, class Row, class RowTrait>
const Element &SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >
	::getEntry (size_t i, size_t j) const
{
	static Element zero;

	const Row &v = _A[i];
	typename Row::const_iterator iter;

	if (v.size () == 0)
		return zero;
	else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

		if (iter == v.end () || iter->first != j)
			return zero;
		else
			return iter->second;
	}
}

template <class Element, class Row, class RowTrait>
const Element &SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag<RowTrait> >
	::getEntry (size_t i, size_t j) const
{
	static Element zero;

	const Row &v = _A[i];
	typename Row::const_iterator iter;

	if (v.size () == 0)
		return zero;
	else {
		iter = v.find (j);

		if (iter == v.end () || iter->first != j)
			return zero;
		else
			return iter->second;
	}
}

template <class Element, class Row, class RowTrait>
void SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	Row &v = _A[i];
	typename Row::first_type::iterator iter;

	if (v.first.size () == 0) {
		v.first.push_back (j);
		v.second.push_back (value);
	} else {
		iter = std::lower_bound (v.first.begin (), v.first.end (), j);

		if (iter == v.first.end () || *iter != j) {
			iter = v.first.insert (iter, j);
			v.second.insert (v.second.begin () + (iter - v.first.begin ()), value);
		}
	}
}

template <class Element, class Row, class RowTrait>
Element &SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >
	::refEntry (size_t i, size_t j) 
{
	static Element zero;

	Row &v = _A[i];
	typename Row::first_type::iterator iter;
	typename Row::second_type::iterator iter_elt;

	if (v.first.size () == 0) {
		v.first.push_back (j);
		v.second.push_back (zero);
		return v.second.front ();
	} else {
		iter = std::lower_bound (v.first.begin (), v.first.end (), j);

		if (iter == v.first.end () || *iter != j) {
			iter = v.first.insert (iter, j);
			iter_elt = v.second.insert (v.second.begin () + (iter - v.first.begin ()), zero);
		}
		else
			iter_elt = v.second.begin () + (iter - v.first.begin ());

		return *iter_elt;
	}
}

template <class Element, class Row, class RowTrait>
const Element &SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >
	::getEntry (size_t i, size_t j) const
{
	static Element zero;

	const Row &v = _A[i];
	typename Row::first_type::const_iterator iter;

	if (v.first.size () == 0)
		return zero;
	else {
		iter = std::lower_bound (v.first.begin (), v.first.end (), j);

		if (iter == v.first.end () || *iter != j)
			return zero;
		else
			return *(v.second.begin () + (iter - v.first.begin ()));
	}
}

template <class Element, class Row, class RowTrait>
template <class Vector>
Vector &SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >::columnDensity (Vector &v) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			++v[j->first];
	}

	return v;
}

template <class Element, class Row, class RowTrait>
template <class Vector>
Vector &SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >::columnDensity (Vector &v) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::first_type::const_iterator j_idx = i->first.begin ();

		for (; j_idx != i->first.end (); ++j_idx)
			++v[*j_idx];
	}

	return v;
}

template <class Element, class Row, class RowTrait>
template <class Vector>
Vector &SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag<RowTrait> >::columnDensity (Vector &v) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			++v[j->first];
	}

	return AT;
}

template <class Element, class Row, class RowTrait>
SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >
	&SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag<RowTrait> >::transpose (SparseMatrixBase &AT) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			AT._A[j->first].push_back (std::pair<size_t, Element> (row, j->second));
	}

	return AT;
}

template <class Element, class Row, class RowTrait>
SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag<RowTrait> >
	&SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag<RowTrait> >::transpose (SparseMatrixBase &AT) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			AT._A[j->first][row] = j->second;
	}

	return AT;
}

template <class Element, class Row, class RowTrait>
SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >
	&SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag<RowTrait> >::transpose (SparseMatrixBase &AT) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::first_type::const_iterator j_idx = i->first.begin ();
		typename Row::second_type::const_iterator j_elt = i->second.begin ();

		for (; j_idx != i->first.end (); ++j_idx, ++j_elt) {
			AT._A[*j_idx].first.push_back (row);
			AT._A[*j_idx].second.push_back (*j_elt);
		}
	}

	return AT;
}

} // namespace LinBox

#endif // __MATRIX_SPARSE_INL
