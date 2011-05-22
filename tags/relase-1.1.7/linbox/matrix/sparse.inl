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

#ifndef __LINBOX_matrix_sparse_INL
#define __LINBOX_matrix_sparse_INL

#include "linbox/linbox-config.h"

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
	A._A.resize (A._m);//cerr<<A.coldim()<<" "<<A.rowdim()<<endl;
		
    Element x;
	while (is >> i) {
		if (i == 0 || i == (size_t) -1) {is >> j; F.read(is, x); break;}
		is >> j;
		if (i > A._m || j > A._n)
			throw InvalidMatrixInput ();
		F.read (is, x);
		if (! F.isZero(x)) A.setEntry (i - 1, j - 1, x);
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
		do is >> c; while (is && !std::isdigit (c));
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

			if (!F.isZero(a_ij))
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
	::readMagmaCpt (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf)
{
	size_t i, j;
	Element a_ij;
	char c;
	const char matrixstart = '[', matrixend = ']';
	const char rowstart = '[', rowend = ']';
	const char pairstart = '[', pairend = ']';

	A._m = A._n = 0;
	A._A.clear ();

    do {is.get(c);} while (c != matrixstart ); // find matrix start
	i = 0;
	while (true)
	{
        do {is.get(c);} while (c != matrixend && c != rowstart);
		if (c == matrixend) return is;
		else
		{   
			A._m++;
			A._A.push_back (Row ());
	        //processrow(i)
			while (true) 
			{
        		do {is.get(c);} while (c != pairstart && c != rowend ); 
				if (c == rowend) break;
				else
				{  //processpair( j v for row i);
					is >> j; 
					if (j > A._n) A._n = j;
    				do {is.get(c);} while (!std::isdigit(c) && c != '-' && c != '+');
					is.unget();
					F.read(is, a_ij);
			        if (!F.isZero(a_ij)) A.setEntry (i, j-1, a_ij);
					do {is.get(c);} while (c != pairend);
				}
			}
			++i;
		}
	}
	//return is; //BB: unreachable
}

template <class Element, class Row, class Trait>
template <class Field>
std::istream &SparseMatrixReadWriteHelper<Element, Row, Trait>
	::read (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F,
		FileFormatTag format)
{
	char buf[80];
	buf[0]=0;
	char c;

	switch (format) {
	    case FORMAT_DETECT: {
		is.getline (buf, 80);
		std::istringstream str (buf);
		do str >> c; while (isspace (c));

		if (c == '[') {
			if (strchr (buf, ';') != NULL)
				readMatlab (A, is, F, buf);
			else
				readPretty (A, is, F, buf);
		} else if (std::isdigit (c)) {
			do str >> c; while (str && (isspace (c) || std::isdigit (c)));

			if (c == 'M')
				return readGuillaume (A, is, F, buf);
			else
				return readTurner (A, is, F, buf);
		} else
			throw InvalidMatrixInput ();
		break;
		}

	    case FORMAT_TURNER:
		return readTurner (A, is, F, buf);
	    case FORMAT_GUILLAUME:
		return readGuillaume (A, is, F, buf);
	    case FORMAT_MATLAB:
		return readMatlab (A, is, F, buf);
	    case FORMAT_PRETTY:
		return readPretty (A, is, F, buf);
	    case FORMAT_MAGMACPT:
		return readMagmaCpt (A, is, F, buf);
	    default:
	    	throw InvalidMatrixInput();
	}

	return is;
}

template <class Element, class Row, class Trait>
template <class Field>
std::ostream &SparseMatrixWriteHelper<Element, Row, Trait>
	::write (const SparseMatrixBase<Element, Row> &A, std::ostream &os, const Field &F, 
		FileFormatTag format)
{
	typename SparseMatrixBase<Element, Row>::Rep::const_iterator i;
	typename Row::const_iterator j;
	typename Field::Element zero;
	size_t i_idx, j_idx;
	//	int col_width;
	integer c;
        bool firstrow;

	// Avoid massive unneeded overhead in the case that this
	// printing is disabled
	if (commentator.isNullStream (os))
		return os;

	switch (format) {
	    case FORMAT_DETECT:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "format != FORMAT_DETECT");
		//break;//BB: unreachable

	    case FORMAT_TURNER:
		// The i j v triples, with zero based indices.
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j = i->begin (), j_idx = 0; j != i->end (); j++, j_idx++) {
				os << i_idx << ' ' << j->first << ' ';
				F.write (os, j->second);
				os << std::endl;
			}
		}
		break;

	    case FORMAT_ONE_BASED:
		// The i j v triples, with zero based indices.
		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j = i->begin (), j_idx = 0; j != i->end (); j++, j_idx++) {
				os << i_idx + 1 << ' ' << j->first + 1 << ' ';
				F.write (os, j->second);
				os << std::endl;
			}
		}
		break;

	    case FORMAT_GUILLAUME:
		// row col 'M' header line followed by the i j v triples, one based, 
		// followed by 0 0 0.
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

		os << "]";

		break;

	    case FORMAT_MAPLE:
		F.init (zero, 0);

		os << "[";
                firstrow=true;

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); ++i, ++i_idx) {
			if (firstrow) {
                            os << "[";
                            firstrow =false;
                        } else 
                             os << ", [";
                           
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

			os << " ]";
		}

		os << "]";

		break;

	    case FORMAT_PRETTY:
		//F.characteristic (c);
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

	    case FORMAT_MAGMACPT: 
		os << "sparse matrix written in MagmaCpt form is not implemented" << std::endl;
		break;
	}

	return os;
}

template <class Element, class Row>
template <class Field>
std::ostream &SparseMatrixWriteHelper<Element, Row, VectorCategories::SparseParallelVectorTag >
	::write (const SparseMatrixBase<Element, Row> &A, std::ostream &os, const Field &F, 
		FileFormatTag format)
{
	typename SparseMatrixBase<Element, Row>::Rep::const_iterator i;
	typename Row::first_type::const_iterator j_idx;
	typename Row::second_type::const_iterator j_elt;
	typename Field::Element zero;
	size_t i_idx, j_idx_1, col_idx;
	//int col_width;
	integer c;
        bool firstrow;
        
	// Avoid massive unneeded overhead in the case that this
	// printing is disabled
	if (commentator.isNullStream (os))
		return os;

	switch (format) {
	    case FORMAT_DETECT:
		throw PreconditionFailed (__FUNCTION__, __LINE__, "format != FORMAT_DETECT");
		//break//BB: unreachable;

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

	    case FORMAT_ONE_BASED:
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

	    case FORMAT_MAPLE:
		F.init (zero, 0);
                firstrow=true;

		os << "[";

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			if (firstrow) {
                            os << "[";
                            firstrow =false;
                        } else 
                             os << ", [";
                           
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

			os << "]";
		}

		os << "]";

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

		os << "]";

		break;

	    case FORMAT_PRETTY:
		//F.characteristic (c);
		//col_width = (int) ceil (log ((double) c) / M_LN10);
		F.init (zero, 0);

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			os << "  [ ";

			j_idx = i->first.begin ();
			j_elt = i->second.begin ();

			for (col_idx = 0; col_idx < A._n; col_idx++) {
				//os.width (col_width);

				if (j_idx == i->first.end () || col_idx != *j_idx)
					F.write (os, zero);
				else {
					F.write (os, *j_elt);
					++j_idx; ++j_elt;
				}

				os << ' ';
			}

			os << ']';
		}

		break;
	    case FORMAT_MAGMACPT:
		os << "sparse matrix written in MagmaCpt form is not implemented" << std::endl;
		break;
	}

	return os;
}

template <class Element, class Row, class Tag>
template <class Field>
SparseMatrixBase<Element,Row,Tag>
	::SparseMatrixBase( MatrixStream<Field>& ms )
	:_A(0), _m(0), _n(0)
{
	Element val;
	size_t i, j;
	while( ms.nextTriple(i,j,val) ) {
		if( i >= _m ) {
			_m = i + 1;
			_A.resize( _m );
		}
		if( j >= _n ) _n = j + 1;
		setEntry(i,j,val);
	}
	if( ms.getError() > END_OF_MATRIX )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( !ms.getDimensions( i, _n ) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( i > _m ) {
		_m = i;
		_A.resize(_m);
	}
}

template <class Element, class Row>
template <class Field>
SparseMatrixBase<Element,Row,VectorCategories::SparseSequenceVectorTag>
	::SparseMatrixBase( MatrixStream<Field>& ms )
	:_A(0), _m(0), _n(0)
{
	Element val;
	size_t i, j;
	while( ms.nextTriple(i,j,val) ) {
		if( i >= _m ) {
			_m = i + 1;
			_A.resize( _m );
		}
		if( j >= _n ) _n = j + 1;
		setEntry(i,j,val);
	}
	if( ms.getError() > END_OF_MATRIX )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( !ms.getDimensions( i, _n ) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( i > _m ) {
		_m = i;
		_A.resize(_m);
	}
}

template <class Element, class Row>
template <class Field>
SparseMatrixBase<Element,Row,VectorCategories::SparseAssociativeVectorTag>
	::SparseMatrixBase( MatrixStream<Field>& ms )
	:_A(0), _m(0), _n(0)
{
	Element val;
	size_t i, j;
	while( ms.nextTriple(i,j,val) ) {
		if( i >= _m ) {
			_m = i + 1;
			_A.resize( _m );
		}
		if( j >= _n ) _n = j + 1;
		setEntry(i,j,val);
	}
	if( ms.getError() > END_OF_MATRIX )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( !ms.getDimensions( i, _n ) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( i > _m ) {
		_m = i;
		_A.resize(_m);
	}
}

template <class Element, class Row>
template <class Field>
SparseMatrixBase<Element,Row,VectorCategories::SparseParallelVectorTag>
	::SparseMatrixBase( MatrixStream<Field>& ms )
	:_A(0), _m(0), _n(0)
{
	Element val;
	size_t i, j;
	while( ms.nextTriple(i,j,val) ) {
		if( i >= _m ) {
			_m = i + 1;
			_A.resize( _m );
		}
		if( j >= _n ) _n = j + 1;
		setEntry(i,j,val);
	}
	if( ms.getError() > END_OF_MATRIX )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( !ms.getDimensions( i, _n ) )
		throw ms.reportError(__FUNCTION__,__LINE__);
	if( i > _m ) {
		_m = i;
		_A.resize(_m);
	}
}

template <class Element, class Row>
void SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
        typedef typename Row::value_type value_type;
	Row &v = _A[i];
	typename Row::iterator iter;
        
	if (v.size () == 0) {
		v.push_back ( value_type(j, value));                
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

		if (iter == v.end () || iter->first != j)
			iter = v.insert (iter, value_type(j, value));
                else
                    	iter->second = value;
 	}
}

template <class Element, class Row>
Element &SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag >
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

template <class Element, class Row>
const Element &SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag >
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

template <class Element, class Row>
const Element &SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag >
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

template <class Element, class Row>
void SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag >
	::setEntry (size_t i, size_t j, const Element &value) 
{
	while (_A.size() < i + 1) _A.push_back(Row());
	_m = _A.size(); 

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
		} else
                    	*(v.second.begin () + (iter - v.first.begin ())) = value;
	}
}

template <class Element, class Row>
Element &SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag >
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

template <class Element, class Row>
const Element &SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag >
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

template <class Element, class Row>
template <class Vector>
Vector &SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag >::columnDensity (Vector &v) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			++v[j->first];
	}

	return v;
}

template <class Element, class Row>
template <class Vector>
Vector &SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag >::columnDensity (Vector &v) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::first_type::const_iterator j_idx = i->first.begin ();

		for (; j_idx != i->first.end (); ++j_idx)
			++v[*j_idx];
	}

	return v;
}

template <class Element, class Row>
template <class Vector>
Vector &SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag >::columnDensity (Vector &v) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			++v[j->first];
	}

	return v;
}

template <class Element, class Row>
SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag >
	&SparseMatrixBase<Element, Row, VectorCategories::SparseSequenceVectorTag >::transpose (SparseMatrixBase &AT) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			AT._A[j->first].push_back (std::pair<size_t, Element> (row, j->second));
	}

	return AT;
}

template <class Element, class Row>
SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag >
	&SparseMatrixBase<Element, Row, VectorCategories::SparseAssociativeVectorTag >::transpose (SparseMatrixBase &AT) const
{
	unsigned int row = 0;

	for (ConstRowIterator i = rowBegin (); i != rowEnd (); ++i, ++row) {
		typename Row::const_iterator j = i.begin ();

		for (; j != i->begin (); ++j)
			AT._A[j->first][row] = j->second;
	}

	return AT;
}

template <class Element, class Row>
SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag >
	&SparseMatrixBase<Element, Row, VectorCategories::SparseParallelVectorTag >::transpose (SparseMatrixBase &AT) const
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

#endif // __LINBOX_matrix_sparse_INL

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
