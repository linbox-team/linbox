/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/sparse-base.inl
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               1999-2001 William J Turner,
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 * Based on sparse-base.h by William J Turner <wjturner@math.ncsu.edu>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __SPARSE_BASE_INL
#define __SPARSE_BASE_INL

#include "linbox-config.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstring>

#include "linbox/blackbox/sparse-base.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/field/vector-domain.h"
#include "linbox/util/debug.h"

namespace LinBox
{

// Versions using field read

template <class Element, class Row>
template <class Field>
istream &SparseMatrix0ReadWriteHelper<Element, Row>
	::readTurner (SparseMatrix0Base<Element, Row> &A, istream &is, const Field &F, char *buf)
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

template <class Element, class Row>
template <class Field>
istream &SparseMatrix0ReadWriteHelper<Element, Row>
	::readGuillaume (SparseMatrix0Base<Element, Row> &A, istream &is, const Field &F, char *buf)
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

template <class Element, class Row>
template <class Field>
istream &SparseMatrix0ReadWriteHelper<Element, Row>
	::readPretty (SparseMatrix0Base<Element, Row> &A, istream &is, const Field &F, char *buf)
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

// Versions using standard extractor operators

template <class Element, class Row>
istream &SparseMatrix0ReadWriteHelper<Element, Row>
	::readTurner (SparseMatrix0Base<Element, Row> &A, istream &is, char *buf)
{
	size_t i, j;

	A._A.clear ();
	A._A.resize (A._m);

	do {
		std::istringstream str (buf);

		str >> i;

		if (i == (size_t) -1) break; // return also if row index is -1
		str >> j;
		str >> A.getEntry (i, j);

		is.getline (buf, 80);
	} while (is);

	return is;
}

template <class Element, class Row>
istream &SparseMatrix0ReadWriteHelper<Element, Row>
	::readGuillaume (SparseMatrix0Base<Element, Row> &A, istream &is, char *buf)
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
		is >> A.getEntry (i - 1, j - 1);
	}

	return is;

}

template <class Element, class Row>
istream &SparseMatrix0ReadWriteHelper<Element, Row>
	::readPretty (SparseMatrix0Base<Element, Row> &A, istream &is, char *buf)
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
			str >> a_ij;

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

// Versions using field read

template <class Element, class Row>
template <class Field>
istream &SparseMatrix0ReadWriteHelper<Element, Row>
	::read (SparseMatrix0Base<Element, Row> &A, istream &is, const Field &F, Format format)
{
	char buf[80];
	size_t i, j;
	Element a_ij;
	char c;

	is.getline (buf, 80);
	std::istringstream str (buf);

	switch (format) {
	    case FORMAT_DETECT:
		do str >> c; while (isspace (c));

		if (c == '[')
			readPretty (A, is, buf);
		else if (isdigit (c)) {
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
	    case FORMAT_PRETTY:
		return readPretty (A, is, F, buf);
	}
}

template <class Element, class Row>
template <class Field>
ostream &SparseMatrix0ReadWriteHelper<Element, Row>
	::write (const SparseMatrix0Base<Element, Row> &A, ostream &os, const Field &F, Format format)
{
	typename SparseMatrix0Base<Element, Row>::Rep::const_iterator i;
	typename Row::const_iterator j;
	typename Field::Element zero;
	size_t i_idx, j_idx;
	int col_width;
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
				os << endl;
			}
		}

		break;

	    case FORMAT_GUILLAUME:
		os << A._m << ' ' << A._n << " M" << endl;

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			for (j = i->begin (), j_idx = 0; j != i->end (); j++, j_idx++) {
				os << i_idx + 1 << ' ' << j->first + 1 << ' ';
				F.write (os, j->second);
				os << endl;
			}
		}

		os << "0 0 0" << endl;

		break;

	    case FORMAT_PRETTY:
		F.characteristic (c);
		col_width = logp (c, 10) + 1;
		F.init (zero, 0);

		for (i = A._A.begin (), i_idx = 0; i != A._A.end (); i++, i_idx++) {
			commentator.indent (os);

			os << "  [";

			j = i->begin ();

			for (j_idx = 0; j_idx < A._n; j_idx++) {
				os.width (col_width);

				if (j == i->end () || j_idx != j->first)
					F.write (os, zero);
				else {
					F.write (os, j->second);
					j++;
				}

				if (j_idx < A._n - 1)
					os << ' ';
			}

			os << ']' << endl;
		}

		break;
	}

	return os;
}

// Versions using standard operators

template <class Element, class Row>
istream &SparseMatrix0ReadWriteHelper<Element, Row>
	::read (SparseMatrix0Base<Element, Row> &A, istream &is, Format format)
{
	char buf[80];
	size_t i, j;
	Element a_ij;
	char c;

	is.getline (buf, 80);

	switch (format) {
	    case FORMAT_DETECT:
		std::istringstream str (buf);

		do str >> c; while (isspace (c));

		if (c == '[')
			readPretty (A, is, buf);
		else if (isdigit (c)) {
			if (buf[strlen (buf) - 2] == 'M')
				return readGuillaume (A, is, buf);
			else
				return readTurner (A, is, buf);
		} else
			throw InvalidMatrixInput ();

	    case FORMAT_TURNER:
		return readTurner (A, is, buf);
	    case FORMAT_GUILLAUME:
		return readGuillaume (A, is, buf);
	    case FORMAT_PRETTY:
		return readPretty (A, is, buf);
	}
}

template <class Element, class Row>
ostream &SparseMatrix0ReadWriteHelper<Element, Row>
	::write (const SparseMatrix0Base<Element, Row> &A, ostream &os, Format format)
{
	typename SparseMatrix0Base<Element, Row>::Rep::const_iterator i;
	typename Row::const_iterator j;
	size_t i_idx, j_idx;
	int col_width;
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
		for (i = _A.begin (), i_idx = 0; i != _A.end (); i++, i_idx++)
			for (j = i->begin (), j_idx = 0; j != i->end (); j++, j_idx++)
				os << i << ' ' << j->first << ' ' << j->second << endl;

		break;

	    case FORMAT_GUILLAUME:
		os << _m << ' ' << _n << " M" << endl;

		for (i = _A.begin (), i_idx = 0; i != _A.end (); i++, i_idx++)
			for (j = i->begin (), j_idx = 0; j != i->end (); j++, j_idx++)
				os << i + 1 << ' ' << j->first + 1 << ' ' << j->second << endl;

		os << "0 0 0" << endl;

		break;

	    case FORMAT_PRETTY:
		col_width = 3;      // FIXME

		for (i = _A.begin (), i_idx = 0; i != _A.end (); i++, i_idx++) {
			commentator.indent (os);

			os << "  [";

			j = i->begin ();

			for (j_idx = 0; j_idx < _n; j_idx++) {
				os.width (col_width);

				if (j == i->end () || j_idx != j->first)
					os << 0;
				else {
					os << j->second;
					j++;
				}

				if (j < _n - 1)
					os << ' ';
			}

			os << ']' << endl;
		}

		break;
	}

	return os;
}

template <class Element, class Row, class VectorTrait>
void SparseMatrix0Base<Element, Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
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

template <class Element, class Row, class VectorTrait>
Element &SparseMatrix0Base<Element, Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
	::refEntry (size_t i, size_t j) 
{
	static Element zero;

	Row &v = _A[i];
	typename Row::iterator iter;

	if (v.size () == 0) {
		v.push_back (pair <size_t, Element> (j, zero));
		return v.front ().second;
	} else {
		iter = std::lower_bound (v.begin (), v.end (), j, VectorWrapper::CompareSparseEntries<Element> ());

		if (iter == v.end () || iter->first != j)
			iter = v.insert (iter, pair <size_t, Element> (j, zero));

		return iter->second;
	}
}

template <class Element, class Row, class VectorTrait>
const Element &SparseMatrix0Base<Element, Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
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

template <class Element, class Row, class VectorTrait>
const Element &SparseMatrix0Base<Element, Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >
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

} // namespace LinBox

#endif // __SPARSE_BASE_INL
