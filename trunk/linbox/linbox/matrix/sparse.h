/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/sparse-base.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               1999-2001 William J Turner,
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-08-06  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed to sparse-base.h from sparse0-base.h
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Refactoring:
 *   - Eliminated SparseMatrix0Aux and moved that functionality into Sparse0
 *   - Made SparseMatrix0Base parameterized only on the element type
 *   - New read/write implementations for SparseMatrix0Base, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, RawIterator, and
 *     RawIndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __SPARSE_BASE_H
#define __SPARSE_BASE_H

#include <vector>
#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox/field/archetype.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"

namespace LinBox
{

/** Exception class for invalid matrix input
 */

class InvalidMatrixInput {};

// Forward declaration
template <class Element, class Row, class Trait = typename VectorTraits<Row>::VectorCategory>
class SparseMatrix0Base;

// Small helper class to make read and write easier
template <class Element, class Row>
class SparseMatrix0ReadWriteHelper
{
	template <class Field>
	static std::istream &readTurner    (SparseMatrix0Base<Element, Row> &A, std::istream &is, const Field &F, char *buf);
	template <class Field>
	static std::istream &readGuillaume (SparseMatrix0Base<Element, Row> &A, std::istream &is, const Field &F, char *buf);
	template <class Field>
	static std::istream &readPretty    (SparseMatrix0Base<Element, Row> &A, std::istream &is, const Field &F, char *buf);

	static std::istream &readTurner    (SparseMatrix0Base<Element, Row> &A, std::istream &is, char *buf);
	static std::istream &readGuillaume (SparseMatrix0Base<Element, Row> &A, std::istream &is, char *buf);
	static std::istream &readPretty    (SparseMatrix0Base<Element, Row> &A, std::istream &is, char *buf);

    public:
	enum Format {
		FORMAT_DETECT, FORMAT_GUILLAUME, FORMAT_TURNER, FORMAT_PRETTY
	};

	template <class Field>
	static std::istream &read (SparseMatrix0Base<Element, Row> &A, std::istream &is, const Field &F, Format format);
	template <class Field>
	static std::ostream &write (const SparseMatrix0Base<Element, Row> &A, std::ostream &os, const Field &F, Format format);

	static std::istream &read (SparseMatrix0Base<Element, Row> &A, std::istream &is, Format format);
	static std::ostream &write (const SparseMatrix0Base<Element, Row> &A, std::ostream &os, Format format);
};

/** Sparse matrix container
 * This class acts as a generic row-wise container for sparse
 * matrices. It is designed to provide various methods to access the
 * entries of the matrix. In particular, it does not meet the black box
 * archetype; see \ref{SparseMatrix0} for an appropriate sparse matrix
 * black box.
 *
 * @param Element Element type
 * @param Row     LinBox sparse vector type to use for rows of matrix
 */
template <class Element, class Row, class Trait>
class SparseMatrix0Base
{
    public:

	typedef typename std::vector<Row> Rep;

	/** Constructor.
	 * Note: the copy constructor and operator= will work as intended
	 *       because of STL's container design
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
	SparseMatrix0Base (size_t m, size_t n)
		: _m (m), _n (n), _A (m)
	{}

	/** Copy constructor.
	 */
	SparseMatrix0Base (const SparseMatrix0Base<Element, Row, Trait> &A)
		: _m (A._m), _n (A._n), _A (A._A) 
	{}

	/** Destructor. */
	~SparseMatrix0Base () {}

	/** Retreive row dimensions of Sparsemat matrix.
	 * @return integer number of rows of SparseMatrix0Base matrix.
	 */
	size_t rowdim () const { return _m; }

	/** Retreive column dimensions of Sparsemat matrix.
	 * @return integer number of columns of SparseMatrix0Base matrix.
	 */
	size_t coldim () const { return _n; }

	/** @name Input and output
	 */
	//@{

	/** Matrix file formats
	 */
	enum Format {
		FORMAT_DETECT, FORMAT_GUILLAUME, FORMAT_TURNER, FORMAT_PRETTY
	};

	/** Read a matrix from the given input stream using field read/write
	 * @param is Input stream from which to read the matrix
	 * @param F Field with which to read
	 * @param format Format of input matrix
	 */
	template <class Field>
	std::istream &read (std::istream &is, const Field &F, Format format = FORMAT_DETECT);

	/** Read a matrix from the given input stream using standard operators
	 * @param is Input stream from which to read the matrix
	 * @param format Format of input matrix
	 */
	std::istream &read (std::istream &is, Format format = FORMAT_DETECT);

	/** Write a matrix to the given output stream using field read/write
	 * @param os Output stream to which to write the matrix
	 * @param F Field with which to write
	 * @param format Format with which to write
	 */
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, Format format = FORMAT_GUILLAUME) const;

	/** Write a matrix to the given output stream using standard operators
	 * @param os Output stream to which to write the matrix
	 * @param format Format with which to write
	 */
	std::ostream &write (std::ostream &os, Format format = FORMAT_GUILLAUME) const;

	//@}

	/** @name Access to matrix elements
	 */
	//@{

	/** Set an individual entry
	 * Setting the entry to 0 will remove it from the matrix
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @value Value of the new entry
	 */
	void setEntry (size_t i, size_t j, const Element &value);

	/** Get a writeable reference to an entry in the matrix
	 * If there is no entry at the position (i, j), then a new entry
	 * with a value of zero is inserted and a reference  to it is
	 * returned.
	 * @param i Row index of entry
	 * @param j Column index of entry
	 * @return Reference to matrix entry
	 */
	Element &refEntry (size_t i, size_t j);

	/** Get a read-only individual entry from the matrix
	 * @param i Row index
	 * @param j Column index
	 * @return Const reference to matrix entry
	 */
	const Element &getEntry (size_t i, size_t j) const;

	/** Get an entry and store it in the given value
	 * This form is more in the Linbox style and is provided for interface
	 * compatibility with other parts of the library
	 * @param x Element in which to store result
	 * @param i Row index
	 * @param j Column index
	 * @return Reference to x
	 */
	Element &getEntry (Element &x, size_t i, size_t j) const
		{ x = getEntry (i, j); return x; }

	/** @name Columns of rows iterator
	 * The columns of row iterator gives each of the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in sparse sequence format
	 */

	typedef typename Rep::iterator ColOfRowsIterator;
	typedef typename Rep::const_iterator ConstColOfRowsIterator;

	ColOfRowsIterator rowsBegin ()
		{ return _A.begin (); }
	ColOfRowsIterator rowsEnd ()
		{ return _A.end (); }
	ConstColOfRowsIterator rowsBegin () const 
		{ return _A.begin (); }
	ConstColOfRowsIterator rowsEnd () const
		{ return _A.end (); }

	/** @name Raw iterator
	 * The raw iterator is a method for accessing all nonzero
	 * entries in the matrix in some unspecified order. This can be
	 * used, e.g. to reduce all matrix entries modulo a prime before
	 * passing the matrix into an algorithm.
	 */

	class RawIterator
	{
	    public:
		typedef Element value_type;

		RawIterator (Rep &A, const typename Rep::iterator &i, const typename Row::iterator &j)
			: _A (A), _i (i), _j (j)
		{}

		RawIterator (const RawIterator &iter)
			: _A (iter._A), _i (iter._i), _j (iter._j)
		{}

		RawIterator &operator = (const RawIterator &iter) 
		{
			linbox_check (&_A == &iter._A);

			_i = iter._i;
			_j = iter._j;

			return *this;
		}

		bool operator == (const RawIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const RawIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		RawIterator &operator ++ ()
		{
			if (++_j == _i->end ())
				if (++_i != _A.end ())
					_j = _i->begin ();

			return *this;
		}

		RawIterator &operator ++ (int)
		{
			RawIterator tmp = *this;
			++(*this);
			return tmp;
		}

		RawIterator &operator -- ()
		{
			if (_j == _i->begin ())
				_j = (--_i)->end ();
			--_j;
			return *this;
		}

		RawIterator &operator -- (int)
		{
			RawIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return _j->second; }
		const value_type &operator * () const
			{ return _j->second; }
		value_type *operator -> ()
			{ return &(_j->second); }
		const value_type *operator -> () const
			{ return &(_j->second); }

	    private:
		typename Rep::iterator  _i;
		typename Row::iterator  _j;
		Rep                    &_A;
	};

	typedef const RawIterator ConstRawIterator;

	RawIterator rawBegin ()
		{ return RawIterator (_A, _A.begin (), _A.front ().begin ()); }
	RawIterator rawEnd ()
		{ return RawIterator (_A, _A.end (), _A.back ().end ()); }
	ConstRawIterator rawBegin () const
		{ return RawIterator (_A, _A.begin (), _A.front ().begin ()); }
	ConstRawIterator rawEnd () const
		{ return RawIterator (_A, _A.end (), _A.back ().end ()); }

	/** @name Index iterator
	 * The index iterator gives the row, column indices of all matrix
	 * elements in the same order as the raw iterator above.
	 */

	class RawIndexIterator
	{
	    public:
		typedef std::pair<size_t, size_t> value_type;

		RawIndexIterator (Rep &A, size_t idx, const typename Rep::iterator &i, const typename Row::iterator &j)
			: _A (A), _i (i), _j (j), _curr (idx, j->second)
		{}

		RawIndexIterator (const RawIndexIterator &iter)
			: _A (iter._A), _i (iter._i), _j (iter._j), _curr (iter._curr)
		{}

		RawIndexIterator &operator = (const RawIndexIterator &iter) 
		{
			linbox_check (&_A == &iter._A);

			_i = iter._i;
			_j = iter._j;
			_curr = iter._curr;

			return *this;
		}

		bool operator == (const RawIndexIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const RawIndexIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		RawIndexIterator &operator ++ ()
		{
			if (++_j == _i->end ()) {
				if (++_i != _A.end ()) {
					_j = _i->begin ();
					_curr.first++;
				}
			}
			else
				_curr.second++;

			return *this;
		}

		RawIndexIterator &operator ++ (int)
		{
			RawIndexIterator tmp = *this;
			++(*this);
			return tmp;
		}

		RawIndexIterator &operator -- ()
		{
			if (_j == _i->begin ()) {
				_j = (--_i)->end ();
				_curr.first--;
			}
			else
				_curr.second--;
			--_j;
			return *this;
		}

		RawIndexIterator &operator -- (int)
		{
			RawIndexIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return _curr; }
		const value_type &operator * () const
			{ return _curr; }
		value_type *operator -> ()
			{ return &(_curr); }
		const value_type *operator -> () const
			{ return &(_curr); }

	    private:
		typename Rep::iterator  _i;
		typename Row::iterator  _j;
		Rep                    &_A;
		value_type              _curr;
	};

	typedef const RawIndexIterator ConstRawIndexIterator;

	RawIndexIterator indexBegin ()
		{ return RawIndexIterator (_A, 0, _A.begin (), _A.front ().begin ()); }
	RawIndexIterator indexEnd ()
		{ return RawIndexIterator (_A, _m, _A.end (), _A.back ().end ()); }
	ConstRawIndexIterator indexBegin () const
		{ return RawIndexIterator (_A, 0, _A.begin (), _A.front ().begin ()); }
	ConstRawIndexIterator indexEnd () const
		{ return RawIndexIterator (_A, _m, _A.end (), _A.back ().end ()); }

	/** Retrieve a row as a writeable reference
	 * @param i Row index
	 */
	Row &getRow (size_t i)
		{ return _A[i]; }

	//@}

    protected:

	friend class SparseMatrix0ReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;
};

/* Specialization for sparse sequence vectors */

template <class Element, class Row, class VectorTrait>
class SparseMatrix0Base<Element, Row, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
{
    public:

	typedef std::vector<Row> Rep;

	SparseMatrix0Base (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrix0Base (const SparseMatrix0Base<Element, Row, VectorTrait> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}
	~SparseMatrix0Base () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	enum Format {
		FORMAT_DETECT, FORMAT_GUILLAUME, FORMAT_TURNER, FORMAT_PRETTY
	};

	template <class Field>
	std::istream &read (std::istream &is, const Field &F, Format format = FORMAT_DETECT)
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::read
			  (*this, is, F, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::read
			  (*this, is, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, Format format = FORMAT_GUILLAUME) const
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::write
			  (*this, os, F, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_GUILLAUME) const
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::write
			  (*this, is, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }

	void           setEntry (size_t i, size_t j, const Element &value);
	Element       &refEntry (size_t i, size_t j);
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const
			{ x = getEntry (i, j); return x; }

	typedef typename Rep::iterator ColOfRowsIterator;
	typedef typename Rep::const_iterator ConstColOfRowsIterator;

	ConstColOfRowsIterator rowsBegin () const 
		{ return _A.begin (); }
	ConstColOfRowsIterator rowsEnd () const
		{ return _A.end (); }
	ColOfRowsIterator rowsBegin ()
		{ return _A.begin (); }
	ColOfRowsIterator rowsEnd ()
		{ return _A.end (); }

	class RawIterator
	{
	    public:
		RawIterator (Rep &A, const typename Rep::iterator &i, const typename Row::iterator &j)
			: _A (A), _i (i), _j (j)
		{}

		RawIterator (const RawIterator &iter)
			: _A (iter._A), _i (iter._i), _j (iter._j)
		{}

		RawIterator &operator = (const RawIterator &iter) 
		{
			linbox_check (&_A == &iter._A);

			_i = iter._i;
			_j = iter._j;

			return *this;
		}

		bool operator == (const RawIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const RawIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		RawIterator &operator ++ ()
		{
			if (++_j == _i->end ())
				if (++_i != _A.end ())
					_j = _i->begin ();
			return *this;
		}

		RawIterator &operator ++ (int)
		{
			RawIterator tmp = *this;
			++(*this);
			return tmp;
		}

		RawIterator &operator -- ()
		{
			if (_j == _i->begin ())
				_j = (--_i)->end ();
			--_j;
			return *this;
		}

		RawIterator &operator -- (int)
		{
			RawIterator tmp = *this;
			--(*this);
			return tmp;
		}

		Element &operator * ()
			{ return _j->second; }
		const Element &operator * () const
			{ return _j->second; }
		Element *operator -> ()
			{ return &(_j->second); }
		const Element *operator -> () const
			{ return &(_j->second); }

	    private:
		typename Rep::iterator  _i;
		typename Row::iterator  _j;
		Rep                    &_A;
	};

	typedef const RawIterator ConstRawIterator;

	RawIterator rawBegin ()
		{ return RawIterator (_A, _A.begin (), _A.front ().begin ()); }
	RawIterator rawEnd ()
		{ return RawIterator (_A, _A.end (), _A.back ().end ()); }
	ConstRawIterator rawBegin () const
		{ return RawIterator (_A, _A.begin (), _A.front ().begin ()); }
	ConstRawIterator rawEnd () const
		{ return RawIterator (_A, _A.end (), _A.back ().end ()); }

	class RawIndexIterator
	{
	    public:
		typedef std::pair<size_t, size_t> value_type;

		RawIndexIterator (Rep &A, size_t idx, const typename Rep::iterator &i, const typename Row::iterator &j)
			: _A (A), _i (i), _j (j), _curr (idx, j->second)
		{}

		RawIndexIterator (const RawIndexIterator &iter)
			: _A (iter._A), _i (iter._i), _j (iter._j), _curr (iter._curr)
		{}

		RawIndexIterator &operator = (const RawIndexIterator &iter) 
		{
			linbox_check (&_A == &iter._A);

			_i = iter._i;
			_j = iter._j;
			_curr = iter._curr;

			return *this;
		}

		bool operator == (const RawIndexIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const RawIndexIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		RawIndexIterator &operator ++ ()
		{
			if (++_j == _i->end ()) {
				if (++_i != _A.end ()) {
					_j = _i->begin ();
					_curr.first++;
				}
			}
			else
				_curr.second++;

			return *this;
		}

		RawIndexIterator &operator ++ (int)
		{
			RawIndexIterator tmp = *this;
			++(*this);
			return tmp;
		}

		RawIndexIterator &operator -- ()
		{
			if (_j == _i->begin ()) {
				_j = (--_i)->end ();
				_curr.first--;
			}

			--_j;
			_curr.second = _j->second;
			return *this;
		}

		RawIndexIterator &operator -- (int)
		{
			RawIndexIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return _curr; }
		const value_type &operator * () const
			{ return _curr; }
		value_type *operator -> ()
			{ return &(_curr); }
		const value_type *operator -> () const
			{ return &(_curr); }

	    private:
		typename Rep::iterator  _i;
		typename Row::iterator  _j;
		Rep                    &_A;
		value_type              _curr;
	};

	typedef const RawIndexIterator ConstRawIndexIterator;

	RawIndexIterator indexBegin ()
		{ return RawIndexIterator (_A, 0, _A.begin (), _A.front ().begin ()); }
	RawIndexIterator indexEnd ()
		{ return RawIndexIterator (_A, _m, _A.end (), _A.back ().end ()); }
	ConstRawIndexIterator indexBegin () const
		{ return RawIndexIterator (_A, 0, _A.begin (), _A.front ().begin ()); }
	ConstRawIndexIterator indexEnd () const
		{ return RawIndexIterator (_A, _m, _A.end (), _A.back ().end ()); }

	Row &getRow (size_t i)
		{ return _A[i]; }

    protected:

	friend class SparseMatrix0ReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;
};

/* Specialization for sparse associative vectors */

template <class Element, class Row, class VectorTrait>
class SparseMatrix0Base<Element, Row, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >
{
    public:

	typedef std::vector<Row> Rep;

	SparseMatrix0Base (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrix0Base (const SparseMatrix0Base<Element, Row, VectorTrait> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}
	~SparseMatrix0Base () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	enum Format {
		FORMAT_DETECT, FORMAT_GUILLAUME, FORMAT_TURNER, FORMAT_PRETTY
	};

	template <class Field>
	std::istream &read (std::istream &is, const Field &F, Format format = FORMAT_DETECT)
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::read
			  (*this, is, F, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::read
			  (*this, is, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, Format format = FORMAT_GUILLAUME) const
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::write
			  (*this, os, F, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_GUILLAUME) const
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::write
			  (*this, is, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }

	void           setEntry (size_t i, size_t j, const Element &value) { _A[i][j] = value; }
	Element       &refEntry (size_t i, size_t j)                       { return _A[i][j]; }
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const     { x = _A[i][j]; return x; }

	typedef typename Rep::iterator ColOfRowsIterator;
	typedef typename Rep::const_iterator ConstColOfRowsIterator;

	ConstColOfRowsIterator rowsBegin () const 
		{ return _A.begin (); }
	ConstColOfRowsIterator rowsEnd () const
		{ return _A.end (); }
	ColOfRowsIterator rowsBegin ()
		{ return _A.begin (); }
	ColOfRowsIterator rowsEnd ()
		{ return _A.end (); }

	class RawIterator
	{
	    public:
		RawIterator (Rep &A, const typename Rep::iterator &i, const typename Row::iterator &j)
			: _A (A), _i (i), _j (j)
		{}

		RawIterator (const RawIterator &iter)
			: _A (iter._A), _i (iter._i), _j (iter._j)
		{}

		RawIterator &operator = (const RawIterator &iter) 
		{
			linbox_check (&_A == &iter._A);

			_i = iter._i;
			_j = iter._j;

			return *this;
		}

		bool operator == (const RawIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const RawIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		RawIterator &operator ++ ()
		{
			if (++_j == _i->end ())
				if (++_i != _A.end ())
					_j = _i->begin ();
			return *this;
		}

		RawIterator &operator ++ (int)
		{
			RawIterator tmp = *this;
			++(*this);
			return tmp;
		}

		RawIterator &operator -- ()
		{
			if (_j == _i->begin ())
				_j = (--_i)->end ();
			--_j;
			return *this;
		}

		RawIterator &operator -- (int)
		{
			RawIterator tmp = *this;
			--(*this);
			return tmp;
		}

		Element &operator * ()
			{ return _j->second; }
		const Element &operator * () const
			{ return _j->second; }
		Element *operator -> ()
			{ return &(_j->second); }
		const Element *operator -> () const
			{ return &(_j->second); }

	    private:
		typename Rep::iterator  _i;
		typename Row::iterator  _j;
		Rep                    &_A;
	};

	typedef const RawIterator ConstRawIterator;

	RawIterator rawBegin ()
		{ return RawIterator (_A, _A.begin (), _A.front ().begin ()); }
	RawIterator rawEnd ()
		{ return RawIterator (_A, _A.end (), _A.back ().end ()); }
	ConstRawIterator rawBegin () const
		{ return RawIterator (_A, _A.begin (), _A.front ().begin ()); }
	ConstRawIterator rawEnd () const
		{ return RawIterator (_A, _A.end (), _A.back ().end ()); }

	class RawIndexIterator
	{
	    public:
		typedef std::pair<size_t, size_t> value_type;

		RawIndexIterator (Rep &A, size_t idx, const typename Rep::iterator &i, const typename Row::iterator &j)
			: _A (A), _i (i), _j (j), _curr (idx, j->second)
		{}

		RawIndexIterator (const RawIndexIterator &iter)
			: _A (iter._A), _i (iter._i), _j (iter._j), _curr (iter._curr)
		{}

		RawIndexIterator &operator = (const RawIndexIterator &iter) 
		{
			linbox_check (&_A == &iter._A);

			_i = iter._i;
			_j = iter._j;
			_curr = iter._curr;

			return *this;
		}

		bool operator == (const RawIndexIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const RawIndexIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		RawIndexIterator &operator ++ ()
		{
			if (++_j == _i->end ()) {
				if (++_i != _A.end ()) {
					_j = _i->begin ();
					_curr.first++;
				}
			}
			else
				_curr.second++;

			return *this;
		}

		RawIndexIterator &operator ++ (int)
		{
			RawIndexIterator tmp = *this;
			++(*this);
			return tmp;
		}

		RawIndexIterator &operator -- ()
		{
			if (_j == _i->begin ()) {
				_j = (--_i)->end ();
				_curr.first--;
			}

			--_j;
			_curr.second = _j->second;
			return *this;
		}

		RawIndexIterator &operator -- (int)
		{
			RawIndexIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return _curr; }
		const value_type &operator * () const
			{ return _curr; }
		value_type *operator -> ()
			{ return &(_curr); }
		const value_type *operator -> () const
			{ return &(_curr); }

	    private:
		typename Rep::iterator  _i;
		typename Row::iterator  _j;
		Rep                    &_A;
		value_type              _curr;
	};

	typedef const RawIndexIterator ConstRawIndexIterator;

	RawIndexIterator indexBegin ()
		{ return RawIndexIterator (_A, 0, _A.begin (), _A.front ().begin ()); }
	RawIndexIterator indexEnd ()
		{ return RawIndexIterator (_A, _m, _A.end (), _A.back ().end ()); }
	ConstRawIndexIterator indexBegin () const
		{ return RawIndexIterator (_A, 0, _A.begin (), _A.front ().begin ()); }
	ConstRawIndexIterator indexEnd () const
		{ return RawIndexIterator (_A, _m, _A.end (), _A.back ().end ()); }

	Row &getRow (size_t i)
		{ return _A[i]; }

    protected:

	friend class SparseMatrix0ReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;
};

/* Specialization for sparse parallel vectors */

template <class Element, class Row, class VectorTrait>
class SparseMatrix0Base<Element, Row, VectorCategories::SparseParallelVectorTag<VectorTrait> >
{
    public:

	typedef std::vector<Row> Rep;

	SparseMatrix0Base (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrix0Base (const SparseMatrix0Base<Element, Row, VectorTrait> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}
	~SparseMatrix0Base () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }

	enum Format {
		FORMAT_DETECT, FORMAT_GUILLAUME, FORMAT_TURNER, FORMAT_PRETTY
	};

	template <class Field>
	std::istream &read (std::istream &is, const Field &F, Format format = FORMAT_DETECT)
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::read
			  (*this, is, F, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	std::istream &read (std::istream &is, Format format = FORMAT_DETECT)
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::read
			  (*this, is, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, Format format = FORMAT_GUILLAUME) const
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::write
			  (*this, os, F, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }
	std::ostream &write (std::ostream &os, Format format = FORMAT_GUILLAUME) const
		{ return SparseMatrix0ReadWriteHelper<Element, Row>::write
			  (*this, is, (typename SparseMatrix0ReadWriteHelper<Element, Row>::Format) format); }

	void           setEntry (size_t i, size_t j, const Element &value);
	Element       &refEntry (size_t i, size_t j);
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const
			{ x = getEntry (i, j); return x; }

	typedef typename Rep::iterator ColOfRowsIterator;
	typedef typename Rep::const_iterator ConstColOfRowsIterator;

	ConstColOfRowsIterator rowsBegin () const 
		{ return _A.begin (); }
	ConstColOfRowsIterator rowsEnd () const
		{ return _A.end (); }
	ColOfRowsIterator rowsBegin ()
		{ return _A.begin (); }
	ColOfRowsIterator rowsEnd ()
		{ return _A.end (); }

	class RawIterator
	{
	    public:
		RawIterator (Rep &A, const typename Rep::iterator &i, const typename Row::second_type::iterator &j)
			: _A (A), _i (i), _j (j)
		{}

		RawIterator (const RawIterator &iter)
			: _A (iter._A), _i (iter._i), _j (iter._j)
		{}

		RawIterator &operator = (const RawIterator &iter) 
		{
			linbox_check (&_A == &iter._A);

			_i = iter._i;
			_j = iter._j;

			return *this;
		}

		bool operator == (const RawIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const RawIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		RawIterator &operator ++ ()
		{
			if (++_j == _i->second.end ())
				if (++_i != _A.end ())
					_j = _i->second.begin ();
			return *this;
		}

		RawIterator &operator ++ (int)
		{
			RawIterator tmp = *this;
			++(*this);
			return tmp;
		}

		RawIterator &operator -- ()
		{
			if (_j == _i->second.begin ())
				_j = (--_i)->second.end ();
			--_j;
			return *this;
		}

		RawIterator &operator -- (int)
		{
			RawIterator tmp = *this;
			--(*this);
			return tmp;
		}

		Element &operator * ()
			{ return *_j; }
		const Element &operator * () const
			{ return *_j; }
		Element *operator -> ()
			{ return &(*_j); }
		const Element *operator -> () const
			{ return &(*_j); }

	    private:
		typename Rep::iterator               _i;
		typename Row::second_type::iterator  _j;
		Rep                                 &_A;
	};

	typedef const RawIterator ConstRawIterator;

	RawIterator rawBegin ()
		{ return RawIterator (_A, _A.begin (), _A.front ().second.begin ()); }
	RawIterator rawEnd ()
		{ return RawIterator (_A, _A.end (), _A.back ().second.end ()); }
	ConstRawIterator rawBegin () const
		{ return RawIterator (_A, _A.begin (), _A.front ().second.begin ()); }
	ConstRawIterator rawEnd () const
		{ return RawIterator (_A, _A.end (), _A.back ().second.end ()); }

	class RawIndexIterator
	{
	    public:
		typedef std::pair<size_t, size_t> value_type;

		RawIndexIterator (Rep &A, size_t idx, const typename Rep::iterator &i, const typename Row::first_type::iterator &j)
			: _A (A), _i (i), _j (j), _curr (idx, *j)
		{}

		RawIndexIterator (const RawIndexIterator &iter)
			: _A (iter._A), _i (iter._i), _j (iter._j), _curr (iter._curr)
		{}

		RawIndexIterator &operator = (const RawIndexIterator &iter) 
		{
			linbox_check (&_A == &iter._A);

			_i = iter._i;
			_j = iter._j;
			_curr = iter._curr;

			return *this;
		}

		bool operator == (const RawIndexIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const RawIndexIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		RawIndexIterator &operator ++ ()
		{
			if (++_j == _i->first.end ()) {
				if (++_i != _A.end ()) {
					_j = _i->first.begin ();
					_curr.first++;
				}
			}
			else
				_curr.second++;

			return *this;
		}

		RawIndexIterator &operator ++ (int)
		{
			RawIndexIterator tmp = *this;
			++(*this);
			return tmp;
		}

		RawIndexIterator &operator -- ()
		{
			if (_j == _i->first.begin ()) {
				_j = (--_i)->first.end ();
				_curr.first--;
			}

			--_j;
			_curr.second = _j->second;
			return *this;
		}

		RawIndexIterator &operator -- (int)
		{
			RawIndexIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return _curr; }
		const value_type &operator * () const
			{ return _curr; }
		value_type *operator -> ()
			{ return &(_curr); }
		const value_type *operator -> () const
			{ return &(_curr); }

	    private:
		typename Rep::iterator            _i;
		typename Row::first_typeiterator  _j;
		Rep                              &_A;
		value_type                        _curr;
	};

	typedef const RawIndexIterator ConstRawIndexIterator;

	RawIndexIterator indexBegin ()
		{ return RawIndexIterator (_A, 0, _A.begin (), _A.front ().first.begin ()); }
	RawIndexIterator indexEnd ()
		{ return RawIndexIterator (_A, _m, _A.end (), _A.back ().first.end ()); }
	ConstRawIndexIterator indexBegin () const
		{ return RawIndexIterator (_A, 0, _A.begin (), _A.front ().first.begin ()); }
	ConstRawIndexIterator indexEnd () const
		{ return RawIndexIterator (_A, _m, _A.end (), _A.back ().first.end ()); }

	Row &getRow (size_t i)
		{ return _A[i]; }

    protected:

	friend class SparseMatrix0ReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;
};

template <class Element, class Row>
std::ostream &operator << (std::ostream &os, const SparseMatrix0Base<Element, Row> &A)
	{ return A.write (os); }

template <class Element, class Row>
std::istream &operator >> (std::istream &is, SparseMatrix0Base<Element, Row> &A)
	{ return A.read (is); }

} // namespace LinBox

#include "linbox/blackbox/sparse-base.inl"

#endif // __SPARSE_BASE_H
