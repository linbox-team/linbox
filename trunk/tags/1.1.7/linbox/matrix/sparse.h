/* linbox/matrix/sparse.h
 * Copyright (C) 2001-2002 Bradford Hovinen
 *               1999-2001 William J Turner,
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 * 
 * --------------------------------------------------------
 * 2003-01-11  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 * Move from blackbox/sparse-base.h to matrix/sparse.h
 * ------------------------------------
 * 2002-11-28  Bradford Hovinen  <bghovinen@math.uwaterloo.ca>
 *
 *   - Renamed ColOfRowsIterator to RowIterator
 *   - Named template argument _Row rather than Row; add a typedef to Row
 *   - Named template argument _Element rather than Row; add a typedef to Element
 *   - Renamed RawIndexIterator as RawIndexedIterator, and adjusted to match
 *     interface in DenseMatrixBase
 * ------------------------------------
 * 2002-08-06  Bradford Hovinen  <hovinen@cis.udel.edu>
 *
 * Renamed to sparse-base.h from sparse0-base.h
 * ------------------------------------
 * Modified by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Refactoring:
 *   - Eliminated SparseMatrixAux and moved that functionality into Sparse0
 *   - Made SparseMatrixBase parameterized only on the element type
 *   - New read/write implementations for SparseMatrixBase, supporting multiple
 *     formats
 *   - Eliminated Gaussian elimination code
 *   - Added iterators, including ColOfRowsIterator, RawIterator, and
 *     RawIndexIterator
 *   - Eliminated operator []; added getEntry; changed put_value to setEntry
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_matrix_sparse_H
#define __LINBOX_matrix_sparse_H

#ifndef _SP_BB_VECTOR_
#include <vector>
#define _SP_BB_VECTOR_ std::vector
#endif

#include <utility>
#include <iostream>
#include <algorithm>

#include "linbox/linbox-config.h"
#include "linbox/blackbox/factory.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/util/matrix-stream.h"

namespace LinBox
{
	

/** Exception class for invalid matrix input
 */

class InvalidMatrixInput {};

// use this in place of NoField member, to avoid code duplication.
// better yet.  Don't do it at all.
/// Dummy field for conceptually unclear io.
template<class _Element>
class FieldIO 
{ public:
	typedef _Element Element;

	std::istream &read (std::istream &stream, Element &elt) const
		{ return stream >> elt; }
	std::ostream &write (std::ostream &stream, const Element &elt) const
		{ return stream << elt; }
};

// made global to avoid duplicate code. 
/// tags for SparseMatrixBase::read() and write()
enum FileFormatTag {
	FORMAT_DETECT, FORMAT_GUILLAUME, FORMAT_TURNER, FORMAT_MATLAB, FORMAT_MAPLE, FORMAT_PRETTY, FORMAT_MAGMACPT, FORMAT_ONE_BASED
};

// Forward declaration
template <class _Element,
	  class _Row   = typename RawVector<_Element>::Sparse,
	  class Trait  = typename VectorTraits<_Row>::VectorCategory>
class SparseMatrixBase;


// Small helper classes to make read and write easier
template <class _Element, class Row, class Trait = typename VectorTraits<Row>::VectorCategory>
class SparseMatrixWriteHelper
{
    public:
	typedef _Element Element;

	// Dummy class to avoid code duplication
	class NoField 
	{
	    public:
		typedef _Element Element;

		std::istream &read (std::istream &stream, Element &elt) const
			{ return stream >> elt; }
		std::ostream &write (std::ostream &stream, const Element &elt) const
			{ return stream << elt; }
	};

	template <class Field>
	static std::ostream &write (const SparseMatrixBase<Element, Row> &A, std::ostream &os, const Field &F, FileFormatTag format);
};

template <class Element, class Row, class Trait = typename VectorTraits<Row>::VectorCategory>
class SparseMatrixReadWriteHelper : public SparseMatrixWriteHelper<Element, Row, Trait>
{
	template <class Field>
	static std::istream &readTurner    (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf);
	template <class Field>
	static std::istream &readGuillaume (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf);
	template <class Field>
	static std::istream &readMatlab    (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf);
	template <class Field>
	static std::istream &readPretty    (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf);
	template <class Field>
	static std::istream &readMagmaCpt  (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, char *buf);

    public:

	template <class Field>
	static std::istream &read (SparseMatrixBase<Element, Row> &A, std::istream &is, const Field &F, FileFormatTag format);
};

// Specialization of the above for sparse parallel vectors
template <class _Element, class Row>
class SparseMatrixWriteHelper<_Element, Row, VectorCategories::SparseParallelVectorTag >
{
    public:
	typedef _Element Element;

	// Dummy class to avoid code duplication
	class NoField 
	{
	    public:
		typedef _Element Element;

		std::istream &read (std::istream &stream, Element &elt) const
			{ return stream >> elt; }
		std::ostream &write (std::ostream &stream, const Element &elt) const
			{ return stream << elt; }
	};

	template <class Field>
	static std::ostream &write (const SparseMatrixBase<Element, Row> &A, std::ostream &os, const Field &F, FileFormatTag format);
};

/** Sparse matrix container
 * This class acts as a generic row-wise container for sparse
 * matrices. It is designed to provide various methods to access the
 * entries of the matrix. In particular, it does not meet the black box
 * archetype; see \ref{SparseMatrix} for an appropriate sparse matrix
 * black box.
 *
 * @param Element Element type
 * @param Row     LinBox sparse vector type to use for rows of matrix
\ingroup matrix
 */
template <class _Element, class _Row, class Trait>
class SparseMatrixBase
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef typename _SP_BB_VECTOR_<Row> Rep;

	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
        struct rebind
        { typedef SparseMatrixBase<typename _Tp1::Element, _R1, Trait> other; };

	/** Constructor.
	 * Note: the copy constructor and operator= will work as intended
	 *       because of STL's container design
	 * @param  m  row dimension
	 * @param  n  column dimension
	 */
        SparseMatrixBase (size_t m, size_t n): _A(m), _m(m), _n(n) {};


	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrixBase ( MatrixStream<Field>& ms );
    

	/** Copy constructor.
	 */
	SparseMatrixBase (const SparseMatrixBase<Element, Row, Trait> &A);

	/** Convert constructor.
	 */
    	template<class VectorType>
	SparseMatrixBase (const SparseMatrixBase<Element, VectorType, Trait> &A);
	/** Destructor. */
	~SparseMatrixBase () {}

	/** Retreive row dimension of the matrix.
	 * @return integer number of rows of SparseMatrixBase matrix.
	 */
	size_t rowdim () const { return _m; }

	/** Retreive column dimension of matrix.
	 * @return integer number of columns of SparseMatrixBase matrix.
	 */
	size_t coldim () const { return _n; }

	/** Retreive number of elements in the matrix.
	 * @return integer number of elements of SparseMatrixBase matrix.
	 */
	size_t size () const { 
            size_t s(0);
            for(typename Rep::const_iterator it = _A.begin(); it != _A.end(); ++it)
                s+= LinBox::RawVector<_Element>::size(*it);
            return s;
        }
	/** Read a matrix from the given input stream using field read/write
	 * @param is Input stream from which to read the matrix
	 * @param F Field with which to read
	 * @param format Format of input matrix
	 */
	template <class Field>
	std::istream &read (std::istream &is, const Field &F, FileFormatTag format = FORMAT_DETECT);
	/** Read a matrix from the given input stream using standard operators
	 * @param is Input stream from which to read the matrix
	 * @param format Format of input matrix
	 */
        std::istream &read (std::istream &is, FileFormatTag format = FORMAT_DETECT);
    

	/** Write a matrix to the given output stream using field read/write
	 * @param os Output stream to which to write the matrix
	 * @param F Field with which to write
	 * @param format Format with which to write
	 */
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, FileFormatTag format = FORMAT_PRETTY) const;

	/** Write a matrix to the given output stream using standard operators
	 * @param os Output stream to which to write the matrix
	 * @param format Format with which to write
	 */
	std::ostream &write (std::ostream &os, FileFormatTag format = FORMAT_PRETTY) const;

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
	Element &getEntry (Element &x, size_t i, size_t j) const;

	/** @name Columns of rows iterator
	 * The columns of row iterator gives each of the rows of the
	 * matrix in ascending order. Dereferencing the iterator yields
	 * a row vector in sparse sequence format
	 */

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	RowIterator rowBegin ();
	RowIterator rowEnd ();
	ConstRowIterator rowBegin () const;
	ConstRowIterator rowEnd () const;

	/** @name Raw iterator
	 * The raw iterator is a method for accessing all nonzero
	 * entries in the matrix in some unspecified order. This can be
	 * used, e.g. to reduce all matrix entries modulo a prime before
	 * passing the matrix into an algorithm.
	 */

	class RawIterator;
	class ConstRawIterator;

	RawIterator rawBegin ();
	RawIterator rawEnd ();
	ConstRawIterator rawBegin () const;
	ConstRawIterator rawEnd () const;

	/** @name Index iterator
	 * The index iterator gives the row, column indices of all matrix
	 * elements in the same order as the raw iterator above. Its value type
	 * is an STL pair with the row and column indices, starting at 0, in the
	 * first and second positions, respectively.
	 */

	class RawIndexedIterator;
	class ConstRawIndexedIterator;

	RawIndexedIterator rawIndexedBegin ();
	RawIndexedIterator rawIndexedEnd ();
	ConstRawIndexedIterator rawIndexedBegin () const;
	ConstRawIndexedIterator rawIndexedEnd () const;

	/** Retrieve a row as a writeable reference
	 * @param i Row index
	 */
	Row &getRow (size_t i);

	/** Retrieve a row as a writeable reference
	 * @param i Row index
	 */
	Row &operator [] (size_t i);

	/** Retrieve a row as a read-only reference
	 * @param i Row index
	 */
	ConstRow &operator [] (size_t i) const;

	/** Compute the column density, i.e. the number of entries per column
	 * @param v Vector in which to store column density 
	 */
	template <class Vector>
	Vector &columnDensity (Vector &v) const;

	/** Construct the transpose of this matrix and place it in the
	 * matrix given
	 */
	SparseMatrixBase &transpose (SparseMatrixBase &AT) const;

    protected:
	
	friend class SparseMatrixWriteHelper<Element, Row>;
	friend class SparseMatrixReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;

    	template<class F, class R, class T> friend class SparseMatrixBase;
};

/* Specialization for sparse sequence vectors */

template <class _Element, class _Row>
class SparseMatrixBase<_Element, _Row, VectorCategories::SparseSequenceVectorTag >
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
        struct rebind
        { typedef SparseMatrixBase<typename _Tp1::Element, _R1, VectorCategories::SparseSequenceVectorTag> other; };

	SparseMatrixBase (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrixBase ( MatrixStream<Field>& ms );

	SparseMatrixBase (const SparseMatrixBase<Element, Row> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}

    	template<class VectorType>
	SparseMatrixBase (const SparseMatrixBase<Element, VectorType> &A)
		: _A(A._m), _m (A._m), _n (A._n) {
            typename Rep::iterator meit = this->_A.begin();
            typename SparseMatrixBase<Element, VectorType>::Rep::const_iterator copit = A._A.begin();
            for( ; meit != this->_A.end(); ++meit, ++copit)
                LinBox::RawVector<Element>::convert(*meit, *copit);
        }

	~SparseMatrixBase () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const { 
            size_t s(0);
            for(typename Rep::const_iterator it = _A.begin(); it != _A.end(); ++it)
                s+= LinBox::RawVector<_Element>::size(*it);
            return s;
        }

	template <class Field>
	std::istream &read (std::istream &is, const Field &F, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row>::read
			  (*this, is, F, format); }
	std::istream &read (std::istream &is, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row>::read
			  (*this, is, SparseMatrixReadWriteHelper<Element, Row>::NoField (),
			   format); }
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row>::write
			  (*this, os, F, format); }
	std::ostream &write (std::ostream &os, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row>::write
			  (*this, os, SparseMatrixReadWriteHelper<Element, Row>::NoField (),
			   format); }

	std::ostream &write(std::ostream &) const;

	void           setEntry (size_t i, size_t j, const Element &value);
	Element       &refEntry (size_t i, size_t j);
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const
			{ return x = getEntry (i, j); }

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	ConstRowIterator rowBegin () const 
		{ return _A.begin (); }
	ConstRowIterator rowEnd () const
		{ return _A.end (); }
	RowIterator rowBegin ()
		{ return _A.begin (); }
	RowIterator rowEnd ()
		{ return _A.end (); }

	template <class RepIterator, class RowIterator, class _I_Element>
	class _RawIterator
	{
	    public:
		typedef _I_Element value_type;

		_RawIterator (const RepIterator &i, const RowIterator &j, const RepIterator &A_end)
			: _i (i), _j (j), _A_end (A_end)
		{
			if( _i == _A_end ) return;
		 	while ( _j == _i->end () ) {
 				if (++_i == _A_end) return;
 				_j = _i->begin ();
			}
		}

		_RawIterator (const _RawIterator &iter)
			: _i (iter._i), _j (iter._j), _A_end (iter._A_end)
		{}

		_RawIterator () {}

		_RawIterator &operator = (const _RawIterator &iter) 
		{
			_i = iter._i;
			_j = iter._j;
			_A_end = iter._A_end;

			return *this;
		}

		bool operator == (const _RawIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const _RawIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		_RawIterator &operator ++ ()
		{
			++_j;
			while( _j == _i->end ()) {
 				if (++_i == _A_end) return *this;
 				_j = _i->begin ();
			}
			
			// if (++_j == _i->end ())
// 				if (++_i != _A_end)
// 					_j = _i->begin ();
			return *this;
		}

		_RawIterator operator ++ (int)
		{
			_RawIterator tmp = *this;
			++(*this);
			return tmp;
		}

		_RawIterator &operator -- ()
		{
			while (_j == _i->begin ())
				_j = (--_i)->end ();
			--_j;
			return *this;
		}

		_RawIterator operator -- (int)
		{
			_RawIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return _j->second; }
			// { return *(new value_type(_j->second)); }
			// Dan Roche 2005-7-7 I believe this was a memory leak.
		value_type *operator -> ()
			{ return &(_j->second); }
		const value_type &operator*() const { return _j->second; }
		const value_type *operator-> () const { return &(_j->second); }

	    private:
		RepIterator _i;
		RowIterator _j;
		RepIterator _A_end;
	};

	typedef _RawIterator<typename Rep::iterator, typename Row::iterator, Element> RawIterator;
	typedef _RawIterator<typename Rep::const_iterator, typename Row::const_iterator, const Element> ConstRawIterator;

	RawIterator rawBegin ()
		{ return RawIterator (_A.begin (), _A.front ().begin (), _A.end ()); }
	RawIterator rawEnd ()
		{ return RawIterator (_A.end (), _A.back ().end (), _A.end ()); }
	ConstRawIterator rawBegin () const
		{ return ConstRawIterator (_A.begin (), _A.front ().begin (), _A.end ()); }
	ConstRawIterator rawEnd () const
		{ return ConstRawIterator (_A.end (), _A.back ().end (), _A.end ()); }



        /* Generic trait for iterators without type */
    template<typename U>
    struct IteratorValueType {
        typedef typename U::value_type value_type;
    };
    
    template<typename X>
    struct IteratorValueType<const X*> {
        typedef X value_type;
    };

        /* Generic trait for iterators without type */



	template <class RepIterator, class RowIdxIterator>
	class _RawIndexedIterator
	{
	    public:
// 		// typedef std::pair<size_t, size_t> value_type;
// 		typedef typename RowIdxIterator/*::value_type*/::second_type value_type;
// 		typedef typename RowIdxIterator::value_type::second_type value_type;
            typedef typename IteratorValueType< RowIdxIterator >::value_type::second_type value_type;

		_RawIndexedIterator (size_t idx, const RepIterator &i, const RowIdxIterator &j, const RepIterator &A_end)
			: _i (i), _j (j), _A_end (A_end), _r_index (idx)
		{
			if( _i == _A_end ) return;
			while(_j == _i->end ()) {
				++_r_index;
				if (++_i == _A_end) return;
				_j = _i->begin ();
			}
			_c_index =_j->first;		
		}

		_RawIndexedIterator (const _RawIndexedIterator &iter)
			: _i (iter._i), _j (iter._j), _A_end (iter._A_end), _r_index (iter._r_index), _c_index (iter._c_index)
		{}

		_RawIndexedIterator ()
		{}

		_RawIndexedIterator &operator = (const _RawIndexedIterator &iter) 
		{
			_A_end = iter._A_end;
			_i = iter._i;
			_j = iter._j;
			_r_index = iter._r_index;
			_c_index = iter._c_index;

			return *this;
		}

		bool operator == (const _RawIndexedIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const _RawIndexedIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		_RawIndexedIterator &operator ++ ()
		{
			++_j;
			while(_j == _i->end ()){
				++_r_index;
				if (++_i == _A_end) return *this;
				_j = _i->begin ();
			}
			_c_index = _j->first;	
						
			// if (++_j == _i->end ()) {
// 				if (++_i != _A_end) {
// 					_j = _i->begin ();
// 					++_r_index;
// 				}
// 			}

// 			_c_index = _j->first;

			return *this;
		}

		_RawIndexedIterator operator ++ (int)
		{
			_RawIndexedIterator tmp = *this;
			++(*this);
			return tmp;
		}

		_RawIndexedIterator &operator -- ()
		{
			while (_j == _i->begin ()) {
				_j = (--_i)->end ();
				--_r_index;
			}

			--_j;
			_c_index = _j->first;
			return *this;
		}

		_RawIndexedIterator operator -- (int)
		{
			_RawIndexedIterator tmp = *this;
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

		size_t rowIndex () const
			{ return _r_index; }
		size_t colIndex () const
			{ return _c_index; }
		const value_type &value() const
			{ return _j->second; }

	    private:
		RepIterator _i;
		RowIdxIterator _j;
		RepIterator _A_end;

		mutable size_t _r_index;
		mutable size_t _c_index;
	};

	typedef _RawIndexedIterator<typename Rep::iterator, typename Row::iterator> RawIndexedIterator;
	typedef _RawIndexedIterator<typename Rep::const_iterator, typename Row::const_iterator> ConstRawIndexedIterator;

	RawIndexedIterator rawIndexedBegin ()
	{ return RawIndexedIterator (0, _A.begin (), _A.front ().begin (), _A.end ()); }
	RawIndexedIterator rawIndexedEnd ()
		{ return RawIndexedIterator (_m, _A.end (), _A.back ().end (), _A.end ()); }
	ConstRawIndexedIterator rawIndexedBegin () const
		{ return ConstRawIndexedIterator (0, _A.begin (), _A.front ().begin (), _A.end ()); }
	ConstRawIndexedIterator rawIndexedEnd () const
		{ return ConstRawIndexedIterator (_m, _A.end (), _A.back ().end (), _A.end ()); }

	Row &getRow (size_t i) { return _A[i]; }
	Row &operator [] (size_t i) { return _A[i]; }
	ConstRow &operator [] (size_t i) const { return _A[i]; }

	template <class Vector> Vector &columnDensity (Vector &v) const;
	SparseMatrixBase &transpose (SparseMatrixBase &AT) const;

    protected:

	friend class SparseMatrixWriteHelper<Element, Row>;
	friend class SparseMatrixReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;

    	template<class F, class R, class T> friend class SparseMatrixBase;
};

/* Specialization for sparse associative vectors */

template <class _Element, class _Row>
class SparseMatrixBase<_Element, _Row, VectorCategories::SparseAssociativeVectorTag >
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
        struct rebind
        { typedef SparseMatrixBase<typename _Tp1::Element, _R1, VectorCategories::SparseAssociativeVectorTag> other; };

	SparseMatrixBase (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrixBase (const SparseMatrixBase<Element, Row> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}

    	template<class VectorType>
	SparseMatrixBase (const SparseMatrixBase<Element, VectorType> &A)
		: _A(A.m), _m (A._m), _n (A._n) {
            typename Rep::iterator meit = this->_A.begin();
            typename SparseMatrixBase<Element, VectorType>::Rep::const_iterator copit = A._A.begin();
            for( ; meit != this->_A.end(); ++meit, ++copit)
                LinBox::RawVector<Element>::convert(*meit, *copit);
        }

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrixBase ( MatrixStream<Field>& ms );
	~SparseMatrixBase () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const { 
            size_t s(0);
            for(typename Rep::const_iterator it = _A.begin(); it != _A.end(); ++it)
                s+= LinBox::RawVector<_Element>::size(*it);
            return s;
        }

	template <class Field>
	std::istream &read (std::istream &is, const Field &F, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row>::read
			  (*this, is, F, format); }
	std::istream &read (std::istream &is, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row>::read
			  (*this, is, SparseMatrixReadWriteHelper<Element, Row>::NoField (),
			   format); }
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row>::write
			  (*this, os, F, format); }
	std::ostream &write (std::ostream &os, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row>::write
			  (*this, os, SparseMatrixReadWriteHelper<Element, Row>::NoField (),
			   format); }

	void           setEntry (size_t i, size_t j, const Element &value) { _A[i][j] = value; }
	Element       &refEntry (size_t i, size_t j)                       { return _A[i][j]; }
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const     { return x = _A[i][j];}

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	ConstRowIterator rowBegin () const 
		{ return _A.begin (); }
	ConstRowIterator rowEnd () const
		{ return _A.end (); }
	RowIterator rowBegin ()
		{ return _A.begin (); }
	RowIterator rowEnd ()
		{ return _A.end (); }

	template <class RepIterator, class RowEltIterator, class _I_Element>
	class _RawIterator
	{
	    public:
		typedef _I_Element value_type;

		_RawIterator (const RepIterator &i, const RowEltIterator &j, const RepIterator &A_end)
			: _i (i), _j (j), _A_end (A_end)
		{
			if( _i == _A_end ) return;
			while( _j == _i->end() ) {
				if( ++_i == _A_end ) return;
				_j = _i->begin();
			}
		}

		_RawIterator (const _RawIterator &iter)
			: _i (iter._i), _j (iter._j), _A_end (iter._A_end)
		{}

		_RawIterator () {}

		_RawIterator &operator = (const _RawIterator &iter) 
		{
			_i = iter._i;
			_j = iter._j;
			_A_end = iter._A_end;

			return *this;
		}

		bool operator == (const _RawIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const _RawIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		_RawIterator &operator ++ ()
		{
			while (++_j == _i->end ()) {
				if (++_i == _A_end ()) return *this;
				_j = _i->begin ();
			}
			return *this;
		}

		_RawIterator operator ++ (int)
		{
			_RawIterator tmp = *this;
			++(*this);
			return tmp;
		}

		_RawIterator &operator -- ()
		{
			while (_j == _i->begin ())
				_j = (--_i)->end ();
			--_j;
			return *this;
		}

		_RawIterator operator -- (int)
		{
			_RawIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return _j->second; }
		value_type *operator -> ()
			{ return &(_j->second); }

	    private:
		RepIterator _i;
		RowEltIterator _j;
		RepIterator _A_end;
	};

	typedef _RawIterator<typename Rep::iterator, typename Row::iterator, Element> RawIterator;
	typedef _RawIterator<typename Rep::const_iterator, typename Row::const_iterator, const Element> ConstRawIterator;

	RawIterator rawBegin ()
		{ return RawIterator (_A.begin (), _A.front ().begin (), _A.end ()); }
	RawIterator rawEnd ()
		{ return RawIterator (_A.end (), _A.back ().end (), _A.end ()); }
	ConstRawIterator rawBegin () const
		{ return ConstRawIterator (_A.begin (), _A.front ().begin (), _A.end ()); }
	ConstRawIterator rawEnd () const
		{ return ConstRawIterator (_A.end (), _A.back ().end (), _A.end ()); }

	template <class RepIterator, class RowIdxIterator>
	class _RawIndexedIterator
	{
	    public:
		typedef std::pair<size_t, size_t> value_type;

		_RawIndexedIterator (size_t idx, const RepIterator &i, const RowIdxIterator &j, const RepIterator &A_end)
			: _i (i), _j (j), _A_end (A_end), _r_index (idx), _c_index (0)
		{
			if( _i == _A_end ) return;
			while( _j == _i->end() ) {
				++_r_index;
				if( ++_i == _A_end ) return;
				_j = _i->begin();
			}
			_c_index = _j->second;
		}

		_RawIndexedIterator (const _RawIndexedIterator &iter)
			: _i (iter._i), _j (iter._j), _A_end (iter._A_end), _r_index (iter._r_index), _c_index (iter._c_index)
		{}

		_RawIndexedIterator ()
		{}

		_RawIndexedIterator &operator = (const _RawIndexedIterator &iter) 
		{
			_A_end = iter._A_end;
			_i = iter._i;
			_j = iter._j;
			_r_index = iter._r_index;
			_c_index = iter._c_index;

			return *this;
		}

		bool operator == (const _RawIndexedIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const _RawIndexedIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		_RawIndexedIterator &operator ++ ()
		{
			++_j;
			while (_j == _i->end ()) {
				++_r_index;
				if (++_i == _A_end ()) return *this;
				_j = _i->begin ();
			}

			_c_index = _j->second;

			return *this;
		}

		_RawIndexedIterator operator ++ (int)
		{
			_RawIndexedIterator tmp = *this;
			++(*this);
			return tmp;
		}

		_RawIndexedIterator &operator -- ()
		{
			while (_j == _i->begin ()) {
				_j = (--_i)->end ();
				--_r_index;
			}

			--_j;
			_c_index = _j->second;
			return *this;
		}

		_RawIndexedIterator operator -- (int)
		{
			_RawIndexedIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return *_j; }
		const value_type &operator * () const
			{ return *_j; }
		value_type *operator -> ()
			{ return &(*_j); }
		const value_type *operator -> () const
			{ return &(*_j); }

		size_t rowIndex () const
			{ return _r_index; }
		size_t colIndex () const
			{ return _c_index; }
		const value_type &value() const
			{ return *_j; }

	    private:
		RepIterator _i;
		RowIdxIterator _j;
		RepIterator _A_end;

		mutable size_t _r_index;
		mutable size_t _c_index;
	};

	typedef _RawIndexedIterator<typename Rep::iterator, typename Row::iterator> RawIndexedIterator;
	typedef _RawIndexedIterator<typename Rep::const_iterator, typename Row::const_iterator> ConstRawIndexedIterator;

	RawIndexedIterator rawIndexedBegin ()
		{ return RawIndexedIterator (0, _A.begin (), _A.front ().begin (), _A.end ()); }
	RawIndexedIterator rawIndexedEnd ()
		{ return RawIndexedIterator (_m, _A.end (), _A.back ().end (), _A.end ()); }
	ConstRawIndexedIterator rawIndexedBegin () const
		{ return ConstRawIndexedIterator (0, _A.begin (), _A.front ().begin (), _A.end ()); }
	ConstRawIndexedIterator rawIndexedEnd () const
		{ return ConstRawIndexedIterator (_m, _A.end (), _A.back ().end (), _A.end ()); }

	Row &getRow (size_t i) { return _A[i]; }
	Row &operator [] (size_t i) { return _A[i]; }
	ConstRow &operator [] (size_t i) const { return _A[i]; }

	template <class Vector> Vector &columnDensity (Vector &v) const;
	SparseMatrixBase &transpose (SparseMatrixBase &AT) const;

    protected:

	friend class SparseMatrixWriteHelper<Element, Row>;
	friend class SparseMatrixReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;

    	template<class F, class R, class T> friend class SparseMatrixBase;
};

/* Specialization for sparse parallel vectors */

template <class _Element, class _Row>
class SparseMatrixBase<_Element, _Row, VectorCategories::SparseParallelVectorTag >
{
    public:

	typedef _Element Element;
	typedef _Row Row;
	typedef const Row ConstRow;
	typedef _SP_BB_VECTOR_<Row> Rep;

	template<typename _Tp1, typename _R1 = typename Rebind<_Row,_Tp1>::other >
        struct rebind
        { typedef SparseMatrixBase<typename _Tp1::Element, _R1, VectorCategories::SparseParallelVectorTag> other; };

	SparseMatrixBase (size_t m, size_t n)
		: _A (m), _m (m), _n (n) {}
	SparseMatrixBase (const SparseMatrixBase<Element, Row> &A)
		: _A (A._A), _m (A._m), _n (A._n) {}
    
    	template<class VectorType>
	SparseMatrixBase (const SparseMatrixBase<Element, VectorType> &A)
		: _A(A._m), _m (A._m), _n (A._n) {
            typename Rep::iterator meit = this->_A.begin();
            typename SparseMatrixBase<Element, VectorType>::Rep::const_iterator copit = A._A.begin();
            for( ; meit != this->_A.end(); ++meit, ++copit)
                LinBox::RawVector<Element>::convert(*meit, *copit);
        }

	/** Constructor from a MatrixStream
	 */
	template <class Field>
	SparseMatrixBase ( MatrixStream<Field>& ms );

	~SparseMatrixBase () {}

	size_t rowdim () const { return _m; }
	size_t coldim () const { return _n; }
	size_t size () const { 
            size_t s(0);
            for(typename Rep::const_iterator it = _A.begin(); it != _A.end(); ++it)
                s+= LinBox::RawVector<_Element>::size(*it);
            return s;
        }

	template <class Field>
	std::istream &read (std::istream &is, const Field &F, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row>::read
			  (*this, is, F, format); }
	std::istream &read (std::istream &is, FileFormatTag format = FORMAT_DETECT)
		{ return SparseMatrixReadWriteHelper<Element, Row>::read
			  (*this, is, SparseMatrixReadWriteHelper<Element, Row>::NoField (),
			   format); }
	template <class Field>
	std::ostream &write (std::ostream &os, const Field &F, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row>::write
			  (*this, os, F, format); }
	std::ostream &write (std::ostream &os, FileFormatTag format = FORMAT_PRETTY) const
		{ return SparseMatrixReadWriteHelper<Element, Row>::write
			  (*this, os, SparseMatrixReadWriteHelper<Element, Row>::NoField (),
			   format); }

	void           setEntry (size_t i, size_t j, const Element &value);
	Element       &refEntry (size_t i, size_t j);
	const Element &getEntry (size_t i, size_t j) const;
	Element       &getEntry (Element &x, size_t i, size_t j) const
			{ return x = getEntry (i, j); }

	typedef typename Rep::iterator RowIterator;
	typedef typename Rep::const_iterator ConstRowIterator;

	ConstRowIterator rowBegin () const 
		{ return _A.begin (); }
	ConstRowIterator rowEnd () const
		{ return _A.end (); }
	RowIterator rowBegin ()
		{ return _A.begin (); }
	RowIterator rowEnd ()
		{ return _A.end (); }

	template <class RepIterator, class RowEltIterator, class _I_Element>
	class _RawIterator
	{
	    public:
		typedef _I_Element value_type;

		_RawIterator (const RepIterator &i, const RowEltIterator &j, const RepIterator &A_end)
			: _i (i), _j (j), _A_end (A_end)
		{
			if( _i == _A_end ) return;
			while( _j == _i->second.end() ) {
				if( ++_i == _A_end ) return;
				_j = _i->second.begin();
			}
		}

		_RawIterator (const _RawIterator &iter)
			: _i (iter._i), _j (iter._j), _A_end (iter._A_end)
		{}

		_RawIterator () {}

		_RawIterator &operator = (const _RawIterator &iter)
		{
			_i = iter._i;
			_j = iter._j;
			_A_end = iter._A_end;

			return *this;
		}

		bool operator == (const _RawIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const _RawIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		_RawIterator &operator ++ ()
		{
			++_j;
			while( _j == _i->second.end() ) {
				if( ++_i == _A_end ) return *this;
				_j = _i->second.begin();
			}
			return *this;
		}

		_RawIterator operator ++ (int)
		{
			_RawIterator tmp = *this;
			++(*this);
			return tmp;
		}

		_RawIterator &operator -- ()
		{
			while (_j == _i->second.begin ())
				_j = (--_i)->second.end ();
			--_j;
			return *this;
		}

		_RawIterator operator -- (int)
		{
			_RawIterator tmp = *this;
			--(*this);
			return tmp;
		}

		value_type &operator * ()
			{ return static_cast<value_type&>(*_j); }
		value_type *operator -> ()
			{ return &(*_j); }

	    private:
		RepIterator _i;
		RowEltIterator _j;
		RepIterator _A_end;
	};

	typedef _RawIterator<typename Rep::iterator, typename Row::second_type::iterator, Element> RawIterator;
	typedef _RawIterator<typename Rep::const_iterator, typename Row::second_type::const_iterator, const Element> ConstRawIterator;

	RawIterator rawBegin ()
		{ return RawIterator (_A.begin (), _A.front ().second.begin (), _A.end ()); }
	RawIterator rawEnd ()
		{ return RawIterator (_A.end (), _A.back ().second.end (), _A.end ()); }
	ConstRawIterator rawBegin () const
		{ return ConstRawIterator (_A.begin (), _A.front ().second.begin (), _A.end ()); }
	ConstRawIterator rawEnd () const
		{ return ConstRawIterator (_A.end (), _A.back ().second.end (), _A.end ()); }

	template <class RepIterator, class RowIdxIterator>
	class _RawIndexedIterator
	{
	    public:
		typedef typename RepIterator::value_type::second_type::value_type value_type;

		// Dan Roche 7-6-05 Fixed a seg fault this code was causing
		_RawIndexedIterator (size_t idx, const RepIterator &i, const RowIdxIterator &j, const RepIterator &A_end)
			: _i (i), _j (j), _A_end (A_end), _r_index (idx), _c_index(0), _value_index(0)
		{
			if( _i == _A_end ) return;
			while( _j == _i->first.end() ) {
				if( ++_i == _A_end ) return;
				++_r_index;
				_j = _i->first.begin();
			}
			
			_c_index = *_j;
			_value_index = _j-_i->first.begin();
                }

		_RawIndexedIterator (const _RawIndexedIterator &iter)
			: _i (iter._i), _j (iter._j), _A_end (iter._A_end), _r_index (iter._r_index), _c_index (iter._c_index), _value_index( iter._value_index )
		{}

		_RawIndexedIterator ()
		{}

		_RawIndexedIterator &operator = (const _RawIndexedIterator &iter) 
		{
			_A_end = iter._A_end;
			_i = iter._i;
			_j = iter._j;
			_r_index = iter._r_index;
			_c_index = iter._c_index;
			_value_index = iter._value_index;

			return *this;
		}

		bool operator == (const _RawIndexedIterator &i) const
			{ return (_i == i._i) && (_j == i._j); }

		bool operator != (const _RawIndexedIterator &i) const
			{ return (_i != i._i) || (_j != i._j); }

		_RawIndexedIterator &operator ++ ()
		{
			if(_j != _i->first.end ()) {
				++_j ;
				++_value_index;
			}
			while(_j == _i->first.end ()) {
				++_r_index;
				if (++_i == _A_end) return *this;
                                _j = _i->first.begin ();
				_value_index = 0;
			}
                        _c_index = *_j;

			return *this;
		}

		_RawIndexedIterator operator ++ (int)
		{
			_RawIndexedIterator tmp = *this;
			++(*this);
			return tmp;
		}

		_RawIndexedIterator &operator -- ()
		{
			while (_j == _i->first.begin ()) {
				_j = (--_i)->first.end ();
				_value_index = _i->first.size();
				--_r_index;
			}

			--_j;
			--_value_index;
			_c_index = *_j;
			return *this;
		}

		_RawIndexedIterator operator -- (int)
		{
			_RawIndexedIterator tmp = *this;
			--(*this);
			return tmp;
		}

//  		const value_type &operator * () const
//  		    { return *(_i->second.begin () + _c_index); }
		value_type &operator * ()
		    { return (value_type&)(_i->second)[_value_index]; }
// 		value_type *operator -> ()
// 			{ return &(*(_i->second.begin () + _c_index)); }
//  		const value_type *operator -> () const
//  			{ return &(*(_i->second.begin () + _c_index)); }

		size_t rowIndex () const
			{ return _r_index; }
		size_t colIndex () const
			{ return _c_index; }
 		const value_type &value () const 
                { return (value_type&)(_i->second)[_value_index]; }

	    private:
		RepIterator _i;
		RowIdxIterator _j;
		RepIterator _A_end;

		mutable size_t _r_index;
		mutable size_t _c_index;
		mutable size_t _value_index;
	};

	typedef _RawIndexedIterator<typename Rep::iterator, typename Row::first_type::iterator> RawIndexedIterator;
	typedef _RawIndexedIterator<typename Rep::const_iterator, typename Row::first_type::const_iterator> ConstRawIndexedIterator;

	RawIndexedIterator rawIndexedBegin ()
		{ return RawIndexedIterator (0, _A.begin (), _A.front ().first.begin (), _A.end ()); }
	RawIndexedIterator rawIndexedEnd ()
		{ return RawIndexedIterator (_m, _A.end (), _A.back ().first.end (), _A.end ()); }
	ConstRawIndexedIterator rawIndexedBegin () const
		{ return ConstRawIndexedIterator (0, _A.begin (), _A.front ().first.begin (), _A.end ()); }
	ConstRawIndexedIterator rawIndexedEnd () const
		{ return ConstRawIndexedIterator (_m, _A.end (), _A.back ().first.end (), _A.end ()); }

	Row &getRow (size_t i) { return _A[i]; }
	Row &operator [] (size_t i) { return _A[i]; }
	ConstRow &operator [] (size_t i) const { return _A[i]; }

	template <class Vector> Vector &columnDensity (Vector &v) const;
	SparseMatrixBase &transpose (SparseMatrixBase &AT) const;

    protected:

	friend class SparseMatrixWriteHelper<Element, Row>;
	friend class SparseMatrixReadWriteHelper<Element, Row>;

	Rep               _A;
	size_t            _m;
	size_t            _n;

    	template<class F, class R, class T> friend class SparseMatrixBase;
};

template <class Element, class Row>
std::ostream &operator << (std::ostream &os, const SparseMatrixBase<Element, Row> &A)
	{ return A.write (os); }

template <class Element, class Row>
std::istream &operator >> (std::istream &is, SparseMatrixBase<Element, Row> &A)
	{ return A.read (is); }

template <class Element, class Row, class Trait>
struct MatrixTraits< SparseMatrixBase<Element, Row, Trait> >
{ 
	typedef SparseMatrixBase<Element, Row, Trait> MatrixType;
	typedef typename MatrixCategories::RowMatrixTag MatrixCategory; 
};

template <class Element, class Row, class Trait>
struct MatrixTraits< const SparseMatrixBase<Element, Row, Trait> >
{ 
	typedef const SparseMatrixBase<Element, Row, Trait> MatrixType;
	typedef typename MatrixCategories::RowMatrixTag MatrixCategory; 
};

} // namespace LinBox

#include "linbox/matrix/sparse.inl"

#endif // __LINBOX_matrix_sparse_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
