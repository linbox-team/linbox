/* -*- mode: c; style: linux -*- */

/* linbox/src/blackbox/sparse-matrix-base.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#ifndef __SPARSE_MATRIX_BASE_H
#define __SPARSE_MATRIX_BASE_H

#include <vector>    // STL vector
#include <utility>   // STL pair
#include <iostream>
#include <algorithm>

#include "linbox/vector/vector-traits.h"
#include "linbox/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{
	/** Auxillary sparse matrix base template.
	 * This is a class of sparse matrices templatized by the field in
	 * which the Elements reside.  The matrix itself is stored as an
	 * STL vector of \Ref{LinBox} sparse vectors of integers and field Elements.
	 * Each sparse vector corresponds to one row of the matrix, 
	 * and each pair (j, a) in sparse vector i corresponds to the (i,j) 
	 * entry of the matrix.
	 *
	 * It is templatized by the \Ref{LinBox} field in which the arithmetic
	 * is done, the sparse vector type with which the rows of the matrix 
	 * are implemented, and the vector category of the row implementation.
	 * This third template parameter is defaulted to be the \Ref{LinBox} 
	 * vector trait of the vector.  This class is then specialized for 
	 * sequence and associative sparse vectors.
	 *
	 * This class does not contain any functions operating on vectors.
	 * These functions must all be implemented in derived classes which must
	 * be able to access the actual data structes stored in this class.  
	 * This separation occurs to avoid implementing all methods twice.  
	 * Now, only the methods that depend on the vectors are implemented twice.
	 *
	 * This class is originally from the Programming Languages for Mathematicians
	 * class taught by Erich Kaltofen at NCSU in Fall 1998.
	 * @param Field LinBox field
	 * @param Row   LinBox sparse vector type to use for rows of matrix
	 */
	template <class Field, class Row, class Trait = VectorTraits<Row>::VectorCategory>
	class SparseMatrixBase
	{
		public:

	        /// Element type
	        typedef typename Field::Element Element;
 
		/// Row iterator type
		typedef typename Row::iterator RowIterator;

		/// Constant row iterator
		typedef typename Row::const_iterator ConstRowIterator;

		/** Constructor.
		 * Note: the copy constructor and operator= will work as intended
		 *       because of STL's container design
		 * @param  F  the field of entries; passed so that a possible paramter 
		 * 	   such as a modulus is known to the matrix.
		 * @param  m  row dimension
		 * @param  n  column dimension
		 */
		SparseMatrixBase (const Field& F, size_t m, size_t n);

		/** Destructor. */
		~SparseMatrixBase () {}

		/** Retreive row dimensions of Sparsemat matrix.
		 * @return integer number of rows of SparseMatrixBase matrix.
		 */
		size_t get_rowdim (void) const { return _m; }

		/** Retreive column dimensions of Sparsemat matrix.
		 * @return integer number of columns of SparseMatrixBase matrix.
		 */
		size_t get_coldim (void) const { return _n; }

		/** Retrieve matrix Element.
		 * If the indices are out of range, an error is printed and
		 * the zero Element is returned.
		 * Unlike STL vector's operator[] this member does not
		 * return a reference to the Field Element in the matrix
		 * it is therefore impossible to write
		 * A[make_pair (3, 5)] = Element (1);
		 * the reason is that zero entries have no memory where
		 * one could place the right side Field Element
		 * @return     a copy of the Element in row i, column j
		 * @param      ind pair of indices (i,j)
		 */
		Element operator[] (const pair<size_t, size_t>& ind) const;

		/** Insert matrix Element.
		 * sets A[ind.first, ind.second] = a.
		 * For example, A.put_value (make_pair (3, 5), Element (-4))
		 * If the indices are out of range, an error is printed.
		 * If the Element attempting to be inserted is the zero Element,
		 * no cell will be inserted and any already existing cell for the
		 * entry will be erased.
		 * @param ind  pair of integers (i,j) for row i and column j
		 * @param a    field Element to insert in matrix
		 */
		void put_value (const pair<size_t, size_t>& ind, const Element& a);

		/** Print matrix.
		 * Prints rows as lists.
		 * Can be used in operator <<.
		 * @param  os  output stream on which to print matrix.
		 */
		ostream& write (ostream& os) const;

		/** Read matrix.
		 * can be called in operator>>.
		 * Works by reading integers for row index and column index
		 * and then Element to insert in matrix.
		 * Stops reading upon non-positive row index or end of file.
		 * @param  is  input stream from which to read matrix.
		 */
		istream& read (istream& is);

		/** Exchange two rows.
		 * Exchanges rows i and j of the matrix.
		 * @param  i first row index
		 * @param  j second row index
		 */
		void swaprow (size_t i, size_t j);

		/** Adds multiple of one row to another.
		 * Adds a*(row i) to (row j).
		 * @param  i first row index
		 * @param  j second row index
		 * @param  a multiple of row i to add to row j
		 */
		void addrow (size_t i, size_t j, const Element& a);

	    protected:

		/* Sparse matrix data structure.
		 * _A[i], the i-th row, is a LinBox sparse vector
		 * Each sparse vector entry is a column index/value pair.
		 */
		std::vector< Row > _A;

		// Field used for all arithmetic
		Field _F;

		// Row dimension from 1..m; the actual dimensions
		size_t _m;

		// Column dimension from 1..n; the actual dimensions
		size_t _n;

	};// end template class SparseMatrixBase

	// Specialization of SparseMatrixBase for sequence rows
	template <class Field, class Row>
	class SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
	{
	    public:
		
		typedef typename Field::Element Element;
		typedef typename Row::iterator RowIterator;
		typedef typename Row::const_iterator ConstRowIterator;

		SparseMatrixBase (const Field& F, size_t m, size_t n);
		~SparseMatrixBase () {}
		size_t get_rowdim (void) const { return _m; }
		size_t get_coldim (void) const { return _n; }
		Element operator[] (const pair<size_t, size_t>& ind) const;
		void put_value (const pair<size_t, size_t>& ind, const Element& a);
		ostream& write (ostream& os) const;
		istream& read (istream& is);
		void swaprow (size_t i, size_t j);
		void addrow (size_t i, size_t j, const Element& a);

	    protected:

		std::vector< Row > _A;
		Field _F;
		size_t _m;
		size_t _n;

	    private:
		// used in lower_bound as function object
		struct comp_w_col_index 
		{
			bool operator ()
			(const pair< size_t, Element >& entry, size_t col_in)

				{ return entry.first < col_in; }
		}; // struct comp_w_col_index
	};// class SparseMatrixBase<SparseSequenceVectorTag>

	// Specialization of SparseMatrixBase for associative rows
	template <class Field, class Row>
		class SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
	{
	    public:

		typedef typename Field::Element Element;
		typedef typename Row::iterator RowIterator;
		typedef typename Row::const_iterator ConstRowIterator;

		SparseMatrixBase (const Field& F, size_t m, size_t n);
		~SparseMatrixBase () {}
		size_t get_rowdim (void) const { return _m; }
		size_t get_coldim (void) const { return _n; }
		Element operator[] (const pair<size_t, size_t>& ind) const;
		void put_value (const pair<size_t, size_t>& ind, const Element& a);
		ostream& write (ostream& os) const;
		istream& read (istream& is);
		void swaprow (size_t i, size_t j);
		void addrow (size_t i, size_t j, const Element& a);

	    protected:

		std::vector< Row > _A;
		Field _F;
		size_t _m;
		size_t _n;
	};// class SparseMatrixBase<SparseAssociativeVectorTag>

	// Implementation of matrix methods for sparse sequence rows

	template <class Field, class Row>
	inline SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
		::SparseMatrixBase (const Field& F, size_t m, size_t n) 
		: // constructs a matrix of the given dimensions with all 0s
		_A (m, Row () ),
		_F (F),
		_m (m), _n (n) // set the dimensions
		{}

	template <class Field, class Row>
	void SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
		::put_value (const pair<size_t, size_t>& ind, const Element& a) 
	{
		size_t i = ind.first;
		size_t j = ind.second;
		RowIterator iter;
		bool found (true);

		// Check row and column indices and print error if they are out of range.
		if ( (i >= _m) || (i < 0) || (j >= _n) || (j < 0) ) {
			cerr << endl
			     << "ERROR:  Element indices exceed matrix dimensions." << endl
			     << "	 Element not inserted in matrix." << endl << endl;
			return;
		}

		// Find appropriate location of Element in sparse vector.
		if ( (_A[i]).begin () == (_A[i]).end () )
			iter = (_A[i]).end ();
		else
			iter = lower_bound ( (_A[i]).begin (), 
					     (_A[i]).end (), 
					     j,
					     comp_w_col_index () );

		// Check to see if Element already exists.
		if ( (_A[i]).end () == iter )
			found = false;
		else
			if ( iter->first != j )
				found = false;


		// If Element is already in row, replace old value with new.
		// Otherwise, insert the Element in the row.
		if (found) 
		{
			if (_F.isZero (a))
				_A[i].erase (iter);
			else
				iter->second = a;
		}
		else
			if (!_F.isZero (a))
				_A[i].insert (iter, make_pair (j,a));

	} // void SparseMatrixBase<SparseSequenceVectorTag>::put_value (...)

	template <class Field, class Row>
	typename Field::Element SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
		::operator[] (const pair<size_t, size_t>& ind) const
	{
		Element zero;

		_F.init (zero, 0);

		size_t i = ind.first;
		size_t j = ind.second;

		// Check row and column indices and print error if they are out of range.
		if ( (i >= _m) || (i < 0) || (j >= _n) || (j < 0) ) {
			cerr << endl
			     << "ERROR:  Element indices exceed matrix dimensions." << endl
			     << endl;
			return zero;
		}

		Row row (_A[i]);
		RowIterator iter;
		bool found (true);

		// Find appropriate location of Element in row.
		if ( row.begin () == row.end () )
			iter = row.end ();
		else
			iter = lower_bound ( row.begin (), row.end (), j, comp_w_col_index () );
		// Check to see if Element exists.
		if ( row.end () == iter )
			found = false;
		else
			if ( iter->first != j )
				found = false;

		// If Element is found, return non-zero value.
		// Otherwise value is zero.
		if (found)
			return iter->second;
		else
			return zero;

	} // Element SparseMatrixBase<SparseSequenceVectorTag>::operator[] (...) const

	template <class Field, class Row>
	inline ostream &SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
		::write (ostream& os) const
	{
		for (size_t i = 0; i <= _m - 1; i++) {
			os << "Row " << i << ": ";

			for (ConstRowIterator iter_i = _A[i].begin (); 
			     iter_i != _A[i].end (); 
			     iter_i++)
			{
				os << "(" << (*iter_i).first << ", ";
				_F.write (os, (*iter_i).second);
				os << ")";
			}
			os << endl;
		}
 
		return os;

	} // ostream& SparseMatrixBase<SparseSequenceVectorTag>::write (...) const

	template <class Field, class Row>
	inline istream& SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
		::read (istream& is)
	{
		size_t i, j;
		Element el;

		_F.init (el, 0);

		while (is >> i) {
			// operator>> returns a reference to an istream
			// then istream::operator void*() is called which
			// returns !basic_ios::fail () [Stroutstrup, p.617, ???]

			if (i == size_t (-1)) break; // return also if row index is -1
			is >> j;
			_F.read (is, el);
			put_value (make_pair (i,j), el);
		}

		return is;

	} // istream& SparseMatrixBase<SparseSequenceVectorTag>::read (...)

	template <class Field, class Row>
	void SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
		::swaprow (size_t i, size_t j) 
	{
		// exchanges row i and j in A
		// note:	  uses row::swap

		// Check row indices and print error if they are out of range.
		if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) {
			cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
			     << "        No rows exchanged." << endl << endl;
			return;
		}

		// Swap rows i and j using row::swap
		_A[i].swap ( _A[j] );

		return;
	} // void SparseMatrixBase<SparseSequenceVectorTag>::swaprow (...)

	/* This implementation is works for lists and deques, but it causes 
	 * segmentation faults for vectors.  Inserting Elements into a vector
	 * invalidates iterators *before* the inserted Element, which is contrary 
	 * to the standard.
	 * 
	 * This may not be a good implementation, anyway, because according
	 * to the C++ standard, insertion and erasure can invalidate all iterators 
	 * and Element references to the sequence.
	 *
	 */
#if 0
	template <class Field, class Row>
	void SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
		::addrow (size_t i, size_t j,const Element& a) 
	{
		// Check row indices and print error if they are out of range.
		if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) {
			cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
			     << "        No row addition preformed." << endl << endl;
			return;
		}

		// Check to see if a is the zero Field Element.
		// If so, no addition is performed.
		if (_F.isZero (a)) return;

		// Check to see if row i is empty.  If so, no addition is preformed.
		if ( (_A[i]).begin () == (_A[i]).end () ) return;

		size_t k;
		Element value;
		_F.init (value, 0);

		LinBox::faxpy<Field> Faxpy (_F, a);

		// iterators to point to place in rows i and j respectively,
		// and extra iterator for erasing from row j;
		RowIterator iter_i, iter_j, iter;

		bool found (true);

		iter_j = (_A[j]).begin (); // start at beginning of second row

		// iterate over Elements in row i
		for ( iter_i = (_A[i]).begin (); iter_i != (_A[i]).end (); iter_i++) {
			found = true;
			k = iter_i->first;  // marks current column.

			// Find where column k occurs in row j.
			while ( ( (_A[j]).end () != iter_j ) && ( iter_j->first < k ) )
				iter_j++;

			// Check if row j has Element for column k.
			if ( ( (_A[j]).end () == iter_j ) || ( iter_j->first != k ) )
				found = false;

			// If row j contains Element for column k, perform sum.
			// Otherwise, sum = a * _A[i,k]
			if (found) {
				if (_F.isZero (Faxpy.applyin (iter_j->second, iter_i->second))) {
					iter = iter_j++;
					_A[j].erase (iter);
				} else
					iter_j++;

			} else
				_A[j].insert (iter_j,
					      make_pair (k, _F.mul (value,a,iter_i->second)));
		}

		return;
	} // void SparseMatrixBase<SparseSequenceVectorTag>::addrow (...)
#endif
  
	/* This implementation works for vectors because it avoids the insert 
	 * method.  It creates a new row, using push_back insert new Elements, and 
	 * then copies it into the _A[j] at the end.  This is less efficient than  
	 * doing an inplace row add like above, but no iterators are invalidated 
	 * through the insert methods.
	 */
	template <class Field, class Row>
	void SparseMatrixBase<Field, Row, VectorCategories::SparseSequenceVectorTag>
		::addrow (size_t i, size_t j,const Element& a) 
	{
		// Check row indices and print error if they are out of range.
		if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) {
			cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
			     << "        No row addition preformed." << endl << endl;
			return;
		}

		// Check to see if a is the zero Field Element.
		// If so, no addition is performed.
		if (_F.isZero (a)) return;

		// Check to see if row i is empty.  If so, no addition is preformed.
		if ( (_A[i]).begin () == (_A[i]).end () ) return;

		// variables used in computation
		Element value;

		_F.init (value, 0);

		Row row;
		RowIterator iter_i, iter_j (_A[j].begin ());

		for (iter_i = _A[i].begin (); iter_i != _A[i].end (); iter_i++) {
			while ( (iter_j != _A[j].end ()) && (iter_j->first < iter_i->first) ) {
				row.push_back (*iter_j);
				iter_j++;
			}

			if ( (iter_j != _A[j].end ()) && (iter_j->first == iter_i->first) ) {
				if (!_F.isZero (_F.axpy (value, a, iter_i->second, iter_j->second)))
					row.push_back (make_pair (iter_i->first, value));
	
				iter_j++;
			}
			else
				row.push_back (make_pair (iter_i->first,
							  _F.mul (value, a, iter_i->second)));
		}

		while (iter_j != _A[j].end ()) {
			row.push_back (*iter_j);
			iter_j++;
		}
    
		_A[j] = row;

		return;
	} //  void SparseMatrixBase<SparseSequenceVectorTag>::addrow (...)

	// Implementation of matrix methods for sparse associative rows

	template <class Field, class Row>
	inline SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
		::SparseMatrixBase (const Field& F, size_t m, size_t n) 
		: // constructs a matrix of the given dimensions with all 0s
		_A (m, Row () ),
		_F (F),
		_m (m), _n (n) // set the dimensions
		{}

	template <class Field, class Row>
	void SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
		::put_value (const pair<size_t, size_t>& ind, const Element& a) 
	{
		size_t i = ind.first;
		size_t j = ind.second;
		RowIterator iter;

		// Check row and column indices and print error if they are out of range.
		linbox_check ( (i >= _m) || (i < 0) || (j >= _n) || (j < 0) );

		// Find Element in map.  
		// If exists, replace value if not zero, or remove if value is zero.
		// If not found, insert non-zero Element
		if ( (iter = _A[i].find (j)) != _A[i].end () ) {
			if (_F.isZero (a))
				_A[i].erase (iter);
			else
				iter->second = a;
		} else
			if (!_F.isZero (a))
				_A[i].insert (make_pair (j, a));
	} // void SparseMatrixBase<SparseAssociativeVectorTag>::put_value (...)

	template <class Field, class Row>
	typename Field::Element SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
		::operator[] (const pair<size_t, size_t>& ind) const
	{
		Element zero;

		_F.init (zero, 0);

		size_t i = ind.first;
		size_t j = ind.second;

		linbox_check ( (i >= _m) || (i < 0) || (j >= _n) || (j < 0) );

		ConstRowIterator iter;

		if ( (iter = _A[i].find (j)) != _A[i].end () )
			return iter->second;
		else
			return zero;
	} // Element SparseMatrixBase<SparseAssociativeVectorTag>::operator[] (...)

	template <class Field, class Row>
	inline ostream &SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
		::write (ostream& os) const
	{
		for (size_t i = 0; i <= _m - 1; i++) {
			os << "Row " << i << ": ";

			for (ConstRowIterator iter_i = _A[i].begin (); 
			     iter_i != _A[i].end (); 
			     iter_i++)
			{
				os << "(" << iter_i->first << ", ";
				_F.write (os, iter_i->second);
				os << ")";
			}
			os << endl;
		}
 
		return os;

	} // ostream& SparseMatrixBase<SparseAssociativeVectorTag>::write (...) const

	template <class Field, class Row>
	inline istream &SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
		::read (istream& is)
	{
		size_t i, j;
		Element el;

		_F.init (el, 0);

		while (is >> i)
			// operator>> returns a reference to an istream
			// then istream::operator void*() is called which
			// returns !basic_ios::fail () [Stroutstrup, p.617, ???]
		{
			if (i == size_t (-1)) break; // return also if row index is -1
			is >> j;
			_F.read (is, el);
			put_value (make_pair (i,j), el);
		} // while-loop

		return is;

	} // istream& SparseMatrixBase<SparseAssociativeVectorTag>::read (...)

	template <class Field, class Row>
	void SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
		::swaprow (size_t i, size_t j) 
	{
		// exchanges row i and j in A
		// note:	  uses row::swap

		// Check row indices and print error if they are out of range.
		if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) 
		{
			cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
			     << "        No rows exchanged." << endl << endl;
			return;
		}

		// Swap rows i and j using row::swap
		_A[i].swap ( _A[j] );

		return;
	} // void SparseMatrixBase<SparseAssociativeVectorTag>::swaprow (...)

	/* This implementation is an inplace row addition, but it isn't clear 
	 * from the standard if insertion and deletion invalidates any iterators
	 * and Element references other than the obvious ones refering to a
	 * deleted entries.
	 *
	 */
#if 0
	template <class Field, class Row>
	void SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
		::addrow (size_t i, size_t j,const Element& a) 
	{
		// Check row indices and print error if they are out of range.
		if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) {
			cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
			     << "        No row addition preformed." << endl << endl;
			return;
		}

		// Check to see if a is the zero Field Element.
		// If so, no addition is performed.
		if (_F.isZero (a)) return;

		// Check to see if row i is empty.  If so, no addition is preformed.
		if ( (_A[i]).begin () == (_A[i]).end () ) return;

		size_t k;
		Element value;

		_F.init (value, 0);

		LinBox::faxpy<Field> Faxpy (_F, a);
		// iterators to point to place in rows i and j respectively,
		// and extra iterator for erasing from row j;
		RowIterator iter_i, iter_j;

		bool found (true);

		iter_j = (_A[j]).begin (); // start at beginning of second row

		// iterate over Elements in row i
		for ( iter_i = (_A[i]).begin (); iter_i != (_A[i]).end (); iter_i++) {
			found = true;
			k = iter_i->first;  // marks current column.

			// Find where column k occurs in row j.
			while ( ( (_A[j]).end () != iter_j ) && ( iter_j->first < k ) )
				iter_j++;

			// Check if row j has Element for column k.
			if ( ( (_A[j]).end () == iter_j ) || ( iter_j->first != k ) )
				found = false;

			// If row j contains Element for column k, perform sum.
			// Otherwise, sum = a * _A[i,k]
			if (found) {
				if (_F.isZero (Faxpy.applyin (iter_j->second, iter_i->second)))
					_A[j].erase (iter_j++);
				else
					iter_j++;

			} else
				_A[j].insert (iter_j,
					      make_pair (k, _F.mul (value,a,iter_i->second)));
		}
	} // void SparseMatrixBase<SparseAssociativeVectorTag>::addrow (...)
#endif
	/* This implementation does not have to worry about any invalidated
	 * iterators and references because the addition is not done inplace.
	 * However, this means it is not as efficient since a new row has
	 * to be created and then assigned to _A[j].
	 */
	template <class Field, class Row>
	void SparseMatrixBase<Field, Row, VectorCategories::SparseAssociativeVectorTag>
		::addrow (size_t i, size_t j,const Element& a) 
	{
		// Check row indices and print error if they are out of range.
		if ( (i >= _m) || (i < 0) || (j >= _m) || (j < 0) ) {
			cerr << endl << "ERROR:  Row indices exceed matrix dimensions." << endl
			     << "        No row addition preformed." << endl << endl;
			return;
		}

		// Check to see if a is the zero Field Element.
		// If so, no addition is performed.
		if (_F.isZero (a)) return;

		// Check to see if row i is empty.  If so, no addition is preformed.
		if ( (_A[i]).begin () == (_A[i]).end () ) return;

		// variables used in computation
		Element value;

		_F.init (value, 0);

		Row row;
		RowIterator iter_i, iter_j (_A[j].begin ());

		for (iter_i = _A[i].begin (); iter_i != _A[i].end (); iter_i++) {
			while ( (iter_j != _A[j].end ()) && (iter_j->first < iter_i->first) ) {
				row.insert (*iter_j);
				iter_j++;
			}

			if ( (iter_j != _A[j].end ()) && (iter_j->first == iter_i->first) ) {
				if (!_F.isZero (_F.axpy (value, a, iter_i->second, iter_j->second)))
					row.insert (make_pair (iter_i->first, value));
	
				iter_j++;
			} else
				row.insert (make_pair (iter_i->first,
						       _F.mul (value, a, iter_i->second)));
		}

		while (iter_j != _A[j].end ()) {
			row.insert (*iter_j);
			iter_j++;
		}
    
		_A[j] = row;
	} //  void SparseMatrixBase<SparseAssociativeVectorTag>::addrow (...)

	// Input/Output Operators.

	template <class Element>
	ostream& operator<<(ostream& os, pair< size_t, Element > entry)
	{
		// Requires operator<<(ostream& Element) which may not be provided.
		os << "(" << entry.first << ", " << entry.second << ")";
		return os;
	}

	template <class Field, class Row>
	ostream& operator<<(ostream& os, const SparseMatrixBase<Field, Row>& A)
		{ return A.write (os); }


	template <class Field, class Row>
	istream& operator>>(istream& is, SparseMatrixBase<Field, Row>& A)
		{ return A.read (is); }

	/* Creates new sparse matrix.
	 * Reads matrix from input stream.
	 * If TRACE is defined, output includes all inputs.
	 * @return pointer to new matrix in dynamic memory
	 * @param  F field in which all arithmetic is done
	 * @param  m row dimensions of matrix (defualt = 0)
	 * @param  n column dimensions of matrix (default = 0)
	 * @param  prompt  boolean for whether to prompt user for input
	 *		      (default = true)
	 * @param  is  istream from which to read input (default = cin)
	 * @param  os  output stream to which to print output (default = cout)
	 * @see SparseMatrixBase
	 */
	template <class Field, class Row>
	SparseMatrixBase<Field, Row> *newSparsemat (const Field &F, 
						    size_t       m      = 0, 
						    size_t       n      = 0,
						    bool         prompt = true,
						    istream     &is     = cin, 
						    ostream     &os     = cout)
	{
		while ( (m <= 0) || (n <= 0) ) {
			if (prompt)
				cout << "What are the matrix's row and column dimenstions? ";
 
			is >> m >> n;

#ifdef TRACE
			os << endl << "The matrix has " << m << " rows and " << n << endl;
#endif
		}
 
		SparseMatrixBase<Field, Row>* A_ptr = new SparseMatrixBase<Field, Row>(F,m,n);

		if (prompt)
			cout << endl << "Input sparse matrix by entering row index, column"
			     << endl << "index, and value.  Remember the matrix is indexed"
			     << endl << "starting at 0.  End with a row index of -1."
			     << endl;

		is >> (*A_ptr);

#ifdef TRACE
		os << endl << "The matrix contains the following elements: "
		   << endl << (*A_ptr) << endl;
#endif

		return A_ptr;
	} // newSparsemat ()
} // namespace LinBox

#endif // __SPARSE_MATRIX_BASE_H
