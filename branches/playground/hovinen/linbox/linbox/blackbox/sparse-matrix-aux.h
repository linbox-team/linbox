/* -*- mode: c; style: linux -*- */

/* linbox/src/blackbox/sparse-matrix-aux.h
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

#ifndef __SPARSE_MATRIX_AUX_H
#define __SPARSE_MATRIX_AUX_H

#include <vector>    // STL vector
#include <utility>   // STL pair
#include <iostream>
#include <algorithm>

#include "linbox/blackbox/sparse-matrix-base.h"
#include "linbox/vector/vector-traits.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Auxillary sparse matrix template.
	 * This class is derived from \Ref{sparsemat_base} which stores the actual
	 * matrix implementation.  This class adds the functions gauss, linsolve, 
	 * and apply to the sparse matrix through two new template parameters.  The
	 * second of these is chosen be default to be the \Ref{LinBox} vector trait
	 * of the vector.  This class is then specialized for dense and sparse 
	 * vectors.  The separation of these methods from the rest of the class is 
	 * done to avoid implementing all methods twice.  Now, only the methods
	 * that depend on the vectors are implemented three times.
	 *
	 * The default class is not implemented.  It's functions should never
	 * be called because partial template specialization should always be
	 * done on the vector traits.
	 *
	 * This class is originally from the Programming Languages for Mathematicians
	 * class taught by Erich Kaltofen at NCSU in Fall 1998.
	 * @param Field LinBox field
	 * @param Row   LinBox sparse vector type to use for rows of matrix
	 * @param Vector LinBox dense or sparse vector
	 * @param Trait  Marker whether to use dense or sparse LinBox vector 
	 *               implementation.  This is chosen by a default parameter 
	 *               and partial template specialization.
	 */
	template <class Field, class Row, class Vector, class Trait = VectorTraits<Vector>::VectorCategory>
	class SparseMatrixAux : public SparseMatrixBase<Field, Row>
	{
	    public:

		/** Constructor.
		 * Note: the copy constructor and operator= will work as intended
		 *       because of STL's container design
		 * @param  F  the field of entries; passed so that a possible paramter 
		 * 	   such as a modulus is known to the matrix.
		 * @param  m  row dimension
		 * @param  n  column dimension
		 */
		SparseMatrixAux (const Field& F, size_t m, size_t n); 

		/** Constructor from sparsemat_base.
		 * @param  A sparsemat_base object of same Field and Row types.
		 */
		SparseMatrixAux (const SparseMatrixBase<Field, Row>& B); 

		/** Destructor. */
		~SparseMatrixAux (); // {}

 		/** Linear Solver.
		 * Solves linear system Ax = b by performing Gaussian elimination
		 * and back substitution.  If the system does not have a uinque
		 * solution, an error is printed and a zero matrix is returned.
		 * Templatized by LinBox vector class.
		 * @param	   b must be encapsulated vector of size m.
		 * @return     encapsualted vector of size n such that A * x = b
		 */
		Vector& linsolve (Vector& b);

		/** Gaussian elimination.
		 * Performs gaussian elimination of matrix A and corresponding
		 * right side vector b.
		 * Templatized by LinBox vector class.
		 * @param  b corresponding right side vector b.  
		 *         (Default is empty vector.)
		 * @return boolean true if succesful, false if not
		 * @see www.math.ncsu.edu/~kaltofen/courses/LinAlgebra/Maple/refpkg/Src/ref.mpl
		 * 	  for a Maple row echelon form algorithm
		 */
		bool gauss (Vector& b = Vector () );
 
		/** Application of BlackBox matrix.
		 * y = A*x.
		 * Templatized by LinBox vector class.
		 * @return reference to output vector y
		 * @param  x input vector
		 */
		Vector& apply (Vector& y, const Vector& x) const;

		/** Application of BlackBox matrix transpose.
		 * y = transpose (A)*x.
		 * Templatized by LinBox vector class.
		 * @return reference to output vector y
		 * @param  x input vector
		 */
		Vector& applyTranspose (Vector& y, const Vector& x) const;
  
	}; // SparseMatrixAux

	// Specialization of SparseMatrixAux for LinBox dense vectors
	template <class Field, class Row, class Vector>
	class SparseMatrixAux<Field, Row, Vector, VectorCategories::DenseVectorTag>
		: public SparseMatrixBase<Field, Row>
	{
	    public:

		SparseMatrixAux (const Field& F, size_t m, size_t n) 
			: SparseMatrixBase<Field, Row>(F, m, n) {}
		SparseMatrixAux (const SparseMatrixBase<Field, Row>& B)
			: SparseMatrixBase<Field, Row>(B) {}
		~SparseMatrixAux () {}
		Vector& linsolve (Vector &x, const Vector &b);
		bool gauss (Vector& b = Vector () );
		Vector& apply (Vector& y, const Vector& x) const;
		Vector& applyTranspose (Vector& y, const Vector& x) const;

	};// SparseMatrixAux<DenseVectorTag>
	  
	// Specialization of SparseMatrixAux for LinBox sparse sequence vectors
	template <class Field, class Row, class Vector>
	class SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseSequenceVectorTag>
		: public SparseMatrixBase<Field, Row>
	{
	    public:

		SparseMatrixAux (const Field& F, size_t m, size_t n) 
			: SparseMatrixBase<Field, Row>(F, m, n) {}
		SparseMatrixAux (const SparseMatrixBase<Field, Row>& B)
			: SparseMatrixBase<Field, Row>(B) {}
		~SparseMatrixAux () {}
		Vector& linsolve (Vector &x, const Vector& b);
		bool gauss (Vector& b = Vector () );
		Vector& apply (Vector& y, const Vector& x) const;
		Vector& applyTranspose (Vector& y, const Vector& x) const;

	    private:
		// used in lower_bound as function object
		struct comp_w_index 
		{
			bool operator ()
			(const pair< size_t, Element >& entry, size_t col_in)

				{ return entry.first < col_in; }
		}; // struct comp_w_index

	};// SparseMatrixAux<SparseSequenceVectorTag>
	  
	// Specialization of SparseMatrixAux for LinBox sparse associative vectors
	template <class Field, class Row, class Vector>
	class SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseAssociativeVectorTag>
		: public SparseMatrixBase<Field, Row>
	{
	    public:

		SparseMatrixAux (const Field& F, size_t m, size_t n) 
			: SparseMatrixBase<Field, Row>(F, m, n) {}
		SparseMatrixAux (const SparseMatrixBase<Field, Row>& B)
			: SparseMatrixBase<Field, Row>(B) {}
		~SparseMatrixAux () {}
		Vector& linsolve (Vector &x, const Vector& b);
		bool gauss (Vector& b = Vector () );
		Vector& apply (Vector& y, const Vector& x) const;
		Vector& applyTranspose (Vector& y, const Vector& x) const;

	};// SparseMatrixAux<SparseAssociativeVectorTag>

	// Implementation of matrix methods for dense vectors

	template <class Field, class Row, class Vector>
	Vector &SparseMatrixAux<Field, Row, Vector, VectorCategories::DenseVectorTag>
		::linsolve (Vector &x, const Vector& b)
	{
		// exceptions: if the system does not have a unique solution,
		//		 an error is printed and a zero vector returned
		size_t bs = b.size ();
    
		// Check to see if dimensions are compatible for solving the system.
		linbox_check (bs != _m);

		// make copies so can restore originals
		std::vector< Row > _AA (_A);
		Vector bb (b);

		// perform gaussian elimination
		if (!gauss (bb)) {
			cerr << endl
			     << "ERROR:  Gaussian elimination not performed." << endl << endl;
			return x;
		}

		// Back substitution to construct solution

		size_t i, j, k;

		// find first pivot element from bottom by scanning rows
		for (i = _m - 1; i >= 0; i--)
			if ( (_A[i]).begin () != (_A[i]).end () ) {
				j = _A[i].begin ()->first;
				break;
			}

		// Check for inconsistancy
		for (k = i+1; k < _m; k++)
			if ( !_F.isZero (b[k]) ) {
				cerr << endl 
				     << "ERROR:  Found inconsistency in row " << k << "." << endl
				     << "	 No solution is possible." << endl << endl;
				_A = _AA;
				return x;
			}

		if ( i < 0 ) { // Zero matrix!
			cerr << endl << "ERROR:  Zero matrix.  No unique solution." << endl
			     << "        No solution found." << endl << endl;
			_A = _AA;
			return x;
		}

		if ( j != _n-1 ) { // Check to see if unique solution
			cerr << endl << "ERROR:  No uniquie solution." << endl
			     << "        No solution found." << endl << endl;
			_A = _AA;
			return x;
		}

		Element value;
		_F.init (value, 0);
		ConstRowIterator iter;

		while ( i >= 0 ) {
			iter = _A[i].begin ();
			x[(iter->first)] = bb[(iter->first)];
			iter++;

			for ( ; iter != _A[i].end (); iter++ ) // Subtract left values.
				_F.subin (x[(_A[i].begin ()->first)],
					 _F.mul (value,
						iter->second,
						x[(_A[iter->first].begin ()->first)]));

			_F.divin (x[(_A[i].begin ()->first)], _A[i].begin ()->second);
			// Divide by pivot

			if ( i > 0 ) { // Find new pivot element.
				i--; // move to previous row

				k = _A[i].begin ()->first;

				if ( k < j-1 ) { // Check to see if unique solution (redundant)
					cerr << endl << "ERROR:  No uniquie solution." << endl
					     << "	   No solution found." << endl
					     << "	   (Redundant Check)" << endl << endl;
					_A = _AA;
					return x;
				}

				j = k;
			}
			else
				break;

		} // while (i > 0)

		// retore original matrix A and vector b.
		_A = _AA;
		return x;

	} // Vector& SparseMatrixAux<DenseVectorTag>::linsolve (Vector&)

	template <class Field, class Row, class Vector>
	bool SparseMatrixAux<Field, Row, Vector, VectorCategories::DenseVectorTag>
		::gauss (Vector& b)
	{
		size_t i, j, k, bs;

		bs = b.size ();

		// Check to see if dimensions are compatible for solving the system.
		// Empty vectors are ignored
		if ( (bs != 0) && (bs != _m) ) {
			cerr << endl
			     << "ERROR:  Matrix and vector indices are incompatible." << endl
			     << "        No solution is possible." << endl << endl;
			return false;
		}

		i = 0;  // current row
		j = 0;  // current column
		Element value, temp;

		_F.init (value, 0);
		_F.init (temp, 0);

		while ( (i < _m - 1) && (j < _n) ) {
#ifdef TRACE
			clog << "In row " << i << " and column " << j << endl;
#endif // TRACE
      
			// look for non zero element in column
			for (k = i; k < _m; k++) 
				if ( _A[k].begin ()->first == j ) break;

			if (k < _m) { // found pivot
#ifdef TRACE
				clog << "  found pivot in row " << k << " and column " 
				     << _A[k].begin ()->first << endl;
#endif // TRACE
	
				if (i != k) { // if current row has zero (i.e., not pivot), swap rows
#ifdef TRACE
					clog << "  pivot is not in current row.  swapping rows " 
					     << i << " and " << k << endl;
#endif // TRACE
	  
					swaprow (i,k);
	  
					if (bs != 0) { // swap vector elements if not empty
						temp = b[i];
						b[i] = b[k];
						b[k] = temp;
					}
				}

#ifdef TRACE
				clog << "  zeroing out rest of column starting at row " 
				     << i + 1 << endl;
#endif // TRACE

				for (k = i+1; k < _m; k++) {
#ifdef TRACE
					clog << "    zeroing out row " << k 
					     << ", which has its first element in column " 
					     <<  _A[k].begin ()->first << endl;
#endif // TRACE

					if ( _A[k].begin ()->first == j ) {
						_F.negin (_F.div (value,
								  _A[k].begin ()->second,
								  _A[i].begin ()->second));

#ifdef TRACE
						clog << "      adding ";
						_F.write (clog, value);
						clog << " times row " << i << " to row " << k << endl;
#endif // TRACE

						addrow (i, k, value);
						if (bs != 0) {
							_F.mul (temp, value, b[i]);
							_F.addin (b[k], temp);
						}

#ifdef TRACE
						clog << "      Matrix is now:" << endl;
						write (clog);
						clog << "      and vector is now:" << endl;
//	    vector_utility<Field, Vector> VU (_F);
//	    VU.write (clog, b);
#endif // TRACE

					} // if ( _A[k].begin ()->first == j )
				} // for (k = i+1; k < _m; k++)
	
				i++;
				j++;
			} // found pivot
			else  // if there is no pivot, go to next column
				j++;
		} // while

		return true;
	} //  bool SparseMatrixAux<DenseVectorTag>::gauss (Vector&)

	template <class Field, class Row, class Vector>
	inline Vector &SparseMatrixAux<Field, Row, Vector, VectorCategories::DenseVectorTag>
		::apply (Vector& y, const Vector& x) const
	{
		if (_n != x.size ()) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
			return *(new Vector);
		}
 
		ConstRowIterator iter;
		typename Vector::iterator y_iter;

		Element temp;

		_F.init (temp, 0);

		y_iter = y.begin ();
		for (size_t i = 0; i < _m; i++, y_iter++)
			for (iter = _A[i].begin (); iter != _A[i].end (); iter++)
				_F.addin (*y_iter, _F.mul (temp, (*iter).second, x[(*iter).first]));

		return y;
 
	} // Vector& SparseMatrixAux<DenseVectorTag>::apply (Vector& y, const Vector&) const

	template <class Field, class Row, class Vector>
	inline Vector &SparseMatrixAux<Field, Row, Vector, VectorCategories::DenseVectorTag>
		::applyTranspose (Vector& y, const Vector& x) const
	{
		if (_m != x.size ()) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
			return *(new Vector);
		}
 
		ConstRowIterator iter;
		typename Vector::iterator y_iter;

		Element temp;

		_F.init (temp, 0);
 
		for (size_t i = 0; i < _m; i++, y_iter++)
			for (iter = _A[i].begin (); iter != _A[i].end (); iter++)
				_F.addin (y[(*iter).first], _F.mul (temp, (*iter).second, x[i]));
    
		return y;
 
	} // Vector& SparseMatrixAux<DenseVectorTag>::applyTranspose (Vector& y, const Vector&) const

	// Implementation of matrix methods for sparse sequence vectors

	template <class Field, class Row, class Vector>
	Vector& SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseSequenceVectorTag>
		::linsolve (Vector &x, const Vector& b)
	{ 
		// exceptions: if the system does not have a unique solution,
		//		 an error is printed and a zero matrix returned
		size_t bs = b.size ();
    
		// Check to see if dimensions are compatible for solving the system.
		if ( (bs != 0) && (_m < b.back ().first) )
		{
			cerr << endl
			     << "ERROR:  Matrix and vector indices are incompatible." << endl
			     << "        No solution is possible." << endl << endl;
			return x;
		} // if ( (bs != 0) && (_m < b.back ().first) )

		// make copies so can restore originals
		std::vector< Row > _AA (_A);
		Vector bb (b);

		// perform gaussian elimination
		if (!gauss (bb))
		{
			cerr << endl
			     << "ERROR:  Gaussian elimination not performed." << endl << endl;
			return x;
		} // if (!gauss (bb))

		// Back substitution to construct solution

		size_t i, j, k;

		// find first pivot element from bottom by scanning rows
		for (i = _m - 1; i >= 0; i--)
			if ( (_A[i]).begin () != (_A[i]).end () )
			{
				j = _A[i].begin ()->first;
				break;
			} // if ( (_A[i]).begin () != (_A[i]).end () )

		// Check for inconsistancy
		for (k = i+1; k < _m; k++)
			if ( (bb.back ().first == k) && (!_F.isZero (bb.back ().second)) )
			{
				cerr << endl 
				     << "ERROR:  Found inconsistency in row " << k << "." << endl
				     << "	 No solution is possible." << endl << endl;
				_A = _AA;
				return x;
			}

		if ( i < 0 ) // Zero matrix!
		{
			cerr << endl << "ERROR:  Zero matrix.  No unique solution." << endl
			     << "        No solution found." << endl << endl;
			_A = _AA;
			return x;
		}

		if ( j != _n-1 ) // Check to see if unique solution
		{
			cerr << endl << "ERROR:  No uniquie solution." << endl
			     << "        No solution found." << endl << endl;
			_A = _AA;
			return x;
		}
   
		Element value, temp;
		RowIterator iter;
		typename Vector::iterator b_iter, x_iter;

		_F.init (value, 0);
		_F.init (temp, 0);

		while ( i >= 0 )
		{
			// Set iterators.
			iter = _A[i].begin ();
			iter++;
			x_iter = x.begin ();

			// Obtain original value to start subtractions
			b_iter = lower_bound (bb.begin (), bb.end (), i, comp_w_index ());
			if ( (b_iter != bb.end ()) && (i == b_iter->first) )
				value = b_iter->second;
			else
				_F.init (value, 0);

			for ( ; iter != _A[i].end (); iter++ ) { // Subtract left values.
				k = _A[iter->first].begin ()->first;
      
				x_iter = lower_bound (x_iter, x.end (), k, comp_w_index ());

				if ( (x_iter != x.end ()) && (k == x_iter->first) )
					_F.subin (value, _F.mul (temp, iter->second, x_iter->second));
			}

			_F.divin (value, _A[i].begin ()->second); // Divide by pivot

			// Insert non-zero element in solution vector
			if (!_F.isZero (value))
				x.insert (x.begin (), make_pair (i, value));

			if ( i > 0 ) { // Find new pivot element.
				i--; // move to previous row
  	
				k = _A[i].begin ()->first;

				if ( k < j-1 ) { // Check to see if unique solution (redundant)
					cerr << endl << "ERROR:  No uniquie solution." << endl
					     << "	   No solution found." << endl
					     << "	   (Redundant Check)" << endl << endl;
					_A = _AA;
					return x;
				}

				j = k;

			}
			else
				break;

		} // while (i >= 0)

		// retore original matrix A and vector b.
		_A = _AA;

		return x;

	} // Vector& SparseMatrixAux<SparseSequenceVectorTag>::linsolve (Vector&)

	template <class Field, class Row, class Vector>
	bool SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseSequenceVectorTag>
		::gauss (Vector& b)
	{
		size_t i, j, k, bs;

		bs = b.size ();

		// Check to see if dimensions are compatible for solving the system.
		// Empty vectors are ignored
		if ( (bs != 0) && (_m < b.back ().first) ) {
			cerr << endl
			     << "ERROR:  Matrix and vector indices are incompatible." << endl
			     << "        No solution is possible." << endl << endl;
			return false;
		}

		i = 0;  // current row
		j = 0;  // current column
		Element value, temp;
		typename Vector::iterator iter_i, iter_k;

		_F.init (value, 0);
		_F.init (temp, 0);

		while ( (i < _m - 1) && (j < _n) ) {
			// look for non zero element in column
			for (k = i; k < _m; k++) 
				if ( _A[k].begin ()->first == j ) break;

			if (k < _m) {  // found pivot
				if (i != k) { // if current row has zero, swap rows
					swaprow (i,k);
	  
					if (bs != 0) { // swap vector elements if not empty
						// find where elements occur in vector
						// element k will always be after element i.
						iter_i = lower_bound (b.begin (), b.end (), i, comp_w_index ());
						iter_k = lower_bound (iter_i, b.end (), k, comp_w_index ());

						// perform swap
						if ( (iter_i != b.end ()) && (iter_k != b.end ()) 
						     && (i == iter_i->first) && (k == iter_k->first) )
						{
							temp = iter_i->second;
							iter_i->second = iter_k->second;
							iter_k->second = temp;
						}
						else if ( ( (iter_i != b.end ()) && (i == iter_i->first) )
							  && ( (iter_k == b.end ()) || (k != iter_k->first) ) )
						{
							b.insert (iter_k, make_pair (k, iter_i->second));
							iter_i = lower_bound (b.begin (), b.end (), i, comp_w_index ());
							b.erase (iter_i);
						}
						else if ( ( (iter_k != b.end ()) && (k == iter_k->first) )
							  && ( (iter_i == b.end ()) || (i != iter_i->first) ) )
						{
							iter_i = b.insert (iter_i, make_pair (i, iter_k->second));
							iter_k = lower_bound (iter_i, b.end (), k, comp_w_index ());
							b.erase (iter_k);
						}
					} // if (bs != 0) // swap vector elements if not empty
				} // if (i != k)  // if current row has zero, swap rows

				for (k = i+1; k < _m; k++)  // zero out rest of column
					if ( _A[k].begin ()->first == j ) {
						_F.negin (_F.div (value,
								_A[k].begin ()->second,
								_A[i].begin ()->second));
						addrow (i, k, value);
						if (bs != 0) {
							// find where elements occur in vector
							// element k will always be after element i.
							iter_i = lower_bound (b.begin (), b.end (), i, comp_w_index ());
							iter_k = lower_bound (iter_i, b.end (), k, comp_w_index ());
	      
							// perform swap
							if ( (iter_i != b.end ()) && (i == iter_i->first) )
							{
								if ( (iter_k != b.end ()) && (k == iter_k->first) )
								{
									_F.addin (iter_k->second, _F.mul (temp, value, iter_i->second));
								} // if (k == iter_k->first)
								else
									b.insert (iter_k, 
										 make_pair (k, 
											   _F.mul (temp, value, iter_i->second)));
							} // if (i == iter_i->first)
						}
					}

				i++;
				j++;
			} // if (k < _m)  // found pivot
			else  // if there is no pivot, go to next column
				j++;
		} // while ( (i < _m - 1) && (j < _n) )

		return true;
	} //  bool SparseMatrixAux<SparseSequenceVectorTag>::gauss (Vector&)

	template <class Field, class Row, class Vector>
	inline Vector &SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseSequenceVectorTag>
		::apply (Vector& y, const Vector& x) const
	{
		if ( (x.size () != 0) && (_n < x.back ().first) ) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
			return y;
		}
 
		ConstRowIterator iter;
		typename Vector::const_iterator x_iter;

		size_t k;

		Element zero;
		_F.init (zero, 0);
		Element value (zero), temp (zero);
 
		for (size_t i = 0; i < _m; i++) {
			value = zero;
			x_iter = x.begin ();

			for (iter = _A[i].begin (); iter != _A[i].end (); iter++) {
				k = (*iter).first;
      
				x_iter = lower_bound (x_iter, x.end (), k, comp_w_index ());

				if ( (x_iter != x.end ()) && (k == (*x_iter).first) )
					_F.addin (value, _F.mul (temp, (*iter).second, (*x_iter).second));
			}

			// Insert non-zero element in solution vector
			if (!_F.isZero (value))
				y.push_back (make_pair (i, value));

		}
    
		return y;
	} // Vector& SparseMatrixAux<SparseSequenceVectorTag>::apply (Vector&, const Vector&) const

	template <class Field, class Row, class Vector>
	inline Vector &SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseSequenceVectorTag>
		::applyTranspose (Vector& y, const Vector& x) const
	{
		if ( (x.size () != 0) && (_m < x.back ().first) ) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
			return y;
		}
 
		ConstRowIterator iter;
		typename Vector::const_iterator x_iter;

		size_t k;

		Element zero;
		_F.init (zero, 0);
		Element value (zero), temp (zero);
 
		std::vector<Element> _y (_n, zero); // temporary vector for calculating output
    
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++) {
			k = (*x_iter).first;       // vector index
			value = (*x_iter).second;  // value in vector

			// apply vector to column k
			for (iter = _A[k].begin (); iter != _A[k].end (); iter++)
				_F.addin (_y[(*iter).first], _F.mul (temp, (*iter).second, value));
             
		}

		// Convert temporary vector to sparse vector for output
		for (size_t i = 0; i < _y.size (); i++)
			if (!_F.isZero (_y[i])) y.push_back (make_pair (i, _y[i]));
 
		return y;
	} // Vector& SparseMatrixAux<SparseSequenceVectorTag>::applyTranspose (...) const

	// Implementation of matrix methods for sparse associative vectors

	template <class Field, class Row, class Vector>
	Vector &SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseAssociativeVectorTag>
		::linsolve (Vector &x, const Vector &b)
	{ 
		// exceptions: if the system does not have a unique solution,
		//		 an error is printed and a zero matrix returned
		size_t bs = b.size ();
    
		// Check to see if dimensions are compatible for solving the system.
		if ( (bs != 0) && (_m < b.rbegin ()->first) ) {
			cerr << endl
			     << "ERROR:  Matrix and vector indices are incompatible." << endl
			     << "        No solution is possible." << endl << endl;
			return x;
		}

		// make copies so can restore originals
		std::vector< Row > _AA (_A);
		Vector bb (b);

		// perform gaussian elimination
		if (!gauss (bb)) {
			cerr << endl
			     << "ERROR:  Gaussian elimination not performed." << endl << endl;
			return x;
		}

		// Back substitution to construct solution

		size_t i, j, k;

		// find first pivot element from bottom by scanning rows
		for (i = _m - 1; i >= 0; i--)
			if ( (_A[i]).begin () != (_A[i]).end () ) {
				j = _A[i].begin ()->first;
				break;
			}

		// Check for inconsistancy
		for (k = i+1; k < _m; k++)
			if ( (bb.rbegin ()->first == k) && (!_F.isZero (bb.rbegin ()->second)) ) {
				cerr << endl 
				     << "ERROR:  Found inconsistency in row " << k << "." << endl
				     << "	 No solution is possible." << endl << endl;
				_A = _AA;
				return x;
			}

		if ( i < 0 ) { // Zero matrix!
			cerr << endl << "ERROR:  Zero matrix.  No unique solution." << endl
			     << "        No solution found." << endl << endl;
			_A = _AA;
			return x;
		}

		if ( j != _n-1 ) { // Check to see if unique solution
			cerr << endl << "ERROR:  No uniquie solution." << endl
			     << "        No solution found." << endl << endl;
			_A = _AA;
			return x;
		}

		Element value, temp;
		RowIterator iter;
		typename Vector::iterator b_iter, x_iter;

		_F.init (value, 0);
		_F.init (temp, 0);

		while ( i >= 0 ) {
			// Set iterators.
			iter = _A[i].begin ();
			iter++;
			x_iter = x.begin ();

			// Obtain original value to start subtractions
			b_iter = bb.find (i);
			if (b_iter != bb.end ())
				value = b_iter->second;
			else
				_F.init (value, 0);

			for ( ; iter != _A[i].end (); iter++ ) { // Subtract left values.
				k = _A[iter->first].begin ()->first;
      
				x_iter = x.find (k);

				if (x_iter != x.end ())
					_F.subin (value, _F.mul (temp, iter->second, x_iter->second));
			}

			_F.divin (value, _A[i].begin ()->second); // Divide by pivot

			// Insert non-zero element in solution vector
			if (!_F.isZero (value))
				x.insert (x.begin (), make_pair (i, value));

			if ( i > 0 ) { // Find new pivot element.
				i--; // move to previous row
  	
				k = _A[i].begin ()->first;

				if ( k < j-1 ) { // Check to see if unique solution (redundant)
					cerr << endl << "ERROR:  No uniquie solution." << endl
					     << "	   No solution found." << endl
					     << "	   (Redundant Check)" << endl << endl;
					_A = _AA;
					return x;
				}

				j = k;
			}
			else
				break;

		} // while (i >= 0)

		// retore original matrix A and vector b.
		_A = _AA;

		return x;
	} // Vector& SparseMatrixAux<SparseAssociativeVectorTag>::linsolve (Vector&)

	template <class Field, class Row, class Vector>
	bool SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseAssociativeVectorTag>::gauss (Vector& b)
	{
		size_t i, j, k, bs;

		bs = b.size ();

		// Check to see if dimensions are compatible for solving the system.
		// Empty vectors are ignored
		if ( (bs != 0) && (_m < b.rbegin ()->first) ) {
			cerr << endl
			     << "ERROR:  Matrix and vector indices are incompatible." << endl
			     << "        No solution is possible." << endl << endl;
			return false;
		}

		i = 0;  // current row
		j = 0;  // current column
		Element value, temp;
		typename Vector::iterator iter_i, iter_k;

		_F.init (value, 0);
		_F.init (temp, 0);

		while ( (i < _m - 1) && (j < _n) )
		{
			// look for non zero element in column
			for (k = i; k < _m; k++) 
				if ( _A[k].begin ()->first == j ) break;

			if (k < _m) { // found pivot
				if (i != k) { // if current row has zero, swap rows
					swaprow (i,k);
	  
					if (bs != 0) { // swap vector elements if not empty
						// find where elements occur in vector
						// element k will always be after element i.
						iter_i = b.find (i);
						iter_k = b.find (k);

						// perform swap
						if ( (iter_i != b.end ()) && (iter_k != b.end ()) ) {
							temp = iter_i->second;
							iter_i->second = iter_k->second;
							iter_k->second = temp;
						}
						else if ( (iter_i != b.end ()) && (iter_k == b.end ()) ) {
							b.insert (make_pair (k, iter_i->second));
							iter_i = b.find (i);
							b.erase (iter_i);
						}
						else if ( (iter_k != b.end ()) && (iter_i == b.end ()) ) {
							b.insert (make_pair (i, iter_k->second));
							iter_k = b.find (k);
							b.erase (iter_k);
						}
					}
				}

				for (k = i+1; k < _m; k++)  // zero out rest of column
					if ( _A[k].begin ()->first == j )
					{
						_F.negin (_F.div (value,
								_A[k].begin ()->second,
								_A[i].begin ()->second));
						addrow (i, k, value);
						if (bs != 0)
						{
							// find where elements occur in vector
							// element k will always be after element i.
							iter_i = b.find (i);
							iter_k = b.find (k);
	      
							// perform addition
							if (iter_i != b.end ())
							{
								if (iter_k != b.end ())
								{
									_F.addin (iter_k->second, _F.mul (temp, value, iter_i->second));
								} // if (iter_k != b.end ())
								else
									b.insert (make_pair (k, 
											   _F.mul (temp, value, iter_i->second)));
							}
						}
					}

				i++;
				j++;
			} // if (k < _m)  // found pivot
			else  // if there is no pivot, go to next column
				j++;
		} // while ( (i < _m - 1) && (j < _n) )

		return true;
	} //  bool SparseMatrixAux<SparseAssociativeVectorTag>::gauss (Vector&)

	template <class Field, class Row, class Vector>
	inline Vector &SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseAssociativeVectorTag>
		::apply (Vector& y, const Vector& x) const
	{
		if ( (x.size () != 0) && (_n < x.rbegin ()->first) ) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
			return y;
		}
 
		ConstRowIterator iter;
		typename Vector::const_iterator x_iter;

		size_t k;

		Element value, temp;

		_F.init (value, 0);
		_F.init (temp, 0);
 
		for (size_t i = 0; i < _m; i++) {
			_F.init (value, 0);
			x_iter = x.begin ();

			for (iter = _A[i].begin (); iter != _A[i].end (); iter++) {
				k = (*iter).first;
      
				x_iter = x.find (k);

				if (x_iter != x.end ())
					_F.addin (value, _F.mul (temp, (*iter).second, (*x_iter).second));
			}

			// Insert non-zero element in solution vector
			if (!_F.isZero (value))
				y.insert (make_pair (i, value));

		} // for (size_t i = 0; i < _m; i++)
    
		return y;
	} // Vector& SparseMatrixAux<SparseAssociativeVectorTag>::apply (...) const

	template <class Field, class Row, class Vector>
	inline Vector &SparseMatrixAux<Field, Row, Vector, VectorCategories::SparseAssociativeVectorTag>
		::applyTranspose (Vector& y, const Vector& x) const
	{
		if ( (x.size () != 0) && (_m < x.rbegin ()->first) ) {
			cerr << endl << "ERROR:  Input vector not of right size." << endl << endl;
			return y;//*(new Vector);
		}
 
		ConstRowIterator iter;
		typename Vector::const_iterator x_iter;

		size_t k;

		Element zero;
		_F.init (zero, 0);
		Element value (zero), temp (zero);
 
		std::vector<Element> _y (_n, zero); // temporary vector for calculating output
    
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++) {
			k = (*x_iter).first;       // vector index
			value = (*x_iter).second;  // value in vector

			// apply vector to column k
			for (iter = _A[k].begin (); iter != _A[k].end (); iter++)
				_F.addin (_y[(*iter).first], _F.mul (temp, (*iter).second, value));
             
		}

		// Convert temporary vector to sparse vector for output
		for (size_t i = 0; i < y.size (); i++)
			if (!_F.isZero (_y[i])) y.insert (make_pair (i, _y[i]));
 
		return y;
	} // Vector& SparseMatrixAux<SparseAssociativeVectorTag>::applyTranspose (...) const

} // namespace LinBox

#endif // __SPARSE_MATRIX_AUX_H
