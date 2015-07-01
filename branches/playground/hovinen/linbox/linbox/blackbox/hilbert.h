/* -*- mode: c; style: linux -*- */

/* linbox/src/blackbox/hilbert.h
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

#ifndef __HILBERT_H
#define __HILBERT_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"

#include "linbox/debug.h"

namespace LinBox
{

	/** Blackbox Hilbert matrix.
	 * This is a class of n by n Hilbert matrices templatized by the 
	 * {@link Fields field} in 
	 * which the elements reside.  The class conforms to the 
	 * {@link Archetypes archetype} for \Ref{BlackBox Matrices}.
	 *
	 * The matrix itself is not stored in memory.  Rather, its apply
	 * methods use a vector of {@link Fields field} elements, which are 
	 * used to "multiply" the matrix to a vector.
	 * 
	 * This class has three template parameters.  The first is the field in 
	 * which the arithmetic is to be done.  The second is the type of 
	 * \Ref{LinBox} vector to which to apply the matrix.  The 
	 * third is chosen be defualt to be the \Ref{LinBox} vector trait
	 * of the vector.  This class is then specialized for dense and sparse 
	 * vectors.
	 * 
	 * The default class is not implemented.  It's functions should never
	 * be called because partial template specialization should always be
	 * done on the vector traits.
	 * @param Field \Ref{LinBox} field
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 * @param Trait  Marker whether to use dense or sparse LinBox vector 
	 *               implementation.  This is chosen by a default parameter 
	 *               and partial template specialization.
	 */
	template <class Field, class Vector, class Trait = VectorTraits<Vector>::VectorCategory>
	class Hilbert : public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;
		typedef typename Field::Element Element;

		/** Constructor from integer and field.
		 * @param n size_t integer number of rows and columns of matrix.
		 */
		Hilbert (Field F, size_t n);

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the Blackbox_archetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Blackbox *clone () const;

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& apply (Vector& y, const Vector& x) const;

		/** Application of BlackBox matrix transpose.
		 * y= transpose (A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * Because the Hilbert matrix is symmetric, this is the same as calling 
		 * the apply function.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& applyTranspose (Vector& y, const Vector& x) const;

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim (void) const;
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim (void) const;

	}; // template <Field, Vector> class hilbert
 
	// Specialization of hilbert for LinBox dense vectors
	template <class Field, class Vector>
	class Hilbert<Field, Vector, VectorCategories::DenseVectorTag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;
		typedef typename Field::Element    Element;

		Hilbert (Field F, size_t n);
		Blackbox *clone () const 
			{ return new Hilbert (*this); }
		Vector& apply (Vector& y, const Vector& x) const;
		Vector& applyTranspose (Vector& y, const Vector& x) const { return apply (y, x); }
		size_t rowdim (void) const { return _n; } 
		size_t coldim (void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _H;
    
	}; // template <Field, Vector> class hilbert<DenseVectorTag>
   
	// Specialization of hilbert for LinBox sparse sequence vectors
	template <class Field, class Vector>
	class Hilbert<Field, Vector, VectorCategories::SparseSequenceVectorTag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;
		typedef typename Field::Element    Element;

		Hilbert (Field F, size_t n);
		Blackbox* clone () const 
			{ return new Hilbert (*this); }
		Vector& apply (Vector& y, const Vector& x) const;
		Vector& applyTranspose (Vector& y, const Vector& x) const { return apply (y, x); }
		size_t rowdim (void) const { return _n; } 
		size_t coldim (void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _H;
    
	}; // template <Field, Vector> class hilbert<SparseSequenceVectorTag>

	// Specialization of hilbert for LinBox sparse associative vectors
	template <class Field, class Vector>
	class Hilbert<Field, Vector, VectorCategories::SparseAssociativeVectorTag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;
		typedef typename Field::Element    Element;

		Hilbert (Field F, size_t n);
		Blackbox_archetype<Vector>* clone () const 
			{ return new Hilbert (*this); }
		Vector& apply (Vector& y, const Vector& x) const;
		Vector& applyTranspose (Vector& y, const Vector& x) const { return apply (y, x); }
		size_t rowdim (void) const { return _n; } 
		size_t coldim (void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _H;
    
	}; // template <Field, Vector> class hilbert<sparse_associative_vector_tag>

	// Method implementations for dense vectors
 
	template <class Field, class Vector>
	inline Hilbert<Field, Vector, VectorCategories::DenseVectorTag>
		::Hilbert (Field F, size_t n) : _F (F), _n (n)
	{
		Element one = F.one (), temp = F.zero ();

		_H = std::vector<Element> (2*_n - 1, temp);

		std::vector<Element>::iterator iter;

		for (iter = _H.begin (); iter!= _H.end (); iter++) {
			_F.addin (temp, one);
			_F.div (*iter, one, temp);
		}
	}

	template <class Field, class Vector>
	inline Vector& Hilbert<Field, Vector, VectorCategories::DenseVectorTag>
		::apply (Vector& y, const Vector& x) const
	{
		// Create zero vector to hold output
		Element temp;

		_F.init (temp, 0);
 
		linbox_check (_n == x.size ());
 
		// Create iterators for input, output, and stored vectors
		std::vector<Element>::const_iterator iter, start_iter;
		typename Vector::const_iterator x_iter;
		typename Vector::iterator y_iter;
 
		// Start at beginning of _H vector for first row
		start_iter = _H.begin ();

		// Set y to the right size
		y.resize (_n);

		// Iterator over elements of output vector.
		// For each element, multiply row of matrix with input vector.
		// Each row of matrix starts one further in _H vector.
		for (y_iter = y.begin (); y_iter != y.end (); y_iter++, start_iter++) {
			// start matrix row at correct place
			iter = start_iter;
 
			// Multiply matrix row and input vector by iterating over both.
			for (x_iter = x.begin (); x_iter != x.end (); x_iter++, iter++) {
				_F.mul (temp, *iter, *x_iter);
				_F.addin (*y_iter, temp);
			}
		}
 
		return y;
	} // Vector& hilbert<DenseVectorTag>::apply (Vector&, const Vector&) const
  
	// Method implementations for sparse sequence vectors
 
	template <class Field, class Vector>
	inline Hilbert<Field, Vector, VectorCategories::SparseSequenceVectorTag>
		::Hilbert (Field F, size_t n) : _F (F), _n (n)
	{
		Element temp = F.zero ();

		_H = std::vector<Element> (2*_n - 1, temp);

		std::vector<Element>::iterator iter;

		for (iter = _H.begin (); iter!= _H.end (); iter++) {
			_F.addin (temp, _F.one ());
			_F.inv (*iter, temp);
		}
	}

	template <class Field, class Vector>
	inline Vector &Hilbert<Field, Vector, VectorCategories::SparseSequenceVectorTag>
		::apply (Vector& y, const Vector& x) const
	{
		linbox_check (x.empty () || _n >= x.back ().first);

		// create field elements to be used in calculations
		Element entry, temp;

		_F.init (entry, 0);
		_F.init (temp, 0);

		// Create iterators for input, output, and stored vectors
		std::vector<Element>::const_iterator iter, start_iter;
		typename Vector::const_iterator x_iter;
 
		// Start at beginning of _H vector for first row
		start_iter = _H.begin ();
 
		// Iterator over indices of output vector.
		// For each element, multiply row of matrix with input vector,
		// and insert non-zero elements into vector
		// Each row of matrix starts one further in _H vector.
		for (size_t i = 0; i < _n; i++, start_iter++) {
			_F.init (entry, 0);
 
			// Multiply matrix row and input vector by iterating over both.
			for (x_iter = x.begin (); x_iter != x.end (); x_iter++, iter++) {
				_F.mul (temp, *(start_iter + (*x_iter).first), (*x_iter).second);
				_F.addin (entry, temp);
			}

			if (!_F.isZero (entry)) y.push_back (make_pair (i, entry));
		}

		return y;
	} // Vector& hilbert<SparseSequenceVectorTag>::apply (Vector&, const Vector&) const

	// Method implementations for sparse associative vectors
 
	template <class Field, class Vector>
	inline Hilbert<Field, Vector, VectorCategories::SparseAssociativeVectorTag>
		::Hilbert (Field F, size_t n) : _F (F), _n (n)
	{
		Element temp = F.zero ();

		_H = std::vector<Element> (2*_n - 1, temp);

		std::vector<Element>::iterator iter;

		for (iter = _H.begin (); iter!= _H.end (); iter++) {
			_F.addin (temp, F.one ());
			_F.inv (*iter, temp);
		}
 
	} // hilbert<sparse_associative_vector_tag>::hilbert (Field, size_t)
	
	template <class Field, class Vector>
	inline Vector& Hilbert<Field, Vector, VectorCategories::SparseAssociativeVectorTag>
		::apply (Vector& y, const Vector& x) const
	{
		linbox_check (x.empty () || _n >= x.rbegin ()->first);

		// Resize y
		y.resize (_n);
 
		// create field elements to be used in calculations
		Element entry, temp;

		_F.init (entry, 0);
		_F.init (temp, 0);

		// Create iterators for input, output, and stored vectors
		std::vector<Element>::const_iterator iter, start_iter;
		typename Vector::const_iterator x_iter;
 
		// Start at beginning of _H vector for first row
		start_iter = _H.begin ();
 
		// Iterator over indices of output vector.
		// For each element, multiply row of matrix with input vector,
		// and insert non-zero elements into vector
		// Each row of matrix starts one further in _H vector.
		for (size_t i = 0; i < _n; i++, start_iter++) {
			entry = F.zero ();
 
			// Multiply matrix row and input vector by iterating over both.
			for (x_iter = x.begin (); x_iter != x.end (); x_iter++, iter++) {
				_F.mul (temp, *(start_iter + x_iter->first), x_iter->second); 
				_F.addin (entry, temp);
			}

			if (!_F.isZero (entry)) y.push_back (make_pair (i, entry));
		}

		return y;
	} // Vector& hilbert<SparseAssociativeVectorTag>::apply (...) const

} // namespace LinBox

#endif // _HILBERT_
