/* -*- mode: c; style: linux -*- */

/* linbox/src/blackbox/diagonal.h
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

#ifndef __DIAGONAL_H
#define __DIAGONAL_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/debug.h"

// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** Blackbox Diagonal matrix.
	 * This is a class of n by n diagonal matrices templatized by the 
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
	class Diagonal : public Blackbox_archetype<Vector>
	{
	    public:

	        typedef Blackbox_archetype<Vector>     Blackbox;
	        typedef typename Field::Element        Element;
		typedef typename Field::RandomIterator RandomIterator;

		/** Constructor from field and dense vector of field elements.
		 * @param F	LinBox field in which to do arithmetic
		 * @param v LinBox dense vector of field elements to be used 
		 * 		as the diagonal of the matrix.
		 */
		Diagonal (const Field F, const std::vector<typename Field::Element>& v);

		/** Constructor from field and random iterator over the field
		 * @param F    LinBox field in which to do arithmetic
		 * @param iter Random iterator from which to get the diagonal elements
		 */
		Diagonal (const Field F, const RandomIterator &iter);

		/** Virtual constructor.
		 * Required because constructors cannot be virtual.
		 * Make a copy of the Blackbox_archetype object.
		 * Required by abstract base class.
		 * @return pointer to new blackbox object
		 */
		Blackbox *clone() const;

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& apply(Vector& y, const Vector& x) const;

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires one vector conforming to the \Ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * Because the diagonal matrix is symmetric, this is the same as calling 
		 * the apply function.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		Vector& applyTranspose(Vector& y, const Vector& x) const;

		/** Retreive row dimensions of BlackBox matrix.
		 * This may be needed for applying preconditioners.
		 * Required by abstract base class.
		 * @return integer number of rows of black box matrix.
		 */
		size_t rowdim(void) const;
    
		/** Retreive column dimensions of BlackBox matrix.
		 * Required by abstract base class.
		 * @return integer number of columns of black box matrix.
		 */
		size_t coldim(void) const;

	}; // template <Field, Vector> class Diagonal
 
	// Specialization of diagonal for LinBox dense vectors
	template <class Field, class Vector>
	class Diagonal<Field, Vector, VectorCategories::DenseVectorTag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;
		typedef typename Field::Element    Element;

		Diagonal(const Field F, const std::vector<typename Field::Element>& v);
		Blackbox *clone() const 
			{ return new Diagonal(*this); }
		Vector& apply(Vector& y, const Vector& x) const;
		Vector& applyTranspose(Vector& y, const Vector& x) const { return apply(y, x); }
		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<DenseVectorTag>
   
	// Specialization of diagonal for LinBox sparse sequence vectors
	template <class Field, class Vector>
	class Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;
		typedef typename Field::Element    Element;

		Diagonal(const Field F, const std::vector<typename Field::Element>& v);
		Blackbox *clone() const 
			{ return new Diagonal(*this); }
		Vector& apply(Vector& y, const Vector& x) const;
		Vector& applyTranspose(Vector& y, const Vector& x) const { return apply(y, x); }
		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<SparseSequenceVectorTag>

	// Specialization of diagonal for LinBox sparse associative vectors
	template <class Field, class Vector>
	class Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag>
		: public Blackbox_archetype<Vector>
	{
	    public:

		typedef Blackbox_archetype<Vector> Blackbox;
		typedef typename Field::Element    Element;

		Diagonal(const Field F, const std::vector<typename Field::Element>& v);
		Blackbox *clone() const 
			{ return new Diagonal(*this); }
		Vector& apply(Vector& y, const Vector& x) const;
		Vector& applyTranspose(Vector& y, const Vector& x) const { return apply(y, x); }
		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<SparseAssociativeVectorTag>

	// Method implementations for dense vectors
 
	template <class Field, class Vector>
	inline Diagonal<Field, Vector, VectorCategories::DenseVectorTag>
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}

	template <class Field, class Vector>
	inline Vector& Diagonal<Field, Vector, VectorCategories::DenseVectorTag>
		::apply(Vector& y, const Vector& x) const
	{
		// Resize y to match the correct dimension
		y.resize (_n);
 
		linbox_check (_n == x.size ());
 
		// Create iterators for input, output, and stored vectors
		std::vector<Element>::const_iterator v_iter;
		typename Vector::const_iterator x_iter;
		typename Vector::iterator y_iter;
 
		// Start at beginning of _v and x vectors
		v_iter = _v.begin ();
		x_iter = x.begin ();

		// Iterate through all three vectors, multiplying input and stored
		// vector elements to create output vector element.
		for (y_iter = y.begin ();
		     y_iter != y.end ();
		     y_iter++, v_iter++, x_iter++)
			_F.mul (*y_iter, *v_iter, *x_iter);
 
		return y;
	} // Vector& Diagonal<DenseVectorTag>::apply(Vector& y, const Vector&) const
  
	// Method implementations for sparse sequence vectors
 
	template <class Field, class Vector>
	inline Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag>
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}

	template <class Field, class Vector>
	inline Vector &Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag>
		::apply(Vector& y, const Vector& x) const
	{
		linbox_check ((!x.empty ()) && (_n < x.back ().first));

		y.clear (); // we'll overwrite using push_backs.

		// create field elements and size_t to be used in calculations
		size_t i;
		Element zero, entry;
		_F.init (zero, 0);
		_F.init (entry, 0);

		// Create iterators for input and stored vectors
		std::vector<Element>::const_iterator v_iter;
		typename Vector::const_iterator x_iter;
 
		// Start at beginning of _v vector
		v_iter = _v.begin ();
 
		// Iterator over indices of input vector.
		// For each element, multiply input element with corresponding element
		// of stored vector and insert non-zero elements into output vector
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++) {
			i = (*x_iter).first;
			_F.mul (entry, *(v_iter + i), (*x_iter).second);
			if (!_F.isZero (entry)) y.push_back (make_pair (i, entry));
		} // for (x_iter = x.begin (); x_iter != x.end (); x_iter++)

		return y;
	} // Vector& Diagonal<SparseSequenceVectorTag>::apply(Vector& y, const Vector&) const

	// Method implementations for sparse associative vectors
 
	template <class Field, class Vector>
	inline Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag>
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}

	template <class Field, class Vector>
	inline Vector& Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag>
		::apply(Vector& y, const Vector& x) const
	{
		linbox_check ((!x.empty ()) && (_n < x.rbegin ()->first));

		y.clear (); // we'll overwrite using inserts

		// create field elements and size_t to be used in calculations
		size_t i;
		Element zero, entry;
		_F.init (zero, 0);
		_F.init (entry, 0);

		// Create iterators for input and stored vectors
		std::vector<Element>::const_iterator v_iter;
		typename Vector::const_iterator x_iter;
 
		// Start at beginning of _v vector
		v_iter = _v.begin ();
 
		// Iterator over indices of input vector.
		// For each element, multiply input element with corresponding element
		// of stored vector and insert non-zero elements into output vector
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++)
		{
			i = x_iter->first;
			_F.mul (entry, *(v_iter + i), (*x_iter).second);
			if (!_F.isZero (entry)) y.insert (y.end (), make_pair (i, entry));
		}

		return y;
	} // Vector& Diagonal<SparseAssociativeVectorTag>::apply(...) const

} // namespace LinBox

#endif // __DIAGONAL_H
