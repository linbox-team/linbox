/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/hilbert.h
 * Copyright (C) 1999-2001 William J Turner,
 *               2001 Bradford Hovinen
 *
 * Written by William J Turner <wjturner@math.ncsu.edu>,
 *            Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * ------------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 28, 2002.
 *
 * Added parametrization of VectorCategory tags by VectorTraits. See 
 * vector-traits.h for more details.
 *  
 * ------------------------------------
 * 
 * See COPYING for license information.
 */

#ifndef __HILBERT_H
#define __HILBERT_H

#include "linbox/vector/vector-traits.h"
#include <linbox/blackbox/blackbox-interface.h>

#include "linbox/util/debug.h"

namespace LinBox
{

	/** \brief Example of a blackbox that is space efficient, though not time efficient.

\ingroup blackbox
	 * This is a class of n by n Hilbert matrices templatized by the 
	 * {@link Fields field} in 
	 * which the elements reside.  The class conforms to the 
	 * {@link Archetypes archetype} for \ref{BlackBox Matrices}.
	 *
	 * The matrix itself is not stored in memory.  Rather, its apply
	 * methods use a vector of {@link Fields field} elements, which are 
	 * used to "multiply" the matrix to a vector.
	 * 
	 * This class has three template parameters.  The first is the field in 
	 * which the arithmetic is to be done.  The second is the type of 
	 * \ref{LinBox} vector to which to apply the matrix.  The 
	 * third is chosen be defualt to be the \ref{LinBox} vector trait
	 * of the vector.  This class is then specialized for dense and sparse 
	 * vectors.
	 * 
	 * The default class is not implemented.  It's functions should never
	 * be called because partial template specialization should always be
	 * done on the vector traits.
	 * @param Field \ref{LinBox} field
	 * @param Vector \ref{LinBox} dense or sparse vector of field elements
	 * @param Trait  Marker whether to use dense or sparse LinBox vector 
	 *               implementation.  This is chosen by a default parameter 
	 *               and partial template specialization.
	 */
	template <class _Field, class Trait = typename VectorTraits<typename Vector<_Field>::Dense>::VectorCategory>
	class Hilbert : public BlackboxInterface
	{
	    public:

		typedef _Field Field;
		typedef typename Field::Element Element;

		/** Constructor from integer and field.
		 * @param n size_t integer number of rows and columns of matrix.
		 */
		Hilbert (Field F, size_t n);


		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		OutVector& apply (OutVector& y, const InVector& x) const;

		/** Application of BlackBox matrix transpose.
		 * y= transpose (A)*x.
		 * Requires one vector conforming to the \ref{LinBox}
		 * vector {@link Archetypes archetype}.
		 * Required by abstract base class.
		 * Because the Hilbert matrix is symmetric, this is the same as calling 
		 * the apply function.
		 * @return reference to vector y containing output.
		 * @param  x constant reference to vector to contain input
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector& y, const InVector& x) const;

            template<typename _Tp1, typename _Trt1 = Trait>
            struct rebind
            { typedef Hilbert<_Tp1, _Trt1> other; };

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
	template <class _Field>
	class Hilbert<_Field,  VectorCategories::DenseVectorTag >
	{
	    public:

		typedef _Field Field;
		typedef typename Field::Element    Element;

		Hilbert (Field F, size_t n);

		template<class OutVector, class InVector>
		OutVector& apply (OutVector& y, const InVector& x) const;


		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector& y, const InVector& x) const { return apply (y, x); }
		size_t rowdim (void) const { return _n; } 
		size_t coldim (void) const { return _n; } 
		const Field& field() const { return _F; }
	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _H;
    
	}; // template <Field, Vector> class hilbert<DenseVectorTag>
   
	// Specialization of hilbert for LinBox sparse sequence vectors
	template <class _Field>
	class Hilbert<_Field, VectorCategories::SparseSequenceVectorTag >
	{
	    public:

		typedef _Field Field;
		typedef typename Field::Element    Element;

		Hilbert (Field F, size_t n);

		template<class OutVector, class InVector>
		OutVector& apply (OutVector& y, const InVector& x) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector& y, const InVector& x) const { return apply (y, x); }

		size_t rowdim (void) const { return _n; } 
		size_t coldim (void) const { return _n; } 

		const Field& field() const { return _F; }

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _H;
    
	}; // template <Field, Vector> class hilbert<SparseSequenceVectorTag>

	// Specialization of hilbert for LinBox sparse associative vectors
	template <class _Field>
	class Hilbert<_Field, VectorCategories::SparseAssociativeVectorTag >
	{
	    public:

		typedef _Field Field;
		typedef typename Field::Element    Element;

		Hilbert (Field F, size_t n);

		template<class OutVector, class InVector>
		OutVector& apply (OutVector& y, const InVector& x) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose (OutVector& y, const InVector& x) const { return apply (y, x); }

		size_t rowdim (void) const { return _n; } 
		size_t coldim (void) const { return _n; } 
		const Field& field() const { return _F; }

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _H;
    
	}; // template <Field, Vector> class hilbert<sparse_associative_vector_tag>

	// Method implementations for dense vectors
 
	template <class Field>
	inline Hilbert<Field, VectorCategories::DenseVectorTag >
		::Hilbert (Field F, size_t n) : _F (F), _n (n)
	{
		Element one, temp;

		_F.init (one, 1);
		_F.init (temp, 0);

		_H = std::vector<Element> (2*_n - 1, temp);

		typename std::vector<Element>::iterator iter;

		for (iter = _H.begin (); iter!= _H.end (); iter++) {
			_F.addin (temp, one);
			_F.div (*iter, one, temp);
		}
	}

	template <class Field>
	template <class OutVector, class InVector>
	inline OutVector& Hilbert<Field, VectorCategories::DenseVectorTag >
		::apply (OutVector& y, const InVector& x) const
	{
		// Create iterators for input, output, and stored vectors
		typename std::vector<Element>::const_iterator iter, start_iter;
		typename InVector::const_iterator x_iter;
		typename OutVector::iterator y_iter;
 

		// Iterator over elements of output vector.
		// For each element, multiply row of matrix with input vector.
		// Start at beginning of _H vector for first row
		// Each row of matrix starts one further in _H vector.
		for (y_iter = y.begin (), start_iter = _H.begin ();
		     y_iter != y.end (); 
		     y_iter++, start_iter++) 
		{
		        _F.init (*y_iter, 0);
			// Multiply matrix row and input vector by iterating over both.
			// start matrix row at correct place
			for (x_iter = x.begin (), iter = start_iter;
			     x_iter != x.end (); 
			     x_iter++, iter++) 
			{ _F.axpyin (*y_iter, *iter, *x_iter); }
		}
 
		return y;
	} // Vector& hilbert<DenseVectorTag>::apply (Vector&, const Vector&) const
  
	// Method implementations for sparse sequence vectors
 
	// Note: sparse vector code has not been fixed.  -bds 03Jan
	template <class Field>
	inline Hilbert<Field, VectorCategories::SparseSequenceVectorTag >
		::Hilbert (Field F, size_t n) : _F (F), _n (n)
	{
		Element temp = F.zero ();

		_H = std::vector<Element> (2*_n - 1, temp);

		typename std::vector<Element>::iterator iter;

		for (iter = _H.begin (); iter!= _H.end (); iter++) {
			_F.addin (temp, _F.one ());
			_F.inv (*iter, temp);
		}
	}

	template <class Field>
	template <class OutVector, class InVector>
	inline OutVector &Hilbert<Field, VectorCategories::SparseSequenceVectorTag >
		::apply (OutVector& y, const InVector& x) const
	{
		linbox_check (x.empty () || _n >= x.back ().first);

		// create field elements to be used in calculations
		Element entry, temp;

		_F.init (entry, 0);
		_F.init (temp, 0);

		// Create iterators for input, output, and stored vectors
		typename std::vector<Element>::const_iterator iter, start_iter;
		typename InVector::const_iterator x_iter;
 
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
 
	template <class Field>
	inline Hilbert<Field, VectorCategories::SparseAssociativeVectorTag >
		::Hilbert (Field F, size_t n) : _F (F), _n (n)
	{
		Element temp = F.zero ();

		_H = std::vector<Element> (2*_n - 1, temp);

		typename std::vector<Element>::iterator iter;

		for (iter = _H.begin (); iter!= _H.end (); iter++) {
			_F.addin (temp, F.one ());
			_F.inv (*iter, temp);
		}
 
	} // hilbert<sparse_associative_vector_tag>::hilbert (Field, size_t)
	
	template <class Field>
	template <class OutVector, class InVector>
	inline OutVector& Hilbert<Field, VectorCategories::SparseAssociativeVectorTag >
		::apply (OutVector& y, const InVector& x) const
	{
		linbox_check (x.empty () || _n >= x.rbegin ()->first);

		// Resize y
		y.resize (_n);
 
		// create field elements to be used in calculations
		Element entry, temp;

		_F.init (entry, 0);
		_F.init (temp, 0);

		// Create iterators for input, output, and stored vectors
		typename std::vector<Element>::const_iterator iter, start_iter;
		typename InVector::const_iterator x_iter;
 
		// Start at beginning of _H vector for first row
		start_iter = _H.begin ();
 
		// Iterator over indices of output vector.
		// For each element, multiply row of matrix with input vector,
		// and insert non-zero elements into vector
		// Each row of matrix starts one further in _H vector.
		for (size_t i = 0; i < _n; i++, start_iter++) {
			entry = _F.zero ();
 
			// Multiply matrix row and input vector by iterating over both.
			for (x_iter = x.begin (); x_iter != x.end (); x_iter++, iter++) {
				_F.mul (temp, *(start_iter + x_iter->first), x_iter->second); 
				_F.addin (entry, temp);
			}

			if (!_F.isZero (entry)) y.push_back (make_pair (i, entry));
		}

		return y;
	}// Vector& hilbert<SparseAssociativeVectorTag>::apply (...) const

} // namespace LinBox

#endif // _HILBERT_
