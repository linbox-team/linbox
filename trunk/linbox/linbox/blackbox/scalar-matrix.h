/* -*- mode: c; style: linux -*- */

/* linbox/blackbox/scalar.h
 * Copyright (C) 2002 by -bds  
 * evolved from diagonal.h written by William J Turner and Bradford Hovinen 
 *
 * -------------------------------
 * Modified by Dmitriy Morozov <linbox@foxcub.org>. May 28, 2002.
 *
 * Added parametrization of VectorCategory tags by VectorTraits. See
 * vector-traits.h for more details.
 *
 * -------------------------------
 */

#ifndef __SCALAR_H
#define __SCALAR_H

#include "linbox/blackbox/archetype.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"

namespace LinBox
{

	/** Blackbox Scalar Matrix.
	 * This is a class of blackbox square scalar matrices.
	 * Each scalar matrix occupies O(scalar-size) memory.
	 * This ScalarMatrix is a subclass of the
	 * {@link Archetypes archetype} for \Ref{BlackBox Matrices}.
	 * The matrix itself is not stored in memory, just the scalar and the dimensions.
	 * 
	 * This class has two template parameters.  The first is the field in 
	 * which the arithmetic is to be done.  The second is the type of 
	 * \Ref{LinBox} vector to which the matrix may be applied and of the vector result.
	 * 
	 * @param Field \Ref{LinBox} field or ring of the entries.
	 * @param Vector \Ref{LinBox} dense or sparse vector of field elements
	 */
	template <class Field, class Vector>//, class Trait = VectorTraits<Vector>::VectorCategory>
	class ScalarMatrix : public BlackboxArchetype<Vector>
	{
	    public:

	        typedef typename Field::Element        Element;

		/*  In each specialization, I must define suitable constructor(s) and
		BlackboxArchetype<Vector> * clone() const;
		Vector& apply(Vector& y, const Vector& x) const;
		Vector& applyTranspose(Vector& y, const Vector& x) const;
		size_t rowdim(void) const;
		size_t coldim(void) const;
		*/

		/** Constructor from field and dense vector of field elements.
		 * @param F	field in which to do arithmetic.
		 * @param n	size of the matrix.
		 * @param s	scalar, a field element, to be used as the diagonal of the matrix.
		 */
		ScalarMatrix (const Field &F, const size_t n, const Element &s)
			: _F(F), _n(n), _v(s) {}

		/** Constructor from field and random iterator over the field
		 * @param F    field in which to do arithmetic.
		 * @param n    size of the matrix.
		 * @param iter Random iterator from which to get the diagonal scalar element.
		 */
		ScalarMatrix (const Field &F, const size_t n, const typename Field::RandIter& iter)
			: _F(F), _n(n), _v(*iter) {}

		BlackboxArchetype<Vector> *clone() const
			{ return new ScalarMatrix(*this); }

		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires time linear in n, the size of the matrix.
		 */
                
		Vector& apply(Vector &y, const Vector &x) const 
			{   VectorTraits<Vector>::VectorCategory t;
			    return _app(y, x, t); }

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires time linear in n, the size of the matrix.
		 */
		Vector& applyTranspose(Vector &y, const Vector &x) const
			{ return apply(y, x); }  // symmetric matrix.

		size_t rowdim(void) const { return _n; }
    
		size_t coldim(void) const { return _n; }

	    protected:

		Field _F;   // Field for arithmetic

		size_t _n;  // Number of rows and columns of square matrix.

		Element _v; // the scalar used in applying matrix.

		// dense vector _app for apply
		template <class VectorTrait> Vector& _app (Vector &y, const Vector &x, VectorCategories::DenseVectorTag<VectorTrait>) const;
		// The third argument is just a device to let overloading determine the method.

		// sparse sequence vector _app for apply
		template <class VectorTrait> Vector& _app (Vector &y, const Vector &x, VectorCategories::SparseSequenceVectorTag<VectorTrait>) const;

		// sparse associative vector _app for apply
		template <class VectorTrait> Vector& _app (Vector &y, const Vector &x, VectorCategories::SparseAssociativeVectorTag<VectorTrait>) const;

	}; // template <Field, Vector> class ScalarMatrix
   
	// dense vector _app
	template <class Field, class Vector>
		template <class VectorTrait>
	inline Vector &ScalarMatrix<Field, Vector>
		::_app(Vector& y, const Vector& x, VectorCategories::DenseVectorTag<VectorTrait> t) const
		{   
		    linbox_check (x.size() >= _n);
		    typename Vector::iterator y_iter = y.begin ();
		    typename Vector::const_iterator y_end = y.begin () + _n;

		    if (_F.isZero(_v)) // just write zeroes
		        for ( ; y_iter != y.end ();  ++y_iter) *y_iter = _v;
                    else if (_F.isOne(_v) ) // just copy 
			copy(x.begin(), x.end(), y.begin());
		    else // use actual muls
		    {   typename Vector::const_iterator x_iter = x.begin ();
		            for (  ; y_iter != y.end () ; ++y_iter, ++x_iter )
		                _F.mul (*y_iter, _v, *x_iter);
		    }
		    return y;

		} // dense vector _app

		
	// sparse sequence vector _app
	template <class Field, class Vector>
		template <class VectorTrait>
	inline Vector &ScalarMatrix<Field, Vector>
		::_app(Vector& y, const Vector& x, VectorCategories::SparseSequenceVectorTag<VectorTrait> t) const
	{
		//linbox_check ((!x.empty ()) && (_n < x.back ().first));
		// neither is required of x ?

		y.clear (); // we'll overwrite using push_backs.

		// field element to be used in calculations
		Element entry;
		_F.init (entry, 0); // needed?

		// For each element, multiply input element with corresponding element
		// of stored scalar and insert non-zero elements into output vector
		for ( typename Vector::const_iterator x_iter = x.begin (); x_iter != x.end (); ++x_iter) 
		{	_F.mul (entry, _v, x_iter->second);
			if (!_F.isZero (entry)) y.push_back (make_pair (x_iter->first, entry));
		} 

		return y;
	} // sparse sequence vector _app

	// sparse associative vector _app
	template <class Field, class Vector>
		template <class VectorTrait>
	inline Vector& ScalarMatrix<Field, Vector>
		::_app(Vector& y, const Vector& x, VectorCategories::SparseAssociativeVectorTag<VectorTrait> t) const
	{
		y.clear (); // we'll overwrite using inserts

		// create field elements and size_t to be used in calculations
		Element entry;
		_F.init (entry, 0);

		// Iterator over indices of input vector.
		// For each element, multiply input element with 
		// stored scalar and insert non-zero elements into output vector
		for ( typename Vector::const_iterator x_iter = x.begin (); x_iter != x.end (); ++x_iter)
		{	_F.mul (entry, _v, x_iter->second);
			if (!_F.isZero (entry)) y.insert (y.end (), make_pair (xiter->first, entry));
		}

		return y;
	} // sparse associative vector _app
} // namespace LinBox

#endif // __ScalarMatrix
