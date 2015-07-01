/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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

#include <algorithm>
#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox-config.h"

#ifdef __LINBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <string>
#include <iostream>

using std::string;
using std::istream;
using std::ostream;

#endif


namespace LinBox
{

	/** @memo Blackbox for aI. Useful particularly for representing 0 and I in example construction.
	 * @doc
	 * This is a class of blackbox square scalar matrices.
	 * Each scalar matrix occupies O(scalar-size) memory.
	 * This ScalarMatrix is a subclass of the
	 * {@link Archetypes archetype} for \Ref{BlackBox Matrices}.
	 * The matrix itself is not stored in memory, just the scalar and the dimensions.
	 * 
	 * This class has two template parameters.  The first is the \ref{Field} in 
	 * which the arithmetic is to be done.  The second is the type of 
	 * LinBox \Ref{Vector} to which the matrix may be applied and of the vector result.
	 */
	template <class _Field>
	class ScalarMatrix 
	{
	    public:
		
		typedef _Field Field;
	        typedef typename Field::Element        Element;

		/*  In each specialization, I must define suitable constructor(s) and
		BlackboxArchetype<Vector> * clone() const;
		Vector& apply(Vector& y, const Vector& x) const;
		Vector& applyTranspose(Vector& y, const Vector& x) const;
		size_t rowdim(void) const;
		size_t coldim(void) const;
		*/

		ScalarMatrix ()	:  _n(0) {}

		/** Scalar matrix Constructor from an element.
		 * @param F	field in which to do arithmetic.
		 * @param n	size of the matrix.
		 * @param s	scalar, a field element, to be used as the diagonal of the matrix.
		 */
		ScalarMatrix (const Field &F, const size_t n, const Element &s)
			: _F(F), _n(n), _v(s) {}

		/** Constructor from a random element.
		 * @param F    field in which to do arithmetic.
		 * @param n    size of the matrix.
		 * @param iter Random iterator from which to get the diagonal scalar element.
		 */
		ScalarMatrix (const Field &F, const size_t n, const typename Field::RandIter& iter)
			: _F(F), _n(n), _v(*iter) {}


#ifdef __LINBOX_XMLENABLED
		ScalarMatrix(Reader &R) : _F(R.Down(1))
		{
			
			size_t i, j;

			R.Up(1);
			if(!R.expectTagName("MatrixOver")) return;
			if(!R.expectAttributeNum("rows", i) || !R.expectAttributeNum("cols", j)) return;
			if(i != j) {
				R.setErrorString("Row and Column dimensions do not match.");
				R.setErrorCode(Reader::OTHER);
				return;
			}


			_n = i;
			if(!R.expectChildTag()) return;

			R.traverseChild();
			if(!R.expectTagName("field")) return;
			R.upToParent();

			if(!R.getNextChild()) {
				R.setErrorString("Got a matrix with a field and no data.");
				R.setErrorCode(Reader::OTHER);
				return;
			}
			if(!R.expectChildTag()) return;

			R.traverseChild();
			if(!R.expectTagName("scalar") || !R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagNum(_v))
				return;
			R.upToParent();

			R.upToParent();

			return;
			
		}

		ScalarMatrix(const ScalarMatrix<Field> &M) : _F(M._F)
		{
			_n = M._n;
			_v = M._v;
		}

#endif


		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires time linear in n, the size of the matrix.
		 */
                template<class OutVector, class InVector>
		OutVector& apply(OutVector &y, const InVector &x) const 
		{
			typename VectorTraits<InVector>::VectorCategory t;
			return _app (y, x, t);
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires time linear in n, the size of the matrix.
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector &y, const InVector &x) const
			{ return apply(y, x); }  // symmetric matrix.

		size_t rowdim(void) const { return _n; }
    
		size_t coldim(void) const { return _n; }

		const Field& field() const {return _F;}

#ifdef __LINBOX_XMLENABLED

		ostream &write(ostream &) const;
		bool toTag(Writer &) const;
#endif



	    protected:

		Field _F;   // Field for arithmetic

		size_t _n;  // Number of rows and columns of square matrix.

		Element _v; // the scalar used in applying matrix.

		// dense vector _app for apply
		template<class OutVector, class InVector, class VectorTrait>
		OutVector& _app (OutVector &y, const InVector &x, VectorCategories::DenseVectorTag<VectorTrait>) const;
		// The third argument is just a device to let overloading determine the method.

		// sparse sequence vector _app for apply

		
		template <class OutVector, class InVector, class VectorTrait>
		OutVector& _app (OutVector &y, const InVector &x, VectorCategories::SparseSequenceVectorTag<VectorTrait>) const;

		// sparse associative vector _app for apply
		template<class OutVector, class InVector, class VectorTrait>
		OutVector& _app (OutVector &y, const InVector &x, VectorCategories::SparseAssociativeVectorTag<VectorTrait>) const;

	}; // template <Field, Vector> class ScalarMatrix
   
	// dense vector _app
	template <class Field>
	template <class OutVector, class InVector, class VectorTrait>
	inline OutVector &ScalarMatrix<Field>
		::_app(OutVector& y, const InVector& x, VectorCategories::DenseVectorTag<VectorTrait> t) const
		{   
		    linbox_check (x.size() >= _n);
		    linbox_check (y.size() >= _n);
		    typename OutVector::iterator y_iter = y.begin ();

		    if (_F.isZero(_v)) // just write zeroes
		        for ( ; y_iter != y.end ();  ++y_iter) *y_iter = _v;
                    else if (_F.isOne(_v) ) // just copy 
			copy(x.begin(), x.end(), y.begin());
		    else // use actual muls
		    {   typename InVector::const_iterator x_iter = x.begin ();
		            for (  ; y_iter != y.end () ; ++y_iter, ++x_iter )
		                _F.mul (*y_iter, _v, *x_iter);
		    }
		    return y;

		} // dense vector _app

		
	// sparse sequence vector _app
	template <class Field>
	template <class OutVector, class InVector, class VectorTrait>
	inline OutVector &ScalarMatrix<Field>
		::_app(OutVector& y, const InVector& x, VectorCategories::SparseSequenceVectorTag<VectorTrait> t) const
	{
		//linbox_check ((!x.empty ()) && (_n < x.back ().first));
		// neither is required of x ?

		y.clear (); // we'll overwrite using push_backs.

		// field element to be used in calculations
		Element entry;
		_F.init (entry, 0); // needed?

		// For each element, multiply input element with corresponding element
		// of stored scalar and insert non-zero elements into output vector
		for ( typename InVector::const_iterator x_iter = x.begin (); x_iter != x.end (); ++x_iter) 
		{	_F.mul (entry, _v, x_iter->second);
			if (!_F.isZero (entry)) y.push_back (make_pair (x_iter->first, entry));
		} 

		return y;
	} // sparse sequence vector _app

	// sparse associative vector _app
	template <class Field>
		template <class OutVector, class InVector, class VectorTrait>
	inline OutVector& ScalarMatrix<Field>
	::_app(OutVector& y, const InVector& x, VectorCategories::SparseAssociativeVectorTag<VectorTrait> t) const
	{
		y.clear (); // we'll overwrite using inserts

		// create field elements and size_t to be used in calculations
		Element entry;
		_F.init (entry, 0);

		// Iterator over indices of input vector.
		// For each element, multiply input element with 
		// stored scalar and insert non-zero elements into output vector
		for ( typename InVector::const_iterator x_iter = x.begin (); x_iter != x.end (); ++x_iter)
		{	_F.mul (entry, _v, x_iter->second);
			if (!_F.isZero (entry)) y.insert (y.end (), make_pair (xiter->first, entry));
		}

		return y;
	} // sparse associative vector _app

#ifdef __LINBOX_XMLENABLED
	
	template<class Field, class Vector>
	ostream &ScalarMatrix<Field, Vector>::write(ostream &out) const
	{
		Writer W;
		if( toTag(W) ) 
			W.write(out);

		return out;
	}

	template<class Field, class Vector>
	bool ScalarMatrix<Field, Vector>::toTag(Writer &W) const
	{

		string s;
		W.setTagName("MatrixOver");
		W.setAttribute("rows", Writer::numToString(s, _n));
		W.setAttribute("cols", s);
		W.setAttribute("implDetail", "scalar");

		W.addTagChild();
		_F.toTag(W);
		W.upToParent();

		W.addTagChild();
		W.setTagName("scalar");
		W.addNum( _v);
		W.upToParent();
		
		return true;
	}
#endif	
		
	       
		   



} // namespace LinBox

#endif // __ScalarMatrix
