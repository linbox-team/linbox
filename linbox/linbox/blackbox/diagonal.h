/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/blackbox/diagonal.h
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

#ifndef __DIAGONAL_H
#define __DIAGONAL_H

#include "linbox/vector/vector-traits.h"
#include "linbox/util/debug.h"
#include "linbox-config.h"

#ifdef __LIBOX_XMLENABLED

#include "linbox/util/xml/linbox-reader.h"
#include "linbox/util/xml/linbox-writer.h"

using LinBox::Reader;
using LinBox::Writer;

#include <string>
#include <iostream>

using std::istream;
using std::ostream;
using std::string;

#endif


// Namespace in which all LinBox library code resides
namespace LinBox
{

	/** @memo Random diagonal matrices are used heavily as preconditioners.
	 * @doc
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
	template <class _Field,
		  class Trait = typename VectorTraits<typename LinBox::Vector<_Field>::Dense>::VectorCategory>
	class Diagonal;
	
 
	// Specialization of diagonal for LinBox dense vectors
	template <class _Field, class VectorTrait>
	class Diagonal<_Field, VectorCategories::DenseVectorTag<VectorTrait> >
	{
	    public:

		typedef _Field Field;
		typedef typename Field::Element    Element;

		Diagonal(const Field F, const std::vector<typename Field::Element>& v);
#ifdef __LIBOX_XMLENABLED
		Diagonal(Reader &);
		Diagonal(const Diagonal<Field, Vector, VectorCategories::DenseVectorTag<VectorTrait> >&);
#endif


		template <class Vector1, class Vector2>
		Vector1 &apply (Vector1 &y, const Vector2 &x) const;

		template <class Vector1, class Vector2>
		Vector1 &applyTranspose (Vector1 &y, const Vector2 &x) const { return apply (y, x); }

		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 
		const Field& field() const{ return _F; }

#ifdef __LIBOX_XMLENABLED
		ostream &write(ostream &) const;
		bool toTag(Writer &) const;
#endif


	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<DenseVectorTag>
   
	// Specialization of diagonal for LinBox sparse sequence vectors
	template <class _Field, class VectorTrait>
	class Diagonal<_Field, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
	{
	    public:

		typedef _Field Field;
		typedef typename Field::Element    Element;

		Diagonal(const Field F, const std::vector<typename Field::Element>& v);
#ifdef __LIBOX_XMLENABLED
		Diagonal(Reader &);
		Diagonal(const Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag<VectorTrait> > &);
#endif
		
		template<class OutVector, class InVector>
		OutVector& apply(OutVector& y, const InVector& x) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const { return apply(y, x); }

		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 
		const Field& field() const {return _F;}

#ifdef __LIBOX_XMLENABLED
		ostream &write(ostream &) const;
		bool toTag(Writer &) const;
#endif


	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<SparseSequenceVectorTag>

	// Specialization of diagonal for LinBox sparse associative vectors
	template <class _Field, class VectorTrait>
	class Diagonal<_Field, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >
	{
	    public:

		typedef _Field Field;

		typedef typename Field::Element    Element;

		Diagonal(const Field F, const std::vector<typename Field::Element>& v);
#ifdef __LIBOX_XMLENABLED
		Diagonal(Reader &);
		Diagonal(const Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >&);
#endif

		template<class OutVector, class InVector>
		OutVector& apply(OutVector& y, const InVector& x) const;

		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector& y, const InVector& x) const { return apply(y, x); } 


		size_t rowdim(void) const { return _n; } 
		size_t coldim(void) const { return _n; } 
		const Field field() const { return _F; }

#ifdef __LIBOX_XMLENABLED
		ostream &write(ostream &) const;
		bool toTag(Writer &) const;
#endif

	    private:

		// Field for arithmetic
		Field _F;

		// Number of rows and columns of square matrix.
		size_t _n;

		// STL vector of field elements used in applying matrix.
		std::vector<Element> _v;
    
	}; // template <Field, Vector> class Diagonal<SparseAssociativeVectorTag>


	// Method implementations for dense vectors
 
	template <class Field, class VectorTrait>
	inline Diagonal<Field, VectorCategories::DenseVectorTag<VectorTrait> >
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}

	template <class Field, class VectorTrait>
	template <class Vector1, class Vector2>
	inline Vector1 &Diagonal<Field, VectorCategories::DenseVectorTag<VectorTrait> >
		::apply (Vector1 &y, const Vector2 &x) const
	{
		linbox_check (_n == x.size ());
 
		// Create iterators for input, output, and stored vectors
		typename std::vector<Element>::const_iterator v_iter;
		typename Vector2::const_iterator x_iter;
		typename Vector1::iterator y_iter;
 
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
 
	template <class Field, class VectorTrait>
	inline Diagonal<Field, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}

	template <class Field, class VectorTrait>
	template<class OutVector, class InVector>
	inline OutVector &Diagonal<Field, VectorCategories::SparseSequenceVectorTag<VectorTrait> >
		::apply(OutVector& y, const InVector& x) const
	{
		linbox_check ((!x.empty ()) && (_n >= x.back ().first));

		y.clear (); // we'll overwrite using push_backs.

		// create field elements and size_t to be used in calculations
		size_t i;
		Element zero, entry;
		_F.init (zero, 0);
		_F.init (entry, 0);

		// Create iterators for input and stored vectors
		typename std::vector<Element>::const_iterator v_iter;
		typename InVector::const_iterator x_iter;
 
		// Start at beginning of _v vector
		v_iter = _v.begin ();
 
		// Iterator over indices of input vector.
		// For each element, multiply input element with corresponding element
		// of stored vector and insert non-zero elements into output vector
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++) {
			i = (*x_iter).first;
			_F.mul (entry, *(v_iter + i), (*x_iter).second);
			if (!_F.isZero (entry)) y.push_back ( std::pair<size_t, Element>(i, entry));
		} // for (x_iter = x.begin (); x_iter != x.end (); x_iter++)

		return y;
	} // Vector& Diagonal<SparseSequenceVectorTag>::apply(Vector& y, const Vector&) const

	// Method implementations for sparse associative vectors
 
	template <class Field, class VectorTrait>
	inline Diagonal<Field, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >
		::Diagonal(const Field F, const std::vector<typename Field::Element>& v)
		: _F(F), _n(v.size()), _v(v)
	{}

	template <class Field, class VectorTrait>
	template<class OutVector, class InVector>
	inline OutVector& Diagonal<Field, VectorCategories::SparseAssociativeVectorTag<VectorTrait> >
		::apply(OutVector& y, const InVector& x) const
	{
		linbox_check ((!x.empty ()) && (_n >= x.rbegin ()->first));

		y.clear (); // we'll overwrite using inserts

		// create field elements and size_t to be used in calculations
		size_t i;
		Element zero, entry;
		_F.init (zero, 0);
		_F.init (entry, 0);

		// Create iterators for input and stored vectors
		typename std::vector<Element>::const_iterator v_iter;
		typename InVector::const_iterator x_iter;
 
		// Start at beginning of _v vector
		v_iter = _v.begin ();
 
		// Iterator over indices of input vector.
		// For each element, multiply input element with corresponding element
		// of stored vector and insert non-zero elements into output vector
		for (x_iter = x.begin (); x_iter != x.end (); x_iter++)
		{
			i = x_iter->first;
			_F.mul (entry, *(v_iter + i), (*x_iter).second);
			if (!_F.isZero (entry)) y.insert (y.end (), std::pair<size_t, Element>(i, entry));
		}

		return y;
	} // Vector& Diagonal<SparseAssociativeVectorTag>::apply(...) const

#ifdef __LIBOX_XMLENABLED

	template<class Field, class Vector, class Trait>
	ostream &Diagonal<Field, Vector, VectorCategories::DenseVectorTag<Trait> >:: write(ostream &out) const
	{
		Writer W;
		if( toTag(W)) 
			W.write(out);

		return out;
	}
	
	template<class Field, class Vector, class Trait>
	Diagonal<Field, Vector, VectorCategories::DenseVectorTag<Trait> >::Diagonal(const Diagonal<Field, Vector, VectorCategories::DenseVectorTag<Trait> >& M) : _F(M._F)
	{
		_n = M._n;
		_v = M._v;
	}


	template<class Field, class Vector, class Trait>
	Diagonal<Field, Vector, VectorCategories::DenseVectorTag<Trait> >::Diagonal(Reader &R) : _F(R.Down(1))
	{
		typedef typename Field::Element Element;
		size_t i, j;
		Element e;
	 
		R.Up(1);
		if(!R.expectTagName("MatrixOver")) return;
		if(!R.expectAttributeNum("rows", i) || !R.expectAttributeNum("cols",j)) {
			return;
		}
		if(i != j) {
			R.setErrorString("Row and column dimensions do not match.");
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
		if(R.checkTagName("scalar") ) {
			if(!R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagNum(e)) return;
			R.upToParent();

			_v.resize(_n, e);
		}
		else if(!R.expectTagName("diag")) 
			return;
		else {
			if(!R.expectTagNumVector(_v)) return;
		}

		R.upToParent();

		return;
	}

	template<class Field, class Vector, class Trait>
	bool Diagonal<Field, Vector, VectorCategories::DenseVectorTag<Trait> >::toTag(Writer &W) const
	{
		string s;
		W.setTagName("MatrixOver");
		W.setAttribute("rows", Writer::numToString(s, _n));
		W.setAttribute("cols", s);
		W.setAttribute("implDetail", "diagonal-dense");

		W.addTagChild();
		_F.toTag(W);
		W.upToParent();

		W.addTagChild();
		W.setTagName("diag");
		W.addNumericalList(_v);

		W.getPrevChild();
		W.upToParent();
		
		return true;
	}

	template<class Field, class Vector, class Trait>
	ostream &Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag<Trait> >:: write(ostream &out) const
	{
		Writer W;
		if( toTag(W)) 
			W.write(out);
		
		return out;
	}

	template<class Field, class Vector, class Trait>
	Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag<Trait> >::Diagonal(const Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag<Trait> >&M) : _F(M._F)
	{
		_n = M._n;
		_v = M._v;
	}

	template<class Field, class Vector, class Trait>
	Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag<Trait> >::Diagonal(Reader &R) : _F(R.Down(1))
	{
		typedef typename Field::Element Element;
		size_t i, j;
		Element e;
	 
		R.Up(1);
		if(!R.expectTagName("MatrixOver")) return;
		if(!R.expectAttributeNum("rows", i) || !R.expectAttributeNum("cols",j)) {
			return;
		}
		if(i != j) {
			R.setErrorString("Row and Column dimensions did not match");
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
		if(R.checkTagName("scalar") ) {
			if(!R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagNum(e)) return;
			R.upToParent();

			_v.resize(_n, e);
		}
		else if(!R.expectTagName("diag")) 
			return;
		else {
			if(!R.expectTagNumVector(_v)) return;
		}
		R.upToParent();

		return;
	}

	template<class Field, class Vector, class Trait>
	bool Diagonal<Field, Vector, VectorCategories::SparseSequenceVectorTag<Trait> >::toTag(Writer &W) const
	{
		string s;
		W.setTagName("MatrixOver");
		W.setAttribute("rows", Writer::numToString(s, _n));
		W.setAttribute("cols", s);
		W.setAttribute("implDetail", "diagonal-sequence");

		W.addTagChild();
		_F.toTag(W);
		W.upToParent();

		W.addTagChild();
		W.setTagName("diag");
		W.addNumericalList(_v);

		W.getPrevChild();
		W.upToParent();
		
		return true;
	}

	template<class Field, class Vector, class Trait>
	ostream &Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag<Trait> >:: write(ostream &out) const
	{
		Writer W;
		if( toTag(W)) 
			W.write(out);

		return out;
	}

	template<class Field, class Vector, class Trait>
	Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag<Trait> >::Diagonal(const Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag<Trait> >&M) : _F(M._F)
	{
		_n = M._n;
		_v = M._v;
	}

	template<class Field, class Vector, class Trait>
	Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag<Trait> >::Diagonal(Reader &R) : _F(R.Down(1))
	{
		typedef typename Field::Element Element;
		size_t i, j;
		Element e;
	 
		R.Up(1);
		if(!R.expectTagName("MatrixOver")) return;
		if(!R.expectAttributeNum("rows", i) || !R.expectAttributeNum("cols",j)) {
			return;
		}
		if(i != j) {
			R.setErrorString("Row and Column dimensions do not match up.");
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
		if(R.checkTagName("scalar") ) {
			if(!R.expectChildTag()) return;
			R.traverseChild();
			if(!R.expectTagNum(e)) return;
			R.upToParent();

			_v.resize(_n, e);
		}
		else if(!R.expectTagName("diag")) 
			return;
		else {
			if(!R.expectTagNumVector(_v)) return;
		}
		R.upToParent();

		return;
	}

	template<class Field, class Vector, class Trait>
	bool Diagonal<Field, Vector, VectorCategories::SparseAssociativeVectorTag<Trait> >::toTag(Writer &W) const
	{
		string s;
		W.setTagName("MatrixOver");
		W.setAttribute("rows", Writer::numToString(s, _n));
		W.setAttribute("cols", s);
		W.setAttribute("implDetail", "diagonal-associative");

		W.addTagChild();
		_F.toTag(W);
		W.upToParent();

		W.addTagChild();
		W.setTagName("diag");
		W.addNumericalList(_v);

		W.getPrevChild();
		W.upToParent();
		
		return true;
	}


#endif



} // namespace LinBox

#endif // __DIAGONAL_H
