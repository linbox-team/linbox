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
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========

 * -------------------------------
 */

#ifndef __LINBOX_scalar_H
#define __LINBOX_scalar_H

#include <algorithm>
#include <iostream>
#include "linbox/field/hom.h"
#include "linbox/vector/vector-traits.h"
#include "linbox/linbox-config.h"
#include "linbox/blackbox/blackbox-interface.h"
#include "linbox/solutions/solution-tags.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/util/write-mm.h"

namespace LinBox
{

	/** \brief Blackbox for <tt>aI</tt>.  Use particularly for representing <tt>0</tt> and <tt>I</tt>.

	 * This is a class of blackbox square scalar matrices.
	 * Each scalar matrix occupies O(scalar-size) memory.
	 * The matrix itself is not stored in memory, just the scalar and the dimensions.
	 * \ingroup blackbox
	 */
	template <class Field_>
	class ScalarMatrix : public  BlackboxInterface {
	public:

		typedef Field_ Field;
		typedef typename Field::Element        Element;
		typedef ScalarMatrix<Field> Self_t;

		/*  In each specialization, I must define suitable constructor(s) and
		 *  BlackboxArchetype<Vector> * clone() const;
		 *  Vector& apply(Vector& y, Vector& x) const;
		 *  Vector& applyTranspose(Vector& y, Vector& x) const;
		 *  size_t rowdim(void) const;
		 *  size_t coldim(void) const;
		 *  Field& field() const;
		 *  ...rebind...
		 */

		/// Constructs an initially 0 by 0 matrix.
		//! @bug this should not be allowed (unknown field)
		ScalarMatrix ()	:
			field_(NULL),
			n_(0)
		{}

		ScalarMatrix( MatrixStream<Field> & ms) :
			field_(&ms.field())
			,n_(0)
		{
			size_t c, i, j;
			if( !ms.getDimensions(n_, c) || c != n_ )
				throw ms.reportError(__FUNCTION__,__LINE__);
			ms.nextTriple(i, j, v_);
			if (i != j) throw ms.reportError(__FUNCTION__,__LINE__);
			// finalize();
		}

		void changeField(const Field &F)
		{
			field_ = &F ;
		}

		/** Constructor of readable scalar matrix.
		 * @param F	field in which to do arithmetic.
		 */
		ScalarMatrix (const Field &F) :
			field_(&F),
			n_(0)
		{}

#if 0
		/** Scalar matrix Constructor from an element.
		 * @param F	field in which to do arithmetic.
		 * @param n	size of the matrix.
		 * @param s	scalar, a field element, to be used as the diagonal of the matrix.
		 * @bug this is a wrong constructor, should be the following...
		 */
		ScalarMatrix (const Field &F, const size_t n, const Element &s) :
			field_(&F), n_(n), v_(s)
		{}
#endif

		ScalarMatrix (const Field &F, const size_t n, const size_t m, const Element &s) :
			field_(&F), n_(n), v_(s)
		{
			linbox_check(n ==m);
		}

		ScalarMatrix (const Field &F, const size_t n, const size_t m) :
			field_(&F), n_(n), v_(0)
		{
			linbox_check(m==n);
		}

		/** Constructor from a random element.
		 * @param F    field in which to do arithmetic.
		 * @param n    size of the matrix.
		 * @param iter Random iterator from which to get the diagonal scalar element.
		 */
		ScalarMatrix (const Field &F, const size_t n, const typename Field::RandIter& iter) :
			field_(&F), n_(n)
		{ iter.random(v_); }

		ScalarMatrix(const ScalarMatrix<Field> &Mat) :
			field_(Mat.field_)
			, n_(Mat.n_), v_(Mat.v_)
		{
			//n_ = Mat.n_;
			//v_ = Mat.v_;
		}

		void setScalar(Element & x)
		{
			field().assign(v_, x) ;
		}


		/** Application of BlackBox matrix.
		 * y= A*x.
		 * Requires time linear in n, the size of the matrix.
		 */
		template<class OutVector, class InVector>
		OutVector& apply(OutVector &y, InVector &x) const
		{
			//typename VectorTraits<InVector>::VectorCategory t;
			//return _app (y, x, t);
			return _app (y, x, VectorCategories::DenseVectorTag());
		}

		/** Application of BlackBox matrix transpose.
		 * y= transpose(A)*x.
		 * Requires time linear in n, the size of the matrix.
		 */
		template<class OutVector, class InVector>
		OutVector& applyTranspose(OutVector &y, InVector &x) const
		{ return apply(y, x); }  // symmetric matrix.


		template<typename _Tp1>
		struct rebind {
			typedef ScalarMatrix<_Tp1> other;

			void operator() (other & Ap, const Self_t& A)
			{
				Hom<typename Self_t::Field, _Tp1> hom(A.field(), Ap.field());
				typename _Tp1::Element e;
				Ap.field().assign(e,Ap.field().zero);
				hom.image (e, A.v_);
				Ap.setScalar(e);
			}
		};

		template<typename _Tp1>
		ScalarMatrix (const ScalarMatrix<_Tp1>& S, const Field &F) :
			field_(&F), n_(S.rowdim())
		{
			typename ScalarMatrix<_Tp1>::template rebind<Field>() (*this, S);
		}


		size_t rowdim(void) const { return n_; }

		size_t coldim(void) const { return n_; }

		const Field& field() const {return *field_;}

		// for specialized solutions

		Element& trace(Element& t) const
		{	Element n; field().init(n, n_);
			return field().mul(t, v_, n);
		}

		Element& getEntry(Element& x, const size_t i, const size_t j) const
		{
			// return (i==j ? field().assign(x,v_) : field().assign(x,field().zero));
			return (i==j ? field().assign(x,v_) : field().assign(x,field().zero));
		}

		Element& det(Element& d) const
		{
			return pow(field(), d, v_, n_);
		}

		long int& rank(long int& r) const
		{
			return r = (field().isZero(v_) ? 0 : n_);
		}

		Element& getScalar(Element& x) const { return this->field().assign(x,this->v_); }
		Element& setScalar(const Element& x) { return this->field().assign(this->v_,x); }
		std::ostream& write(std::ostream& os) const {
			writeMMCoordHeader(os, *this, 1, "ScalarMatrix");
			field().write(os << "1 1 ", v_) << std::endl;
			return os;
		}

		std::istream& read(std::istream& is) {
			MatrixStream<Field> ms(field(), is);
			size_t c, i, j;
			if( !ms.getDimensions(n_, c) || c != n_ )
				throw ms.reportError(__FUNCTION__,__LINE__);
			ms.nextTriple(i, j, v_);
			if (i != j) throw ms.reportError(__FUNCTION__,__LINE__);
			return is;
		}

	protected:

		const Field *field_;   // Field for arithmetic

		size_t n_;  // Number of rows and columns of square matrix.

		Element v_; // the scalar used in applying matrix.

		// dense vector _app for apply
		template<class OutVector, class InVector>
		OutVector& _app (OutVector &y, const InVector &x, VectorCategories::DenseVectorTag) const;
		// The third argument is just a device to let overloading determine the method.

		// sparse sequence vector _app for apply


		template <class OutVector, class InVector>
		OutVector& _app (OutVector &y, const InVector &x, VectorCategories::SparseSequenceVectorTag) const;

		// sparse associative vector _app for apply
		template<class OutVector, class InVector>
		OutVector& _app (OutVector &y, const InVector &x, VectorCategories::SparseAssociativeVectorTag) const;

		// p <- a^e.  Really should be a field op
		Element& pow(Field& F, Element& p, const Element& a, const size_t e) {
			Element x; F.init(x);
			if (e == 0) return F.assign(p, F.one);
			if (e%2 == 0) return pow(F, p, F.mul(x, a, a), e/2);
			else /* (e%2 == 1)*/ return F.mul(p, a, pow(F, p, a, e-1));
		}
	}; // template <Field> class ScalarMatrix

	// dense vector _app
	template <class Field>
	template <class OutVector, class InVector>
	inline OutVector &ScalarMatrix<Field>::
	_app(OutVector& y, const InVector& x, VectorCategories::DenseVectorTag t) const
	{
		linbox_check (x.size() >= n_);
		linbox_check (y.size() >= n_);
		typename OutVector::iterator y_iter = y.begin ();

		if (field().isZero(v_)) // just write zeroes
			for ( ; y_iter != y.end ();  ++y_iter) *y_iter = v_;
		else if (field().isOne(v_) ) // just copy
			std::copy(x.begin(), x.end(), y.begin());
		else // use actual muls
		{   typename InVector::const_iterator x_iter = x.begin ();
			for (  ; y_iter != y.end () ; ++y_iter, ++x_iter )
				field().mul (*y_iter, v_, *x_iter);
		}
		return y;

	} // dense vector _app


	// sparse sequence vector _app
	template <class Field>
	template <class OutVector, class InVector>
	inline OutVector &ScalarMatrix<Field>::
	_app(OutVector& y, const InVector& x, VectorCategories::SparseSequenceVectorTag t) const
	{
		//linbox_check ((!x.empty ()) && (n_ < x.back ().first));
		// neither is required of x ?

		y.clear (); // we'll overwrite using push_backs.

		// field element to be used in calculations
		Element entry;
		field().assign(entry, field().zero);

		// For each element, multiply input element with corresponding element
		// of stored scalar and insert non-zero elements into output vector
		for ( typename InVector::const_iterator x_iter = x.begin (); x_iter != x.end (); ++x_iter)
		{	field().mul (entry, v_, x_iter->second);
			if (!field().isZero (entry)) y.push_back (make_pair (x_iter->first, entry));
		}

		return y;
	} // sparse sequence vector _app

	// sparse associative vector _app
	template <class Field>
	template <class OutVector, class InVector>
	inline OutVector& ScalarMatrix<Field> ::
	_app(OutVector& y, const InVector& x, VectorCategories::SparseAssociativeVectorTag t) const
	{
		y.clear (); // we'll overwrite using inserts

		// create field elements and size_t to be used in calculations
		Element entry;
		field().assign(entry, field().zero);

		// Iterator over indices of input vector.
		// For each element, multiply input element with
		// stored scalar and insert non-zero elements into output vector
		for ( typename InVector::const_iterator x_iter = x.begin (); x_iter != x.end (); ++x_iter)
		{	field().mul (entry, v_, x_iter->second);
			if (!field().isZero (entry)) y.insert (y.end (), make_pair (x_iter->first, entry));
		}

		return y;
	} // sparse associative vector _app

	// let solutions know we have getEntry() and trace().
	template <class Field>
	struct GetEntryCategory<ScalarMatrix<Field> >
	{ typedef SolutionTags::Local Tag; };

	template <class Field>
	struct TraceCategory<ScalarMatrix<Field> >
	{ typedef SolutionTags::Local Tag; };

	template <class Field>
	struct DetCategory<ScalarMatrix<Field> >
	{ typedef SolutionTags::Local Tag; };

	template <class Field>
	struct RankCategory<ScalarMatrix<Field> >
	{ typedef SolutionTags::Local Tag; };

} // namespace LinBox

#endif // __LINBOX_scalar_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
