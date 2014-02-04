/* Copyright (C) 2010 LinBox
 *
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file field/Modular/modular-int32.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c int32_t .
 */
#ifndef __LINBOX_modular_int32_H
#define __LINBOX_modular_int32_H


#include <math.h>
#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/write-mm.h"

#include <fflas-ffpack/field/modular-int32.h>

#ifndef LINBOX_MAX_INT /* 2147483647 */
#define LINBOX_MAX_INT INT32_MAX
#endif

// Namespace in which all LinBox code resides
namespace LinBox
{

	template< class Element >
	class Modular;
	template< class Element >
	class ModularRandIter;
	template< class Field, class RandIter >
	class NonzeroRandIter;

	template<class Field>
	class DotProductDomain;
	template<class Field>
	class FieldAXPY;
	template<class Field>
	class MVProductDomain;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<Modular<Element> >;

	template <>
	struct ClassifyRing<Modular<int32_t> > {
		typedef RingCategories::ModularTag categoryTag;
	};



	/** \brief Specialization of Modular to int32_t element type with efficient dot product.
	 *
	 * Efficient element operations for dot product, mul, axpy, by using floating point
	 * inverse of modulus (borrowed from NTL) and some use of non-normalized intermediate values.
	 *
	 * For some uses this is the most efficient field for primes in the range from half word
	 * to 2^30.
	 *
	 * Requires: Modulus < 2^30.
	 * Intended use: 2^15 < prime modulus < 2^30.
	 \ingroup field
	 */
	template <>
	class Modular<int32_t> : public FieldInterface ,
	      public ::FFPACK::Modular<int32_t> {

	protected:

	public:

		friend class FieldAXPY<Modular<int32_t> >;
		friend class DotProductDomain<Modular<int32_t> >;
		friend class MVProductDomain<Modular<int32_t> >;

		typedef FFPACK::Modular<int32_t> Father_t;

		typedef int32_t Element;
		typedef ModularRandIter<int32_t> RandIter;

		Modular (integer &p) :
			Father_t((unsigned long)p)
		{
		}

	       	Modular (int32_t value, int32_t exp=1) :
			Father_t(value,exp)
		      {}

		Modular (long value) :
			Father_t(value)
		      {}

		Modular (unsigned long value) :
			Father_t(value)
		      {}

		using Father_t ::cardinality;
		 integer &cardinality (integer &c) const
		{
			return c = modulus;
		}

		 using Father_t ::characteristic;
		 integer &characteristic (integer &c) const
		{
		       	return c = modulus;
		}

		 using Father_t ::convert;
		 integer &convert (integer &x, const Element &y) const
		{
			return x = y;
		}

		 using Father_t ::init;
		 Element &init (Element &x, const integer &y) const
		{
			x = Element (y % lmodulus);
			if (x < 0) x += modulus;
			return x;
		}

		unsigned long AccBound(const Element&r) const
		{
			double max_double = (double) (INT_MAX) - modulus ;
			double p = modulus-1 ;
			if (areEqual(zero,r))
				return (unsigned long) (max_double/p) ;
			else if (areEqual(one,r))
			{
				if (modulus>= getMaxModulus())
					return 0 ;
				else
					return (unsigned long) max_double/(unsigned long)(modulus*modulus) ;
			}
			else
				throw LinboxError("Bad input, expecting 0 or 1");
			return 0;
		}

		/*- Print field as a constructor call.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 * @param  F  optional name to give the field in the description.  IF F is the null string, only the class name is written.
		 * Example: For element type double and modulus 101,
		 * write(os) produces      "Modular< double > ( 101 )"  on os,
		 * write(os, "F") produces "Modular< double > F( 101 )" on os, and
		 * write(os, "") produces  "Modular< double >"          on os.
		 */
		std::ostream &write (std::ostream &os) const
		{
			integer p = cardinality();
			return os << "Modular<" << eltype( Element() ) << " >( " << p << " )";
		}

		std::ostream &write (std::ostream &os, std::string F) const
		{
			os << "Modular<" << eltype( Element() ) << " > "; // class name
			if (F != "") {
				integer p = cardinality();
				os << F << "( " << p << " )"; // show constuctor args
			}
			return os;
		}

		std::ostream &write (std::ostream & os, const Element & x) const
		{
			return Father_t::write(os,x);
		}


	private:

	};

	template <>
	class FieldAXPY<Modular<int32_t> > {
	public:

		typedef int32_t Element;
		typedef int64_t Abnormal;
		typedef Modular<int32_t> Field;

		FieldAXPY (const Field &F) :
			_field (&F),_y(0)
		{ }


		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field), _y (0)
		{}

		FieldAXPY<Modular<int32_t> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;
			return *this;
		}

		 uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (uint64_t) a * (uint64_t) x;
			_y += t;
			if (_y < t) {
				 _y += (uint64_t)field()._two64;
				 return _y ;
			}
			else
				return _y;
		}

		 uint64_t& accumulate (const Element &t)
		{
			_y += (uint64_t) t;
			if (_y < (uint64_t)t)
				return _y += (uint64_t) field()._two64;
			else
				return _y;
		}

		 Element& get (Element &y)
		{
			y = Element (_y % (uint64_t) field().modulus);
			return y;
		}

		 FieldAXPY &assign (const Element y)
		{
			_y = (uint64_t) y;
			return *this;
		}

		 void reset()
		{
			_y = 0;
		}

		inline const Field & field() { return *_field; }

	protected:
		const Field * _field;
		uint64_t _y;
	};


	template <>
	class DotProductDomain<Modular<int32_t> > : public virtual VectorDomainBase<Modular<int32_t> > {

	public:
		typedef int32_t Element;
		DotProductDomain(){}
		DotProductDomain (const Modular<int32_t> &F) :
			VectorDomainBase<Modular<int32_t> > (F)
		{}

		using VectorDomainBase<Modular<int32_t> >::field;

	protected:
		template <class Vector1, class Vector2>
		 Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			typename Vector1::const_iterator i;
			typename Vector2::const_iterator j;

			uint64_t y = 0;
			uint64_t t;

			for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j)
			{
				t = ( (uint64_t) *i ) * ( (uint64_t) *j );
				y += t;

				if (y < t)
					y += (uint64_t) field()._two64;
			}

			y %= (uint64_t) field().modulus;
			return res = Element(y);

		}

		template <class Vector1, class Vector2>
		 Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{
			typename Vector1::first_type::const_iterator i_idx;
			typename Vector1::second_type::const_iterator i_elt;

			uint64_t y = 0;
			uint64_t t;

			for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt)
			{
				t = ( (uint64_t) *i_elt ) * ( (uint64_t) v2[*i_idx] );
				y += t;

				if (y < t)
					y += (uint64_t) field()._two64;
			}


			y %= (uint64_t) field().modulus;

			return res = (Element) y;
		}
	};

	// Specialization of MVProductDomain for int32_t modular field

	template <>
	class MVProductDomain<Modular<int32_t> > {
	public:

		typedef int32_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		 Vector1 &mulColDense
		(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized
			(VD, w, A, v, typename VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
	};

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<Modular<int32_t> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::DenseVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<uint64_t>::iterator l;

		uint64_t t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());

		std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i)
		{
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
			{
				t = ((uint64_t) *k) * ((uint64_t) *j);

				*l += t;

				if (*l < t)
					*l += (uint64_t) VD.field ()._two64;
			}
		}

		typename Vector1::iterator w_j;
		typedef typename Vector1::value_type elements ;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = elements(*l % VD.field ().modulus);

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<Modular<int32_t> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseSequenceVectorTag) const
	{
		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator       i = A.colBegin ();
		typename Vector2::const_iterator        j;
		typename Matrix::Column::const_iterator k;
		std::vector<uint64_t>::iterator         l;

		uint64_t t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());

		std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint64_t) k->second) * ((uint64_t) *j);

				_tmp[k->first] += t;

				if (_tmp[k->first] < t)
					_tmp[k->first] += (uint64_t)VD.field ()._two64;
			}
		}

		typename Vector1::iterator w_j;
		typedef typename Vector1::value_type val_t;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = (val_t)( (int32_t)(*l) % VD.field ().modulus );

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<Modular<int32_t> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseAssociativeVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<uint64_t>::iterator l;

		uint64_t t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());

		std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i)
		{
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
			{
				t = ((uint64_t) k->second) * ((uint64_t) *j);

				_tmp[k->first] += t;

				if (_tmp[k->first] < t)
					_tmp[k->first] += VD.field ()._two64;
			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<Modular<int32_t> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int32_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseParallelVectorTag) const
	{

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());

		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::first_type::const_iterator k_idx;
		typename Matrix::Column::second_type::const_iterator k_elt;
		std::vector<uint64_t>::iterator l;

		uint64_t t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());

		std::fill (_tmp.begin (), _tmp.begin () +(ptrdiff_t) w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i)
		{
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)
			{
				t = ((uint64_t) *k_elt) * ((uint64_t) *j);

				_tmp[*k_idx] += t;

				if (_tmp[*k_idx] < t)
					_tmp[*k_idx] += VD.field ()._two64;
			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;

		return w;
	}


}

#include "linbox/randiter/modular.h"

#endif //__LINBOX_modular_int32_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
