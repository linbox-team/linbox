/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2010 LinBox
 *
 *
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

/*! @file field/modular-int32_t.h
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

#include "fflas-ffpack/field/modular-int32.h"

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

		typedef int32_t Element;
		typedef ModularRandIter<int32_t> RandIter;
		typedef NonzeroRandIter<Modular<int32_t>, ModularRandIter<int32_t> > NonZeroRandIter;

		Modular (integer &p) :
			FFPACK::Modular<int32_t>((unsigned long)p)
		{}

	       	Modular (int32_t value, int32_t exp=1) :
			FFPACK::Modular<int32_t>(value,exp)
		      {}

		 integer &cardinality (integer &c) const
		{
			return c = modulus;
		}

		 integer &characteristic (integer &c) const
		{
		       	return c = modulus;
		}


		 template<class T>T&convert(T&x,const Element&y)const{return x=T(y);}
		 template<class T>T&characteristic(T&x)const{return x=T(lmodulus);}
		 unsigned long characteristic()const{return FFPACK::Modular<int32_t>::characteristic();}
		 unsigned long cardinality()const{return FFPACK::Modular<int32_t>::cardinality();}

		 integer &convert (integer &x, const Element &y) const
		{
			return x = y;
		}


		 Element &init (Element &x, const integer &y) const
		{
			x = Element (y % modulus);
			if (x < 0) x += modulus;
			return x;
		}

		 Element init(Element&x) const { return FFPACK::Modular<int32_t>::init(x) ; }

		unsigned long AccBound(const Element&r) const
		{
			Element one, zero ; init(one,1UL) ; init(zero,0UL);
			double max_double = (double) (INT_MAX) - modulus ;
			double p = modulus-1 ;
			if (areEqual(zero,r))
				return (unsigned long) (max_double/p) ;
			else if (areEqual(one,r))
			{
				if (modulus>= getMaxModulus())
					return 0 ;
				else
					return (unsigned long) max_double/(modulus*modulus) ;
			}
			else
				throw LinboxError("Bad input, expecting 0 or 1");
			return 0;
		}

	private:

	};

	template <>
	class FieldAXPY<Modular<int32_t> > {
	public:

		typedef int32_t Element;
		typedef Modular<int32_t> Field;

		FieldAXPY (const Field &F) :
			_F (F),_y(0)
		{ }


		FieldAXPY (const FieldAXPY &faxpy) :
			_F (faxpy._F), _y (0)
		{}

		FieldAXPY<Modular<int32_t> > &operator = (const FieldAXPY &faxpy)
		{
			_F = faxpy._F;
			_y = faxpy._y;
			return *this;
		}

		 uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = (uint64_t) a * (uint64_t) x;
			_y += t;
			if (_y < t)
				return _y += _F._two64;
			else
				return _y;
		}

		 uint64_t& accumulate (const Element &t)
		{
			_y += t;
			if (_y < (uint64_t)t)
				return _y += _F._two64;
			else
				return _y;
		}

		 Element& get (Element &y)
		{
			y = Element (_y % (uint64_t) _F.modulus);
			return y;
		}

		 FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		 void reset()
		{
			_y = 0;
		}

	protected:
		Field _F;
		uint64_t _y;
	};


	template <>
	class DotProductDomain<Modular<int32_t> > : private virtual VectorDomainBase<Modular<int32_t> > {

	public:
		typedef int32_t Element;
		DotProductDomain (const Modular<int32_t> &F) :
			VectorDomainBase<Modular<int32_t> > (F)
		{}


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
					y += _F._two64;
			}

			y %= (uint64_t) _F.modulus;
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
					y += _F._two64;
			}


			y %= (uint64_t) _F.modulus;

			return res = y;
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

		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i)
		{
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l)
			{
				t = ((uint64_t) *k) * ((uint64_t) *j);

				*l += t;

				if (*l < t)
					*l += VD.field ()._two64;
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

		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint64_t) k->second) * ((uint64_t) *j);

				_tmp[k->first] += t;

				if (_tmp[k->first] < t)
					_tmp[k->first] += VD.field ()._two64;
			}
		}

		typename Vector1::iterator w_j;
		typedef typename Vector1::value_type Element;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = (Element)( *l % VD.field ().modulus );

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

		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

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

		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);

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

