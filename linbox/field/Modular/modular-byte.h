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

/*! @file field/Modular/modular-byte.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c byte .
 */
#ifndef __LINBOX_modular_bit_H
#define __LINBOX_modular_bit_H


#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"

#ifndef LINBOX_MAX_INT8 /* 127 */
#define LINBOX_MAX_INT8 INT8_MAX
#endif

#ifdef __ICC
#pragma warning(disable:2259)
#endif

// Namespace in which all LinBox code resides
namespace LinBox
{

	template<class Element>
	class Modular;

	template<class Element>
	class ModularRandIter;

	template<class Field>
	class FieldAXPY;

	template<class Field>
	class DotProductDomain;

	template<class Field>
	class MVProductDomain;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<Modular<Element> >;

	template <>
	struct ClassifyRing<Modular<int8_t> >{
		typedef RingCategories::ModularTag categoryTag;
	};

	/** \brief Specialization of Modular to signed 8 bit element type with efficient dot product.
	 *
	 * Efficient element operations for dot product, mul, axpy, by using floating point
	 * inverse of modulus (borrowed from NTL) and some use of non-normalized intermediate values.
	 *
	 * Requires: modulus < 2^7.
	 * Intended use: prime modulus < 2^7.
	 \ingroup field
	 */
	template <>
	class Modular<int8_t> : public FieldInterface {
	public:
		typedef int8_t Element;
	protected:
		Element modulus;
		unsigned long lmodulus ;
		double modulusinv;
	public:
		const Element one,zero,mOne;
		friend class FieldAXPY<Modular<Element> >;
		friend class DotProductDomain<Modular<Element> >;
		friend class MVProductDomain<Modular<Element> >;

		typedef ModularRandIter<Element> RandIter;

		//default modular field,taking 65521 as default modulus
		Modular () :
			modulus(13),lmodulus(13)
			,one((Element)1),zero((Element)0),mOne(Element(modulus-(Element)1))
		{
			modulusinv=1/(double)13;
		}

		Modular (int value, int exp = 1)  :
			modulus(Element(value)),lmodulus((unsigned int)value)
			,one((Element)1),zero((Element)0),mOne(Element(modulus-(Element)1))
		{
			modulusinv = 1 / ((double) value);
#ifdef DEBUG
			if(exp != 1) throw PreconditionFailed(__func__,__FILE__,__LINE__,"exponent must be 1");
			if(value <= 1) throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			integer max;
			if(value > FieldTraits< Modular<Element> >::maxModulus(max)) throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
		}

		Modular(const Modular<Element>& mf) :
			modulus(mf.modulus),lmodulus(mf.lmodulus),modulusinv(mf.modulusinv)
			,one(mf.one),zero(mf.zero),mOne(mf.mOne)
		{}

		Modular &operator=(const Modular<Element> &F)
		{
			modulus    = F.modulus;
			lmodulus   = F.lmodulus;
			modulusinv = F.modulusinv;
			F.assign(const_cast<Element&>(one),F.one);
			F.assign(const_cast<Element&>(zero),F.zero);
			F.assign(const_cast<Element&>(mOne),F.mOne);

			return *this;
		}


		inline integer &cardinality (integer &c) const
		{
			return c = modulus;
		}

		inline integer &characteristic (integer &c) const
		{
			return c = modulus;
		}

		inline unsigned long cardinality () const
		{
			return  lmodulus;
		}

		inline unsigned  long characteristic () const
		{
			return  lmodulus;
		}


		inline integer &convert (integer &x, const Element &y) const
		{
			return x = y;
		}

		inline std::ostream &write (std::ostream &os) const
		{
			return os << "Element mod " << (int)modulus;
		}

		inline std::istream &read (std::istream &is)
		{
			int prime;
			is >> prime;
			modulus = (Element) prime;
			modulusinv = 1 /((double) modulus );
#ifdef DEBUG
			if(prime <= 1) throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			integer max;
			if(prime > FieldTraits< Modular<Element> >::maxModulus(max)) throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus is too big");
#endif

			return is;
		}

		inline std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << x;
		}

		inline std::istream &read (std::istream &is, Element &x) const
		{
			integer tmp;
			is >> tmp;
			init(x,tmp);
			return is;
		}


		inline Element &init (Element &x, const integer &y) const
		{
			x =(Element)((int16_t) (y % (long) (modulus)));
			if (x < 0) x=Element(x+modulus);
			return x;
		}

		inline Element& init(Element& x, int y ) const
		{
			x = Element(y % int(modulus));
			if ( x < 0 ) x=Element(x+modulus);
			return x;
		}

		inline Element& init(Element& x, long y) const
		{
			x = Element(y % long(modulus));
			if ( x < 0 ) x=Element(x+modulus);
			return x;
		}

	inline Element& init(Element& x, long unsigned y) const
		{
			x = Element (y % lmodulus);
			if ( x < 0 ) x=Element(x+modulus);
			return x;
		}


	Element& init(Element&x) const
	{
		return x = 0 ;
	}

		inline Element& assign(Element& x, const Element& y) const
		{
			return x=y;
		}


		inline bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		inline  bool isZero (const Element &x) const
		{
			return x == 0;
		}

		inline bool isOne (const Element &x) const
		{
			return x == 1;
		}

		inline Element &add (Element &x, const Element &y, const Element &z) const
		{
			x = Element(y + z);
			if ( (uint8_t)x >= modulus )
				x =Element(( (uint8_t)x )- modulus);
			return x;
		}

		inline Element &sub (Element &x, const Element &y, const Element &z) const
		{
			x = Element(y - z);
			if (x < 0) x=Element(x+modulus);
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			Element q;

			double ab=((double) y)* ((double) z);
			q  = (Element)(ab*modulusinv);  // q could be off by (+/-) 1
			x = (Element) (ab - ((double) q )* ((double) modulus));


			if (x >= modulus)
				x=Element(x-modulus);
			else if (x < 0)
				x=Element(x+modulus);

			return x;
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		inline Element &neg (Element &x, const Element &y) const
		{
			if(y==0) return x=0;
			else return x=Element(modulus-y);
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			Element d, t;
			XGCD(d, x, t, y, modulus);
#ifdef DEBUG
			if (d != 1)
				throw PreconditionFailed(__func__,__FILE__,__LINE__,"InvMod: inverse undefined");
#endif
			if (x < 0)
				return x=Element(x+modulus);
			else
				return x;
		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			Element q;

			double ab = ((double) a)* ((double) x) + y;
			q  = (Element)(ab*modulusinv);  // q could be off by (+/-) 1
			r = (Element) (ab - ((double) q )* ((double) modulus));


			if (r >= modulus)
				r=Element(r-modulus);
			else if (x < 0)
				r=Element(r+modulus);

			return r;

		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x = Element(x+y);
			if ( ((uint8_t) x) >= modulus )
				x = Element( ((uint8_t) x)-modulus );
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x = Element(x-y);
			if (x < 0) x=Element(x+modulus);
			return x;
		}

		inline Element &mulin (Element &x, const Element &y) const
		{
			return mul(x,x,y);
		}

		inline Element &divin (Element &x, const Element &y) const
		{
			return div(x,x,y);
		}

		inline Element &negin (Element &x) const
		{
			if (x == 0) return x;
			else return x = Element(modulus - x);
		}

		inline Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const
		{


			Element q;

			double ab = ((double) a)* ((double) x) + r;
			q  = (Element)(ab*modulusinv);  // q could be off by (+/-) 1
			r = (Element) (ab - ((double) q )* ((double) modulus));


			if (r >= modulus)
				r=Element(r-modulus);
			else if (x < 0)
				r=Element(r+modulus);

			return r;
		}

		static inline Element getMaxModulus()
		{
		       	return INT8_MAX;
		} // 2^7-1


	private:

		static void XGCD(Element& d, Element& s, Element& t, Element a, Element b)
		{
			int32_t u, v, u0, v0, u1, v1, u2, v2, q, r;

			Element aneg = 0, bneg = 0;

			if (a < 0) {
#ifdef DEBUG
				if (a < -LINBOX_MAX_INT8) throw PreconditionFailed(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				a = Element(-a);
				aneg = 1;
			}

			if (b < 0) {
#ifdef DEBUG
				if (b < -LINBOX_MAX_INT8) throw PreconditionFailed(__func__,__FILE__,__LINE__,"XGCD: integer overflow");
#endif
				b = (Element)-b;
				bneg = 1;
			}

			u1 = 1; v1 = 0;
			u2 = 0; v2 = 1;
			u = a; v = b;

			while (v != 0) {
				q = u / v;
				r = u % v;
				u = v;
				v = r;
				u0 = u2;
				v0 = v2;
				u2 =  u1 - q*u2;
				v2 = v1- q*v2;
				u1 = u0;
				v1 = v0;
			}

			if (aneg)
				u1 = -u1;

			if (bneg)
				v1 = -v1;

			d = Element(u);
			s = Element(u1);
			t = Element(v1);
		}

	};

	template <>
	class FieldAXPY<Modular<int8_t> > {
	public:

		typedef int8_t Element;
		typedef Modular<int8_t> Field;

		FieldAXPY (const Field &F) :
			_field (F),_y(0)
		{
		}

		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field), _y (0)
		{}

		FieldAXPY<Modular<int8_t> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y = faxpy._y;

			return *this;
		}

		inline uint64_t& mulacc (const Element &a, const Element &x)
		{
			uint64_t t = ( (uint16_t) a ) * ( (uint16_t) x );
			return _y +=t;
		}

		inline uint64_t& accumulate (const Element &t)
		{
			return _y += (uint64_t)t;
		}

		inline Element& get (Element &y)
		{
			y = Element(_y % (uint64_t) _field.modulus);
			return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = (uint64_t) y;
			return *this;
		}

		inline void reset()
		{
			_y = 0;
		}

	private:

		Field _field;
		uint64_t _y;
		uint8_t _two_64;
	};


	template <>
	class DotProductDomain<Modular<int8_t> > : private virtual VectorDomainBase<Modular<int8_t> > {

	public:
		typedef int8_t Element;
		DotProductDomain (const Modular<int8_t> &F) :
			VectorDomainBase<Modular<int8_t> > (F)
		{ }

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			typename Vector1::const_iterator i;
			typename Vector2::const_iterator j;

			uint64_t y = 0;
			// uint64_t t;

			for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
				y  += ( (uint16_t) *i ) * ( (uint16_t) *j );
			}


			y %= (uint64_t) _field.modulus;

			return res = (Element) y;

		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{
			typename Vector1::first_type::const_iterator i_idx;
			typename Vector1::second_type::const_iterator i_elt;

			uint64_t y = 0;

			for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {
				y += ( (uint16_t) *i_elt ) * ( (uint16_t) v2[*i_idx] );
			}

			y %= (uint64_t) _field.modulus;

			return res = (Element)y;

		}

	};


	template <>
	class MVProductDomain<Modular<int8_t> >
	{
	public:

		typedef int8_t Element;

	protected:
		template <class Vector1, class Matrix, class Vector2>
		inline Vector1 &mulColDense
		(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
		{
			return mulColDenseSpecialized
			(VD, w, A, v, VectorTraits<typename Matrix::Column>::VectorCategory ());
		}

	private:
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const;
		template <class Vector1, class Matrix, class Vector2>
		Vector1 &mulColDenseSpecialized
		(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const;

		mutable std::vector<uint64_t> _tmp;
	};

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<Modular<int8_t> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
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

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint16_t) *k) * ((uint16_t) *j);

				*l += t;

			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<Modular<int8_t> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
	 VectorCategories::SparseSequenceVectorTag) const
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

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint16_t) k->second) * ((uint16_t) *j);

				_tmp[k->first] += t;

			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<Modular<int8_t> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
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

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint16_t) k->second) * ((uint16_t) *j);

				_tmp[k->first] += t;

			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;

		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
	Vector1 &MVProductDomain<Modular<int8_t> >::mulColDenseSpecialized
	(const VectorDomain<Modular<int8_t> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
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

		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)
			{
				t = ((uint16_t) *k_elt) * ((uint16_t) *j);

				_tmp[*k_idx] += t;

			}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;

		return w;
	}

}

#ifdef __ICC
#pragma warning(enable:2259)
#endif

#include "linbox/randiter/modular.h"
#endif //__LINBOX_modular_bit_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

