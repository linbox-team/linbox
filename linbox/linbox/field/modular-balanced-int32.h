/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* Copyright (C) 2009 LinBox
 * Written by C Pernet
 * updated to compilable condition by <brice.boyer@imag.fr>
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


/*! @file field/modular-balanced-int32_t.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c int32_t .
 */

#ifndef __LINBOX_modular_balanced_int32_H
#define __LINBOX_modular_balanced_int32_H


#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"
#include "linbox/field/modular-int32.h"

#include "fflas-ffpack/field/modular-balanced-int32.h"

#ifndef LINBOX_MAX_INT /* 2147483647 */
#define LINBOX_MAX_INT INT32_MAX
#endif


// Namespace in which all LinBox code resides
namespace LinBox
{

	template< class Element >
	class ModularBalanced;
	template< class Element >
	class ModularBalancedRandIter;
	template< class Field, class RandIter >
	class NonzeroRandIter;


	template <class Ring>
	struct ClassifyRing;

	template<class Element>
	struct ClassifyRing<ModularBalanced<Element> >;

	template<>
	struct ClassifyRing<ModularBalanced<int32_t> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	/// \ingroup field
	template <>
	class ModularBalanced<int32_t> : public FieldInterface,
	      public FFPACK::ModularBalanced<int32_t>	{
	protected:
		// int32_t modulus;
		// int32_t halfmodulus;
		// int32_t nhalfmodulus;
		// double modulusinv;

	public:

		friend class FieldAXPY<ModularBalanced<int32_t> >;
		friend class DotProductDomain<ModularBalanced<int32_t> >;

		typedef int32_t Element;
		typedef ModularBalancedRandIter<int32_t> RandIter;
		typedef NonzeroRandIter<ModularBalanced<int32_t>,RandIter> NonZeroRandIter;

		ModularBalanced(int32_t p, int32_t e=1) :
			FFPACK::ModularBalanced<int32_t>(p,e)
		      {}

		integer &cardinality (integer &c) const
		{
			return c = modulus;
		}

		integer &characteristic (integer &c) const
		{
		       	return c = modulus;
		}

		// this function converts an int to a natural number ?
		integer &convert (integer &x, const Element &y) const
		{
			if(y >= 0)
				return x = y;
			else
				return x = y + modulus;
		}

		Element &init (Element &x, const integer &y) const
		{
			x = y % (long)modulus;

			if (x < nhalfmodulus)
				x += modulus;
			else if (x > halfmodulus)
				x -= modulus;

			return x;
		}


	};

	template <>
	class FieldAXPY<ModularBalanced<int32_t> > {
	public:

		typedef int32_t Element;
		typedef ModularBalanced<int32_t> Field;

		FieldAXPY (const Field &F) :
			_F (F),_y(0),_times(0)
		{ }


		FieldAXPY (const FieldAXPY &faxpy) :
			_F (faxpy._F), _y (0),_times(0)
		{}

		FieldAXPY<ModularBalanced<int32_t> > &operator = (const FieldAXPY &faxpy)
		{
			_F = faxpy._F;
			_y = faxpy._y;
			_times = faxpy._times;
			return *this;
		}

		inline int64_t& mulacc (const Element &a, const Element &x)
		{
			int64_t t = (int64_t) a * (int64_t)   x;
			if (_times < blocksize) {
				++_times;
				return _y += t;
			}

			else {
				_times = 1;
				normalize();
				return _y += t;
			}
		}

		inline int64_t& accumulate (const Element &t)
		{
			if (_times < blocksize) {
				++_times;
				return _y += t;
			}

			else {
				_times = 1;
				normalize();
				return _y += t;
			}
		}

		inline Element& get (Element &y)
		{

			normalize();

			y = _y;

			if (y > _F.halfmodulus)
				y -= _F.modulus;
			else if (y < _F.nhalfmodulus)
				y += _F.modulus;

			return y;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset()
		{
			_y = 0;
		}

	private:

		Field _F;
		int64_t _y;
		int32_t _times;
		static const int32_t blocksize = 32;

		inline void normalize() {
			_y = (int32_t)_y -(int32_t)(int64_t)((double) _y * _F.modulusinv) * (int32_t)_F.modulus;
		}

	};


	template <>
	class DotProductDomain<ModularBalanced<int32_t> > : private virtual VectorDomainBase<ModularBalanced<int32_t> > {

	private:
		const int32_t blocksize;

	public:
		typedef int32_t Element;
		DotProductDomain (const ModularBalanced<int32_t> &F) :
			VectorDomainBase<ModularBalanced<int32_t> > (F) ,blocksize(32)
		{ }

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			typename Vector1::const_iterator pv1,pv1e;
			typename Vector2::const_iterator pv2;

			int64_t y = 0;
			int64_t t;
			// int32_t times = 0;

			pv1 = pv1e = v1.begin();
			pv2 = v2.begin();

			for(size_t i = 0; i < v1.size() / blocksize ;++i) {
				pv1e = pv1e + blocksize;
				for(;pv1 != pv1e;++pv1,++pv2) {
					t = (((int64_t) *pv1 ) * ((int64_t) *pv2 ));
					y += t;
				}
				normalize(y);
			}

			for(;pv1 != v1.end(); ++pv1, ++pv2) {
				t = (((int64_t) *pv1 ) * ((int64_t) *pv2 ));
				y += t;
			}

			normalize(y);
			res = y;

			if (res > _F.halfmodulus) res -= _F.modulus;
			else if(res < _F.nhalfmodulus) res += _F.modulus;

			return res;

		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			typename Vector1::first_type::const_iterator i_idx, i_idxe;
			typename Vector1::second_type::const_iterator i_elt;

			int64_t y = 0;
			int64_t t;

			i_idx = i_idxe = v1.first.begin();
			i_elt = v1.second.begin();

			for(size_t i = 0; i < v1.first.size() / blocksize ; ++i) {
				i_idxe = i_idxe + blocksize;
				for(;i_idx!= i_idxe;++i_idx, ++i_elt) {
					t = ( (int64_t) *i_elt ) * ( (int64_t) v2[*i_idx] );
					y += t;
				}
				normalize(y);
			}


			for(;i_idx!= v1.first.end();++i_idx, ++i_elt) {
				t = ( (int64_t) *i_elt ) * ( (int64_t) v2[*i_idx] );
				y += t;
			}

			normalize(y);

			res = y;
			if (res > _F.halfmodulus) res -= _F.modulus;
			else if(res < _F.nhalfmodulus) res += _F.modulus;

			return res;
		}

		inline void normalize(int64_t& _y) const
		{
			_y = (int32_t)_y -(int32_t)(int64_t)((double) _y * _F.modulusinv) * (int32_t)_F.modulus;
		}

	};
}

#include "linbox/randiter/modular-balanced.h"
#endif //__LINBOX_modular_balanced_int32_H

