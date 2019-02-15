/* linbox/field/modular-float.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2007 Clement Pernet
 * Written by Clement Pernet <cpernet@uwaterloo.ca>
 *
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file field/modular/modular-float.h
 * @ingroup field
 * @brief  representation of <code>Z/mZ</code> over \c float .
 */

#ifndef __LINBOX_modular_float_H
#define __LINBOX_modular_float_H

#ifdef __INTEL_COMPILER
#define FmodF fmodf
#else
#define FmodF fmod
#endif

#include <cmath>
#include <givaro/modular-floating.h>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/ring/modular.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include "linbox/field/field-traits.h"

// Namespace in which all LinBox code resides
namespace LinBox
{
	class MultiModFloat;

	template <>
	class FieldAXPY<Givaro::Modular<float> > {
	public:

		typedef float Element;
		typedef float Abnormal;
		typedef Givaro::Modular<float> Field;

		FieldAXPY (const Field &F) :
			_field (&F) , //_invmod(1./field().fcharacteristic()),
			_y(0.) , _bound( (float) ( (1_i32 << 23) - (uint32_t) (field().characteristic()*field().characteristic())))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field),// _invmod(faxpy._invmod) ,
			_y(faxpy._y), _bound(faxpy._bound)
		{}

		inline Element& mulacc (const Element &a, const Element &x)
		{
			Element tmp= a*x;
			return accumulate(tmp);
		}

		inline Element& accumulate (const Element &tmp)
		{
			_y += tmp;
			if (_y > _bound)
				return _y = fmodf (_y, field().fcharacteristic());
			else
				return _y;
		}

		inline Element& get (Element &y)
		{
			_y = fmodf (_y, field().fcharacteristic());
			return y=_y ;
		}

		inline FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		inline void reset()
		{
			_y = 0.;
		}

		inline const Field & field() const { return *_field; }

	private:

		const Field *_field;
		//float _invmod;
		float _y;
		float _bound;
	};


	template <>
	class DotProductDomain<Givaro::Modular<float> > : public VectorDomainBase<Givaro::Modular<float> > {
	private:
		float _bound;
		size_t _nmax;
		//float _invmod;

	public:
		typedef float Element;
		DotProductDomain (const Givaro::Modular<float> &F) :
			VectorDomainBase<Givaro::Modular<float> > (F)
			, _bound( (float) ( (1<<23) - (uint32_t) (F.characteristic()*F.characteristic())))
			//, _invmod(1./field().fcharacteristic())
		{
			_nmax= (size_t)floor((double(1<<11)* double(1<<12))/ double(F.fcharacteristic() * F.fcharacteristic()));
		}

		using VectorDomainBase<Givaro::Modular<float> >::field;
	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			float y = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y = fmodf(y, field().fcharacteristic());
			}
			else {
				float t = 0.;
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax)
				{
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=fmodf(y, field().fcharacteristic());
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=fmodf(y, field().fcharacteristic());
				y = fmodf(t, field().fcharacteristic());
			}
			return res = y;
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			float y = 0.;


			if (v1.first.size() < _nmax)
			{
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmodf(y, field().fcharacteristic());
			}
			else
			{
				float t =0.;
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax)
				{
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmodf(y, field().fcharacteristic());
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmodf(y, field().fcharacteristic());
				y = fmodf(t, field().fcharacteristic());
			}
			return res = y;
		}
	};
}

#undef FmodF

#endif //__LINBOX_modular_float_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
