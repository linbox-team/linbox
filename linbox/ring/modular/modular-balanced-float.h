/* field/modular-balanced-float.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2005,2008 Clement Pernet
 * Written by Clement Pernet <clement.pernet@gmail.com>
 *            Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * Modified   Brice Boyer (briceboyer) <boyer.brice@gmail.com>
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

/*! @file field/modular/modular-balanced-float.h
 * @ingroup field
 * @brief Balanced  representation of <code>Z/mZ</code> over \c float .
 */

#ifndef __LINBOX_modular_balanced_float_H
#define __LINBOX_modular_balanced_float_H

#ifdef __INTEL_COMPILER
#define FmodF fmodf
#else
#define FmodF fmod
#endif

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/ring/modular.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <cmath>
#include "linbox/field/field-traits.h"
#include "linbox/randiter/modular-balanced.h"

#include <givaro/modular-balanced-float.h>


// Namespace in which all LinBox code resides
namespace LinBox
{

	class MultiModFloat;

	template <>
	class FieldAXPY<Givaro::ModularBalanced<float> > {
	public:
		typedef float Element;
		typedef float Abnormal;
		typedef Givaro::ModularBalanced<Element> Field;

		FieldAXPY (const Field &F) :
			_field (&F),
			_y(0.) , _bound( (Element) (((1UL << 24) - (unsigned long) (field().characteristic()*field().characteristic()))))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field),
			_y(faxpy._y), _bound(faxpy._bound)
		{}

		FieldAXPY<Givaro::ModularBalanced<Element> > &operator = (const FieldAXPY &faxpy) {
			_field = faxpy._field;
			_y= faxpy._y;
			_bound= faxpy._bound;
			return *this;
		}

		inline Element& mulacc (const Element &a, const Element &x) {
			//                 Element tmp= a*x;
			//                 return accumulate(tmp);
			return accumulate(a*x);
		}

		inline Element& accumulate (const Element &tmp) {
			_y += tmp;
			if (_y > _bound)
				return _y = fmodf (_y, field().characteristic());
			else
				return _y;
		}
		inline Element& subumulate (const Element &tmp) {
			_y -= tmp;
			if (_y < 0)
				return _y += field().characteristic();
			else
				return _y;
		}

		inline Element& get (Element &y) {
			_y =  fmodf (_y, field().characteristic());
			return y=_y ;
		}

		inline FieldAXPY &assign (const Element y) {
			_y = y;
			return *this;
		}

		inline void reset() {
			_y = 0.;
		}

		inline Element& set (const Element &tmp) {
			_y = tmp;
			if (_y > _bound)
				return _y =  fmodf (_y, field().characteristic());
			else
				return _y;
		}

		inline const Field & field() const { return *_field; }

	private:
		const Field *_field;
		Element _y;
		Element _bound;
	};


	template <>
	class DotProductDomain<Givaro::ModularBalanced<float> > : public  VectorDomainBase<Givaro::ModularBalanced<float> > {
	public:
		typedef float Element;
		DotProductDomain(){}
		DotProductDomain (const Givaro::ModularBalanced<Element> &F) :
			VectorDomainBase<Givaro::ModularBalanced<Element> > (F), _bound( (Element) ( (1UL<<24) - (unsigned long) (field().characteristic()*field().characteristic())))
		{
			_nmax= (size_t)floor((Element(1<<11)* Element(1<<11)*2.)/ (field().characteristic() * field().characteristic()));
		}

		using VectorDomainBase<Givaro::ModularBalanced<Element> >::field;
	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			Element y = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y =  fmodf(y, field().characteristic());
			}
			else{
			Element t = 0.;
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+= fmodf(y, field().characteristic());
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+= fmodf(y, field().characteristic());
				y = fmodf(t, field().characteristic());
			}
			//!@bug should not be neccessary (use assign)
			return field().init(res, y);
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			Element y = 0.;


			if (v1.first.size() < _nmax) {
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmodf(y, field().characteristic());
			}
			else {
			Element t =0.;
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmodf(y, field().characteristic());
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmodf(y, field().characteristic());
				y = fmodf(t, field().characteristic());
			}
			//!@bug should not be neccessary (use assign)
			return field().init(res, y);
		}
	private:
		Element _bound;
		size_t _nmax;

	};
} // Namespace LinBox

#include "linbox/randiter/modular-balanced.h"

#undef FmodF

#endif //__LINBOX_modular_balanced_float_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
