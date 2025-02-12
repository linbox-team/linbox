/* linbox/field/modular-balanced-double.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2005,2008 Clement Pernet
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * and Clement Pernet <Clement.Pernet@imag.fr>
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

/*! @file field/modular/modular-balanced-double.h
 * @ingroup field
 * @brief Balanced representation of <code>Z/mZ</code> over \c double .
 */

#ifndef __LINBOX_modular_balanced_double_H
#define __LINBOX_modular_balanced_double_H

#ifdef __INTEL_COMPILER
#define FmodF fmodf
#else
#define FmodF fmod
#endif

#include "linbox/linbox-config.h"
#include "linbox/ring/modular.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <cmath>
#include "linbox/field/field-traits.h"
#include "linbox/randiter/modular-balanced.h"

#include <givaro/modular-balanced-double.h>


// Namespace in which all LinBox code resides
namespace LinBox
{

	class MultiModDouble;

	//! Specialization  of FieldAXPY.
	template <>
	class FieldAXPY<Givaro::ModularBalanced<double> > {
	public:

		typedef double Element;
		typedef double Abnormal;
		typedef Givaro::ModularBalanced<double> Field;

		FieldAXPY (const Field &F) :
			_field (&F),
			_y(0.) , _bound( (double) ((1ULL << 53) - (unsigned long) (field().characteristic()*field().characteristic())))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field),
			_y(faxpy._y), _bound(faxpy._bound)
		{}

		FieldAXPY<Givaro::ModularBalanced<double> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			_y= faxpy._y;
			_bound= faxpy._bound;
			return *this;
		}

		inline Element& mulacc (const Element &a, const Element &x)
		{
			//                 Element tmp= a*x;
			//                 return accumulate(tmp);
			return accumulate(a*x);
		}

		inline Element& accumulate (const Element &tmp)
		{
			_y += tmp;
			if (_y > _bound)
				return _y = fmod (_y, field().characteristic());
			else
				return _y;
		}
		inline Element& subumulate (const Element &tmp)
		{
			_y -= tmp;
			if (_y < 0)
				return _y += field().characteristic();
			else
				return _y;
		}

		inline Element& get (Element &y) {
			_y = fmod (_y, field().characteristic());
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
				return _y = fmod (_y, field().characteristic());
			else
				return _y;
		}

		inline const Field & field() const { return *_field; }
	private:

		const Field *_field;
		double _y;
		double _bound;
	};


	//! Specialization  of DotProductDomain.
	template <>
	class DotProductDomain<Givaro::ModularBalanced<double> > : public  VectorDomainBase<Givaro::ModularBalanced<double> > {
	private:
		double _bound;
		size_t _nmax;

	public:
		typedef double Element;
		DotProductDomain(){}
		DotProductDomain (const Givaro::ModularBalanced<double> &F) :
			VectorDomainBase<Givaro::ModularBalanced<double> > (F), _bound( (double) ( (1ULL<<53) - (unsigned long) (field().characteristic()*field().characteristic())))
		{
			_nmax= (size_t)floor((double(1<<26)* double(1<<26)*2.)/ (field().characteristic() * field().characteristic()));
		}

		using VectorDomainBase<Givaro::ModularBalanced<double> >::field;
	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				// y = fmod(y, field().characteristic());
				field().init(res, y);
			}
			else{
				size_t i=0;
				double t = 0.;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=fmod(y, field().characteristic());
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=fmod(y, field().characteristic());
				// y = fmod(t, field().characteristic());
				field().init(res, t);
			}
			return res;
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;


			if (v1.first.size() < _nmax) {
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				// y = fmod(y, field().characteristic());
				field().init(res, y);
			}
			else {
				double t =0.;
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmod(y, field().characteristic());
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmod(y, field().characteristic());
				// y = fmod(t, field().characteristic());
				field().init(res, t);
			}
			return res ;
		}
	};
}

#endif //__LINBOX_modular_balanced_double_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
