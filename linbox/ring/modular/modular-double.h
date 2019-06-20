/* linbox/field/modular-double.h
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
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

/*! @file field/modular/modular-double.h
 * @ingroup field
 * @brief Standard representation of <code>Z/mZ</code> over \c double .
 */

#ifndef __LINBOX_modular_double_H
#define __LINBOX_modular_double_H

#include <cmath>
#include <givaro/modular-floating.h>

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/ring/modular.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"

#include "linbox/util/write-mm.h"


// Namespace in which all LinBox code resides
namespace LinBox
{
	class MultiModDouble;
	
}


// FieldAXPY/DotProductDomain
namespace LinBox
{
	template <>
	class FieldAXPY<Givaro::Modular<double> > {
	public:

		typedef double Element;
		typedef double Abnormal;
		typedef Givaro::Modular<double> Field;

		FieldAXPY (const Field &F) :
			_field (&F) , //_invmod(1./field().fcharacteristic()),
			_y(0.) , _bound( (double) ((1_i64 << 53) - (uint64_t) (field().characteristic()*field().characteristic())))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field),// _invmod(faxpy._invmod) ,
			_y(faxpy._y), _bound(faxpy._bound)
		{}

#if 0
		FieldAXPY<Givaro::Modular<double> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			//_invmod= faxpy._invmod;
			_y= faxpy._y;
			_bound= faxpy._bound;
			return *this;
		}
#endif

		 Element& mulacc (const Element &a, const Element &x)
		{
			//                 Element tmp= a*x;
			//                 return accumulate(tmp);
			return accumulate(a*x);
		}

		 Element& accumulate (const Element &tmp)
		{
			_y += tmp;
			if (_y > _bound)
				return _y = fmod (_y, field().fcharacteristic());
			else
				return _y;
		}

		 Element& subumulate (const Element &tmp)
		{
			_y -= tmp;
			if (_y < 0)
				return _y += field().fcharacteristic();
			else
				return _y;
		}

		 Element& get (Element &y)
		{
			_y = fmod (_y, field().fcharacteristic());
			return y=_y ;
		}

		 FieldAXPY &assign (const Element y)
		{
			_y = y;
			return *this;
		}

		 void reset()
		{
			_y = 0.;
		}

		 Element& set (const Element &tmp)
		{
			_y = tmp;
			if (_y > _bound)
				return _y = fmod (_y, field().fcharacteristic());
			else
				return _y;
		}

		inline const Field & field() const { return *_field; }

	protected:

		const Field *_field;
		//double _invmod;
		double _y;
		double _bound;
	};

	template <>
	class DotProductDomain<Givaro::Modular<double> > : public  VectorDomainBase<Givaro::Modular<double> > {
	private:
		// double _bound; // BB : not used
		size_t _nmax;
		//double _invmod;

	public:
		//DotProductDomain () { /*std::cerr << "DPD-Md def cstor" << std::endl;*/ }
		typedef double Element;
		DotProductDomain (const Givaro::Modular<double> &F) :
			VectorDomainBase<Givaro::Modular<double> > (F)
			// , _bound( (double) ( (1ULL<<53) - (unsigned long int) (F.characteristic()*F.characteristic())))
			, _nmax(0)//, _invmod(1./field().characteristic())
		{
			_nmax= (size_t)floor((double(1<<26)* double(1<<26)*2.)/ (F.fcharacteristic() * F.fcharacteristic()));
			_nmax = (_nmax>0?_nmax:1);
		}

		using VectorDomainBase<Givaro::Modular<double> >::field;

	protected:
		template <class Vector1, class Vector2>
		 Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y = fmod(y, field().fcharacteristic());
			}
			else{
				double t = 0.;
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=fmod(y, field().fcharacteristic());
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=fmod(y, field().fcharacteristic());
				y = fmod(t, field().fcharacteristic());
			}
			return res = y;
		}


		template <class Vector1, class Vector2>
		 Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;

			if (v1.first.size() < _nmax) {
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmod(y, field().fcharacteristic());
			} else {
				double t = 0.;
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmod(y, field().fcharacteristic());
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmod(y, field().fcharacteristic());
				y = fmod(t, field().fcharacteristic());
			}
			return res = y;
		}
	};
}


#endif //__LINBOX_modular_double_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
