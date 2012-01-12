/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
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
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file field/Modular/modular-balanced-double.h
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
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <math.h>
#include "linbox/field/field-traits.h"
#include "linbox/randiter/modular-balanced.h"
#include "linbox/randiter/nonzero.h"

#include <fflas-ffpack/field/modular-balanced-double.h>


// Namespace in which all LinBox code resides
namespace LinBox
{

	template< class Element >
	class ModularBalanced;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<ModularBalanced<Element> >;

	template <>
	struct ClassifyRing<ModularBalanced<double> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	class MultiModDouble;

	/*! \ingroup modular
	 * Centered representation of \f$\mathbf{Z}/m\mathbf{Z}\f$.
	 * If \c m is the modulus, then elements are represented in \f[ \left
	 * \llbracket \left \lceil -\frac{m-1}{2} \right \rceil, \left \lceil
	 * \frac{m-1}{2} \right \rceil \right \rrbracket.\f] This
	 * representation allows more accumulations before a reduction is
	 * necessary, at the cost of a more expensive reduction.
	 */
	template<>
	class ModularBalanced<double> : public FieldInterface,
	      public FFPACK::ModularBalanced<double> {

	      protected:

	      public:
		      friend class FieldAXPY<ModularBalanced<double> >;
		      friend class DotProductDomain<ModularBalanced<double> >;
		      friend class MultiModDouble;

		      typedef double Element;
		      typedef ModularBalancedRandIter<double> RandIter;

		      static ClassifyRing <ModularBalanced<double> >::categoryTag getCategory()
		      {
			      return ClassifyRing<ModularBalanced<double> >::categoryTag();
		      }

		      ModularBalanced (const integer& p) :
			      FFPACK::ModularBalanced<double>((unsigned long)p)
		      {
#ifdef DEBUG
			      if (p > (integer) ULONG_MAX)
				      throw PreconditionFailed(__func__,__FILE__,__LINE__,"prime too big");
			      if(modulus <= 1)
				      throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			      if(modulus > getMaxModulus())
				      throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus is too big");
#endif

		      }

		      integer &cardinality (integer &c) const
		      {
			      return c = integer(modulus);
		      }

		      integer &characteristic (integer &c) const
		      {
			      return c = integer(modulus);
		      }

		      long unsigned characteristic(long unsigned int&p) const { return FFPACK::ModularBalanced<double>::characteristic(p) ; }
		      double & convert(double &x, const Element &y) const { return FFPACK::ModularBalanced<double>::convert(x,y) ; }
		      float & convert(float&x, const Element &y) const { return FFPACK::ModularBalanced<double>::convert(x,y) ; }
		      unsigned long characteristic()const{return FFPACK::ModularBalanced<double>::characteristic();}
		      unsigned long cardinality()const{return FFPACK::ModularBalanced<double>::cardinality();}

		      integer &convert (integer &x, const Element &y) const
		      {
			      return x = integer (y);
		      }

		      Element &init (Element &x, const integer &y) const
		      {
			      x = (Element)(y%lmodulus);
			      if (x<mhalf_mod) return x += modulus ;
			      else if (x>half_mod) return x -= modulus ;
			      return  x ;
		      }

		      Element &init(Element &x) const
		      {
			      return x = 0 ;
		      }

		      //! @bug faux si modulus==2
		      inline bool isMinusOne (const Element &x) const
		      {
			      return (x == -1.);
		      }

		      unsigned long AccBound(const Element&r) const
		      {
			      // Element one, zero ; init(one,1UL) ; init(zero,0UL);
			      double max_double = (double) (1ULL<<DBL_MANT_DIG) - modulus ;
			      double p = std::max(half_mod,-mhalf_mod) ;
			      if (areEqual(zero,r))
				      return (unsigned long) (double(max_double)/p) ;
			      else if (areEqual(one,r))
			      {
				      if (modulus>= getMaxModulus())
					      return 0 ;
				      else
					      return (unsigned long) (double(max_double)/(p*p)) ;
			      }
			      else
				      throw LinboxError("Bad input, expecting 0 or 1");
			      return 0;
		      }

	      };

	//! Specialization  of FieldAXPY.
	template <>
	class FieldAXPY<ModularBalanced<double> > {
	public:

		typedef double Element;
		typedef ModularBalanced<double> Field;

		FieldAXPY (const Field &F) :
			_field (F),
			_y(0.) , _bound( (double) ((1ULL << 53) - (int) (_field.modulus*_field.modulus)))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field),
			_y(faxpy._y), _bound(faxpy._bound)
		{}

		FieldAXPY<ModularBalanced<double> > &operator = (const FieldAXPY &faxpy)
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
				return _y = fmod (_y, _field.modulus);
			else
				return _y;
		}
		inline Element& subumulate (const Element &tmp)
		{
			_y -= tmp;
			if (_y < 0)
				return _y += _field.modulus;
			else
				return _y;
		}

		inline Element& get (Element &y) {
			_y = fmod (_y, _field.modulus);
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
				return _y = fmod (_y, _field.modulus);
			else
				return _y;
		}

	private:

		Field _field;
		double _y;
		double _bound;
	};


	//! Specialization  of DotProductDomain.
	template <>
	class DotProductDomain<ModularBalanced<double> > : private virtual VectorDomainBase<ModularBalanced<double> > {
	private:
		double _bound;
		size_t _nmax;

	public:
		typedef double Element;
		DotProductDomain (const ModularBalanced<double> &F) :
			VectorDomainBase<ModularBalanced<double> > (F), _bound( (double) ( (1ULL<<53) - (int) (_field.modulus*_field.modulus)))
		{
			_nmax= (size_t)floor((double(1<<26)* double(1<<26)*2.)/ (_field.modulus * _field.modulus));
		}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;
			double t = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y = fmod(y, _field.modulus);
			}
			else{
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=fmod(y, _field.modulus);
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=fmod(y, _field.modulus);
				y = fmod(t, _field.modulus);
			}
			return res = y;
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;
			double t =0.;


			if (v1.first.size() < _nmax) {
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmod(y, _field.modulus);
			}
			else {
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmod(y, _field.modulus);
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmod(y, _field.modulus);
				y = fmod(t, _field.modulus);
			}
			return res = y;
		}
	};
}

#endif //__LINBOX_modular_balanced_double_H

