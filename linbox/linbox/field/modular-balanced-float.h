/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* field/modular-balanced-float.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2005,2008 Clement Pernet
 * Written by Clement Pernet <clement.pernet@gmail.com>
 *            Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 * Modified   Brice Boyer <bboyer@imag.fr>
 * ------------------------------------
 *
 * See COPYING for license information.
 */

/*! @file field/modular-balanced-float.h
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
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <math.h>
#include <linbox/field/field-traits.h>
#include "linbox/randiter/modular-balanced.h"
#include "linbox/randiter/nonzero.h"

#include "fflas-ffpack/field/modular-balanced-float.h"


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
	struct ClassifyRing<ModularBalanced<float> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	class MultiModFloat;

	/// \ingroup field
	template <>
	class ModularBalanced<float> : public FieldInterface,
	      public FFPACK::ModularBalanced<float> {
	      public :
		      typedef float Element;

	      protected:

	      public:
		      friend class FieldAXPY<ModularBalanced<Element> >;
		      friend class DotProductDomain<ModularBalanced<Element> >;
		      friend class MultiModFloat;

		      typedef ModularBalancedRandIter<Element> RandIter;
		      typedef NonzeroRandIter<ModularBalanced<Element>, RandIter > NonZeroRandIter;

		      static ClassifyRing <ModularBalanced<Element> >::categoryTag
		      getCategory()
		      {
			      return ClassifyRing<ModularBalanced<Element> >::categoryTag();
		      }

		      ModularBalanced (const integer& p) :
			      FFPACK::ModularBalanced<float>((unsigned long)p)
		      {
#ifdef DEBUG
			      if(modulus <= 1)
				      throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			      if(modulus > getMaxModulus())
				      throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus is too big");
			      // check integer not too big.
#endif

		      }


		      inline integer &cardinality (integer &c) const
		      {
			      return c = integer(modulus);
		      }

		      inline integer &characteristic (integer &c) const
		      {
			      return c = integer(modulus);
		      }

		      long unsigned characteristic(long unsigned int&p)const{return FFPACK::ModularBalanced<float>::characteristic(p) ; }
		      unsigned long characteristic()const{return FFPACK::ModularBalanced<float>::characteristic();}
		      double&convert(double &x,const Element&y)const{return FFPACK::ModularBalanced<float>::convert(x,y) ; }
		      float&convert(float&x,const Element&y)const{return FFPACK::ModularBalanced<float>::convert(x,y) ; }


		      inline integer &convert (integer &x, const Element &y) const
		      {
			      // if ( y < 0. )
				      // return x = integer (y + modulus) ;
			      // else
				      return x = integer (y);
		      }


		      inline Element &init (Element &x, const integer &y) const
		      {
			      x = (Element)(y%lmodulus);
			      if (x > half_mod) return   x -= modulus;
			      else if (x < mhalf_mod) return x += modulus;

			      return x;
		      }

		      inline bool isMinusOne (const Element &x) const
		      {
			      return (x == -1.);
		      }

		      unsigned long AccBound(const Element&r) const
		      {
			      Element one, zero ; init(one,1UL) ; init(zero,0UL);
			      double max_double = (double) (1ULL<<FLT_MANT_DIG) - modulus ;
			      double p = std::max(half_mod,-mhalf_mod) ;
			      if (areEqual(zero,r))
				      return (unsigned long) (double(max_double)/p) ;
			      else if (areEqual(one,r))
			      {
				      if (modulus>= getMaxModulus())
					      return 0 ;
				      else
					      return (unsigned long) (double(max_double)/(p*p)) ;
			      } else
				      throw LinboxError("Bad input, expecting 0 or 1");
			      return 0;
		      }

	      };

	template <>
	class FieldAXPY<ModularBalanced<float> > {
	public:
		typedef float Element;
		typedef ModularBalanced<Element> Field;

		FieldAXPY (const Field &F) :
			_F (F),
			_y(0.) , _bound( (Element) (((1ULL << 24) - (int) (_F.modulus*_F.modulus))))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_F (faxpy._F),
			_y(faxpy._y), _bound(faxpy._bound)
		{}

		FieldAXPY<ModularBalanced<Element> > &operator = (const FieldAXPY &faxpy) {
			_F = faxpy._F;
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
				return _y = fmodf (_y, _F.modulus);
			else
				return _y;
		}
		inline Element& subumulate (const Element &tmp) {
			_y -= tmp;
			if (_y < 0)
				return _y += _F.modulus;
			else
				return _y;
		}

		inline Element& get (Element &y) {
			_y =  fmodf (_y, _F.modulus);
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
				return _y =  fmodf (_y, _F.modulus);
			else
				return _y;
		}

	private:
		Field _F;
		Element _y;
		Element _bound;
	};


	template <>
	class DotProductDomain<ModularBalanced<float> > : private virtual VectorDomainBase<ModularBalanced<float> > {
	public:
		typedef float Element;
		DotProductDomain (const ModularBalanced<Element> &F) :
			VectorDomainBase<ModularBalanced<Element> > (F), _bound( (Element) ( (1ULL<<24) - (int) (_F.modulus*_F.modulus)))
		{
			_nmax= (size_t)floor((Element(1<<11)* Element(1<<11)*2.)/ (_F.modulus * _F.modulus));
		}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			Element y = 0.;
			Element t = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y =  fmodf(y, _F.modulus);
			}
			else{
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+= fmodf(y, _F.modulus);
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+= fmodf(y, _F.modulus);
				y = fmodf(t, _F.modulus);
			}
			return res = y;
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			Element y = 0.;
			Element t =0.;


			if (v1.first.size() < _nmax) {
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmodf(y, _F.modulus);
			}
			else {
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmodf(y, _F.modulus);
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmodf(y, _F.modulus);
				y = fmodf(t, _F.modulus);
			}
			return res = y;
		}
	private:
		Element _bound;
		size_t _nmax;

	};
} // Namespace LinBox

#include "linbox/randiter/modular-balanced.h"
#include "linbox/randiter/nonzero.h"

#undef FmodF

#endif //__LINBOX_modular_balanced_float_H

