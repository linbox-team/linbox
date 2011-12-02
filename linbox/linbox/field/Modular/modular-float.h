/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/modular-float.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2007 Clement Pernet
 * Written by Clement Pernet <cpernet@uwaterloo.ca>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

/*! @file field/Modular/modular-float.h
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



#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <math.h>
#include "linbox/field/field-traits.h"
#include "linbox/randiter/nonzero.h"

#include <fflas-ffpack/field/modular-float.h>

// Namespace in which all LinBox code resides
namespace LinBox
{

	template< class Element >
	class Modular;
	template< class Element >
	class ModularRandIter;

	template< class Field, class RandIter >
	class NonzeroRandIter;

	template <class Ring>
	struct ClassifyRing;
	template <class Element>
	struct ClassifyRing<Modular<Element> >;
	template <>
	struct ClassifyRing<Modular<float> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	class MultiModFloat;

	/// \ingroup field
	template <>
	class Modular<float> : public FieldInterface,
	      public FFPACK::Modular<float>	{

	      public :
		      typedef float Element;
		      using FFPACK::Modular<float>::one ;
		      using FFPACK::Modular<float>::zero ;
		      using FFPACK::Modular<float>::mOne ;

	      public:
		      friend class FieldAXPY<Modular<Element> >;
		      friend class DotProductDomain<Modular<Element> >;
		      friend class MultiModFloat;

		      typedef ModularRandIter<Element> RandIter;
		      typedef NonzeroRandIter<Modular<Element>, ModularRandIter<Element> > NonZeroRandIter;

		      static ClassifyRing<Modular<Element> >::categoryTag getCategory()
		      {
			      return ClassifyRing<Modular<Element> >::categoryTag();
		      }

		      Modular (const integer& p) :
			      FFPACK::Modular<float>((unsigned long)p)
		      {
#ifdef DEBUG
			      if(modulus <= 1)
				      throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			      integer max;
			      if(modulus > (Element) FieldTraits<Modular<Element> >::maxModulus(max))
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

		      integer &convert (integer &x, const Element &y) const
		      {
			      return x = integer(y);
		      }

		      template<class T>T&convert(T&x,const Element&y)const{return x=T(y);}
		      template<class T>T&characteristic(T&x)const{return x=T(lmodulus);}
		      unsigned long characteristic(void)const{return FFPACK::Modular<float>::characteristic();}
		      unsigned long cardinality(void)const{return FFPACK::Modular<float>::cardinality();}

		      Element &init (Element &x, const integer &y) const
		      {
			      x = (Element)(y%lmodulus);

			      if (x<0) return x+=modulus ;
			      return x;
		      }

		      Element &init(Element &x) const
		      {
			      return x = 0 ;
		      }

		      unsigned long AccBound(const Element&r) const
		      {
			      double max_double = (double) (1ULL<<FLT_MANT_DIG) - modulus ;
			      double p = modulus-1 ;
			      if (areEqual(zero,r))
				      return (unsigned long) (double(max_double)/p) ;
			      else if (areEqual(one,r))
			      {
				      if (modulus>= getMaxModulus())
					      return 0 ;
				      else
					      return (unsigned long) (double(max_double)/(modulus*modulus)) ;
			      }
			      else
				      throw LinboxError("Bad input, expecting 0 or 1");
			      return 0;
		      }

	      };

	template <>
	class FieldAXPY<Modular<float> > {
	public:

		typedef float Element;
		typedef Modular<float> Field;

		FieldAXPY (const Field &F) :
			_field (F) , //_invmod(1./_field.modulus),
			_y(0.) , _bound( (float) ( (1UL << 23) - (unsigned long int) (_field.modulus*_field.modulus)))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field),// _invmod(faxpy._invmod) ,
			_y(faxpy._y), _bound(faxpy._bound)
		{}

#if 0
		FieldAXPY<Modular<float> > &operator = (const FieldAXPY &faxpy)
		{
			_field    = faxpy._field ;
			//_invmod= faxpy._invmod;
			_y    = faxpy._y;
			_bound= faxpy._bound;
			return *this;
		}
#endif

		inline Element& mulacc (const Element &a, const Element &x)
		{
			Element tmp= a*x;
			return accumulate(tmp);
		}

		inline Element& accumulate (const Element &tmp)
		{
			_y += tmp;
			if (_y > _bound)
				return _y = fmodf (_y, _field.modulus);
			else
				return _y;
		}

		inline Element& get (Element &y)
		{
			_y = fmodf (_y, _field.modulus);
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

	private:

		Field _field;
		//float _invmod;
		float _y;
		float _bound;
	};


	template <>
	class DotProductDomain<Modular<float> > : private virtual VectorDomainBase<Modular<float> > {
	private:
		float _bound;
		size_t _nmax;
		//float _invmod;

	public:
		typedef float Element;
		DotProductDomain (const Modular<float> &F) :
			VectorDomainBase<Modular<float> > (F)
			, _bound( (float) ( (1<<23) - (int) (_field.modulus*_field.modulus)))
			//, _invmod(1./_field.modulus)
		{
			_nmax= (size_t)floor((float(1<<11)* float(1<<12))/ (_field.modulus * _field.modulus));
		}

	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			float y = 0.;
			float t = 0.;
			if (v1.size() < _nmax)
			{
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y = fmodf(y, _field.modulus);
			}
			else
			{
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax)
				{
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=fmodf(y, _field.modulus);
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=fmodf(y, _field.modulus);
				y = fmodf(t, _field.modulus);
			}
			return res = y;
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			float y = 0.;
			float t =0.;


			if (v1.first.size() < _nmax)
			{
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmodf(y, _field.modulus);
			}
			else
			{
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax)
				{
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmodf(y, _field.modulus);
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmodf(y, _field.modulus);
				y = fmodf(t, _field.modulus);
			}
			return res = y;
		}
	};
}

#include "linbox/randiter/modular.h"

#undef FmodF

#endif //__LINBOX_modular_float_H

