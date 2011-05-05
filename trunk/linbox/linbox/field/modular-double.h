/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/modular-double.h
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

/*! @file field/modular-double.h
 * @ingroup field
 * @brief Standard representation of <code>Z/mZ</code> over \c double .
 */

#ifndef __LINBOX_modular_double_H
#define __LINBOX_modular_double_H


#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <math.h>
#include <linbox/field/field-traits.h>
#include "linbox/randiter/nonzero.h"
#include "linbox/randiter/modular.h"

#include "fflas-ffpack/field/modular-double.h"


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
	struct ClassifyRing<Modular<double> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	class MultiModDouble;

	/*! \ingroup modular
	 * Standard representation of \f$\mathbf{Z}/m\mathbf{Z}\f$.
	 * If \c m is the modulus, then elements are represented in \f[ \left
	 * \llbracket 0, m-1  \right \rrbracket.\f]
	 */
	template <>
	class Modular<double> :
	      public FFPACK::Modular<double>,public FieldInterface	{
	      public:
		      typedef double Element;

	      protected:

	      public:
		      friend class FieldAXPY<Modular<Element> >;
		      friend class DotProductDomain<Modular<Element> >;
		      friend class MultiModDouble;

	      public:

		      typedef ModularRandIter<Element> RandIter;
		      typedef NonzeroRandIter<Modular<Element>, ModularRandIter<Element> > NonZeroRandIter;

		      static ClassifyRing<Modular<Element> >::categoryTag getCategory() {return ClassifyRing<Modular<Element> >::categoryTag();}

		      Modular (const integer& p, int e=1) :
			      FFPACK::Modular<double>((unsigned long) p)
		      {
			      linbox_check(e==1);
#ifdef DEBUG
			      if(modulus <= 1)
				      throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			      if(modulus > 94906265)
				      throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus is too big");
#endif
		      }

		      Modular () : FFPACK::Modular<double>() {};

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
		      unsigned long characteristic()const{return FFPACK::Modular<double>::characteristic();}
		      unsigned long cardinality()const{return FFPACK::Modular<double>::cardinality();}
		      template<class T>T&init(T&x)const{return init(x,0);}

		      Element &init (Element &x, const integer &y) const
		      {
			      x = (Element)(y%lmodulus);
			      if (x<0) x+= lmodulus ;
			      linbox_check(x < lmodulus);
			      linbox_check(!(x < 0));
			      return x  ;
		      }

		       bool isMinusOne (const Element &x) const
		      {
			      return (x == modulus-1.);
		      }

		      /** Max number of operations before reducing
		       * @param r if \c r=0, we consider how many \c += are performable.
		       * if \c r=1, then we look for the maximum \c axpy operations doable.
		       * @return \p 0 if the field is too big, a positive number otherwise, \p -1 if infinity
		       * on general fields, it is \p 1.
		       */
		      unsigned long AccBound(const Element r) const
		      {
			      Element one, zero ; init(one,1UL) ; init(zero,0UL);
			      Element max_Element = (Element) (1ULL<<DBL_MANT_DIG) - modulus ; /* other wise 2^52+(2^52-1) */
			      Element p = modulus-1 ;
			      if (areEqual(zero,r))
				      return (unsigned long) (Element(max_Element)/p) ;
			      else if (areEqual(one,r))
			      {
				      if (modulus>= getMaxModulus())
					      return 0 ;
				      else
					      return (unsigned long) (Element(max_Element)/(modulus*modulus)) ;
			      } else
				      throw LinboxError("Bad input, expecting 0 or 1");
			      return 0;
		      }

	      };

	template <>
	class FieldAXPY<Modular<double> > {
	public:

		typedef double Element;
		typedef Modular<double> Field;

		FieldAXPY (const Field &F) :
			_F (F) , //_invmod(1./_F.modulus),
			_y(0.) , _bound( (double) ((1ULL << 53) - (int) (_F.modulus*_F.modulus)))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_F (faxpy._F),// _invmod(faxpy._invmod) ,
			_y(faxpy._y), _bound(faxpy._bound)
		{}

		FieldAXPY<Modular<double> > &operator = (const FieldAXPY &faxpy)
		{
			_F = faxpy._F;
			//_invmod= faxpy._invmod;
			_y= faxpy._y;
			_bound= faxpy._bound;
			return *this;
		}

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
				return _y = fmod (_y, _F.modulus);
			else
				return _y;
		}

		 Element& subumulate (const Element &tmp)
		{
			_y -= tmp;
			if (_y < 0)
				return _y += _F.modulus;
			else
				return _y;
		}

		 Element& get (Element &y)
		{
			_y = fmod (_y, _F.modulus);
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
				return _y = fmod (_y, _F.modulus);
			else
				return _y;
		}

	private:

		Field _F;
		//double _invmod;
		double _y;
		double _bound;
	};

	template <>
	class DotProductDomain<Modular<double> > : private virtual VectorDomainBase<Modular<double> > {
	private:
		double _bound;
		size_t _nmax;
		//double _invmod;

	public:
		typedef double Element;
		DotProductDomain (const Modular<double> &F) :
			VectorDomainBase<Modular<double> > (F), _bound( (double) ( (1ULL<<53) - (int) (F.modulus*F.modulus)))//, _invmod(1./_F.modulus)
		{
			_nmax= (size_t)floor((double(1<<26)* double(1<<26)*2.)/ (F.modulus * F.modulus));
			_nmax = (_nmax>0?_nmax:1);
		}

	protected:
		template <class Vector1, class Vector2>
		 Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;
			double t = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y = fmod(y, _F.modulus);
			}
			else{
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=fmod(y, _F.modulus);
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=fmod(y, _F.modulus);
				y = fmod(t, _F.modulus);
			}
			return res = y;
		}

		template <class Vector1, class Vector2>
		 Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;
			double t =0.;


			if (v1.first.size() < _nmax) {
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmod(y, _F.modulus);
			}
			else {
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmod(y, this->_F.modulus);
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmod(y, _F.modulus);
				y = fmod(t, _F.modulus);
			}
			return res = y;
		}
	};
}


#endif //__LINBOX_modular_double_H

