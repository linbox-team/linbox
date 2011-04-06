/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/modular-float.h
 * Copyright (C) 2003 Pascal Giorgi
 *               2008 Clement Pernet
 * Written by Clement Pernet <clement.pernet@gmail.com>
 *            Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
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




// Namespace in which all LinBox code resides
namespace LinBox
{

	template< class Element >
	class ModularBalanced;

	template< class Element >
	class ModularBalancedRandIter;

	template <class Field, class RandIter>
	class NonZeroRandIter;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<ModularBalanced<Element> >;

	template <>
	struct ClassifyRing<ModularBalanced<float> >
	{
		typedef RingCategories::ModularTag categoryTag;
	};

	class MultiModFloat;

	/// \ingroup field
	template <>
	class ModularBalanced<float> : public FieldInterface {
	private :

	public :
		typedef float Element;

	protected:

		Element  modulus;
		Element half_mod;
		Element mhalf_mod;
		unsigned long   lmodulus;

	public:
		friend class FieldAXPY<ModularBalanced<Element> >;
		friend class DotProductDomain<ModularBalanced<Element> >;
		friend class MultiModFloat;

		typedef ModularBalancedRandIter<Element> RandIter;
		typedef NonzeroRandIter<ModularBalanced<Element>,ModularBalancedRandIter<Element> > NonZeroRandIter;

		static ClassifyRing <ModularBalanced<Element> >::categoryTag
		getCategory()
		{
			return ClassifyRing<ModularBalanced<Element> >::categoryTag();
		}

		ModularBalanced () {}

		ModularBalanced (int32_t p, int exp = 1) :
			modulus((Element)p),
			half_mod( Element((p-1)/2)),
			mhalf_mod( half_mod-p+1),
			lmodulus (p)
		{
			if(modulus <= 1)
				throw PreconditionFailed(__func__,
							 __LINE__,
							 "modulus must be > 1");
			if( exp != 1 ) throw PreconditionFailed(__func__,
								__LINE__,
								"exponent must be 1");
			integer max;
			if (modulus > (Element) FieldTraits<ModularBalanced<Element> >::maxModulus(max))
				throw PreconditionFailed (__func__,
							  __LINE__,
							  "modulus is too big");
		}

		ModularBalanced (Element p) :
			modulus (p),
			half_mod( Element((p-1)/2)),
			mhalf_mod( half_mod-p+1),
			lmodulus ((unsigned long)p)
		{
			if (modulus <= 1)
				throw PreconditionFailed(__func__,
							 __LINE__,
							 "modulus must be > 1");
			integer max;
			if (modulus > (Element) FieldTraits<ModularBalanced<Element> >::maxModulus(max))
				throw PreconditionFailed (__func__,
							  __LINE__,
							  "modulus is too big");
		}

		ModularBalanced (long int p) :
			modulus((Element)p),
			half_mod( Element((p-1)/2)),
			mhalf_mod( half_mod-p+1),
			lmodulus(p)
		{
			if ((Element) modulus <= 1)
				throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			integer max;
			if ((Element) modulus > (Element) FieldTraits<ModularBalanced<Element> >::maxModulus(max))
				throw PreconditionFailed (__func__,
							  __LINE__,
							  "modulus is too big");
		}

		ModularBalanced (const integer& p) :
			modulus((Element) p),
			half_mod( Element((p-1)/2)),
			mhalf_mod( half_mod-(Element)p+1),
			lmodulus(p)
		{
			if(modulus <= 1)
				throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if(modulus > getMaxModulus())
				throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus is too big");

		}

		ModularBalanced (const ModularBalanced<Element>& mf) :
			modulus (mf.modulus),
			half_mod (mf.half_mod),
			mhalf_mod( mf.mhalf_mod),
			lmodulus (mf.lmodulus) {}

		const ModularBalanced &operator= (const ModularBalanced<Element> &F) {
			modulus = F.modulus;
			half_mod = F.half_mod;
			mhalf_mod = F.mhalf_mod;
			lmodulus= F.lmodulus;
			return *this;
		}


		inline integer &cardinality (integer &c) const{
			return c = integer(modulus);
		}

		inline integer &characteristic (integer &c) const
		{
			return c = integer(modulus);
		}

		inline size_t characteristic () const
		{
			return modulus;
		}

		inline integer &convert (integer &x, const Element &y) const
		{
			if ( y < 0. ) return x = integer (y + modulus) ;
			else return x = integer (y);
		}

		inline Element &convert (Element &x, const Element& y) const
		{ return x=y; }

		inline double &convert (double &x, const Element& y) const
		{ return x=y; }


		std::ostream &write (std::ostream &os) const
		{
			return os << "balanced Element mod " << int(modulus);
		}

		std::istream &read (std::istream &is) {
			is >> modulus;
			if(modulus <= 1)
				throw PreconditionFailed (__func__,
							  __LINE__,
							  "modulus must be > 1");
			if(modulus > getMaxModulus())
				throw PreconditionFailed (__func__,
							  __LINE__,
							  "modulus is too big");
			return is;
		}

		std::ostream &write (std::ostream &os, const Element &x) const
		{
			return os << int(x);
		}

		std::istream &read (std::istream &is, Element &x) const
		{
			integer tmp;
			// JGD : should'nt it be Element tmp ???
			is >> tmp;
			init(x,tmp);
			return is;
		}


		inline Element &init (Element &x, const integer &y) const
		{
			x = (Element)(y%lmodulus);
			if (x > half_mod) return   x -= modulus;
			else if (x < mhalf_mod) return x += modulus;

			return x;
		}

		inline Element& init(Element& x, const Element y =0) const
		{

			x = fmodf (y, modulus);
			if (x > half_mod) return   x -= modulus;
			else if (x < mhalf_mod) return x += modulus;
			else return x;
		}

		inline Element& assign(Element& x, const Element& y) const
		{
			return x = y;
		}

		inline bool areEqual (const Element &x, const Element &y) const
		{
			return x == y;
		}

		inline  bool isZero (const Element &x) const
		{
			return x == 0.;
		}

		inline bool isOne (const Element &x) const
		{
			return x == 1.;
		}

		inline bool isMinusOne (const Element &x) const
		{
			return (x == -1.);
		}

		inline Element &add (Element &x,
				     const Element &y,
				     const Element &z) const
		{
			x = y + z;
			if ( x > half_mod ) return x -= modulus;
			if ( x < mhalf_mod ) return x += modulus;
			return x;
		}

		inline Element &sub (Element &x,
				     const Element &y,
				     const Element &z) const
		{
			x = y - z;
			if (x > half_mod) return x -= modulus;
			if (x < mhalf_mod) return x += modulus;
			return x;
		}

		inline Element &mul (Element &x, const Element &y, const Element &z) const
		{
			x = y * z;
			return init (x,x);
		}

		inline Element &div (Element &x, const Element &y, const Element &z) const
		{
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}

		inline Element &neg (Element &x, const Element &y) const
		{
			return x = -y;
		}

		inline Element &inv (Element &x, const Element &y) const
		{
			// The extended Euclidean algoritm
			int x_int, y_int, q, tx, ty, temp;
			x_int = int (modulus);
			y_int = (y < 0.) ? int(y + modulus) : int(y);
			tx = 0;
			ty = 1;

			while (y_int != 0) {
				// always: gcd (modulus,residue) = gcd (x_int,y_int)
				//         sx*modulus + tx*residue = x_int
				//         sy*modulus + ty*residue = y_int
				q = x_int / y_int; // integer quotient
				temp = y_int; y_int = x_int - q * y_int;
				x_int = temp;
				temp = ty; ty = tx - q * ty;
				tx = temp;
			}
			if (tx > half_mod ) return x = tx - modulus;
			if ( tx < mhalf_mod ) return x = tx + modulus;
			return x = (Element) tx;
		}

		inline Element &axpy (Element &r,
				      const Element &a,
				      const Element &x,
				      const Element &y) const
		{
			r = a * x + y;
			return init (r, r);
		}

		inline Element &addin (Element &x, const Element &y) const
		{
			x += y;
			if ( x > half_mod ) return x -= modulus;
			if ( x < mhalf_mod ) return x += modulus;
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if ( x > half_mod ) return x -= modulus;
			if ( x < mhalf_mod ) return x += modulus;
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
			return x = -x;
		}

		inline Element &invin (Element &x) const
		{
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const
		{
			r += a * x;
			return init (r, r);
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

		static inline Element getMaxModulus()
		{ return 2048.; } // 2^11

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
}

#include "linbox/randiter/modular-balanced.h"
#include "linbox/randiter/nonzero.h"

#undef FmodF

#endif //__LINBOX_modular_balanced_float_H

