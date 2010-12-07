/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

/* linbox/field/modular-crooked-double.h
 * Copyright (C) 2010 LinBox
 *
 * adapted from field/modular-balanced-double.h
 * by Brice Boyer <brice.boyer@imag.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_modular_crooked_double_H
#define __LINBOX_modular_crooked_double_H


#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <cmath>
#include "linbox/field/field-traits.h"
#include "linbox/randiter/modular-crooked.h"
#include "linbox/randiter/nonzero.h"


// Namespace in which all LinBox code resides
namespace LinBox
{
	template< class Element >
	class ModularCrooked;

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<ModularCrooked<Element> >;

	template <>
	struct ClassifyRing<ModularCrooked<double> >
	{
		typedef RingCategories::ModularTag categoryTag;
	};

	class MultiModDouble;

	/// \ingroup field
	template <>
	class ModularCrooked<double> : public FieldInterface {

	protected:

		double  modulus;
		double up_mod;
		double lo_mod;
		unsigned long   lmodulus;

	public:

		friend class FieldAXPY<ModularCrooked<double> >;
		friend class DotProductDomain<ModularCrooked<double> >;
		friend class MultiModDouble;

		typedef double Element;
		typedef ModularCrookedRandIter<double> RandIter;
		typedef NonzeroRandIter<ModularCrooked<double>,RandIter> NonZeroRandIter;

		static ClassifyRing <ModularCrooked<double> >::categoryTag
		getCategory()
		{
			return ClassifyRing<ModularCrooked<double> >::categoryTag();
		}

		ModularCrooked () {}

		ModularCrooked (int32 p, float f = 0.5, int exp = 1) :
			modulus((double)p), up_mod( std::ceil((p-1.)*f) ), lo_mod( up_mod-modulus+1 ),lmodulus (p)
		{
			if(modulus <= 1)
				throw PreconditionFailed(__FUNCTION__,
							 __LINE__,
							 "modulus must be > 1");
			if( exp != 1 ) throw PreconditionFailed(__FUNCTION__,
								__LINE__,
								"exponent must be 1");
			integer max;
			if (modulus > (double) FieldTraits<ModularCrooked<double> >::maxModulus(max))
				throw PreconditionFailed (__FUNCTION__,
							  __LINE__,
							  "modulus is too big");
		}

		ModularCrooked (double p, float f = 0.5) :
			modulus((double)p), up_mod( std::ceil((p-1.)*f) ), lo_mod( up_mod-modulus+1 ),lmodulus (p)
		{
			if (modulus <= 1)
				throw PreconditionFailed(__FUNCTION__,
							 __LINE__,
							 "modulus must be > 1");
			integer max;
			if (modulus > (double) FieldTraits<ModularCrooked<double> >::maxModulus(max))
				throw PreconditionFailed (__FUNCTION__,
							  __LINE__,
							  "modulus is too big");
		}

		ModularCrooked (long int p, float f = 0.5) :
			modulus((double)p), up_mod( std::ceil((p-1.)*f) ), lo_mod(  up_mod-modulus+1 ),lmodulus (p)
		{
			if ((double) modulus <= 1)
				throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			integer max;
			if ((double) modulus > (double) FieldTraits<ModularCrooked<double> >::maxModulus(max))
				throw PreconditionFailed (__FUNCTION__,
							  __LINE__,
							  "modulus is too big");
		}

		ModularCrooked (const integer& p, float f = 0.5)  :
			modulus((double)p), up_mod( std::ceil((double)(p-1)*f) ), lo_mod(  up_mod-modulus+1 ),lmodulus (p)
		{
			if(modulus <= 1)
				throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus must be > 1");
			if(modulus > getMaxModulus())
				throw PreconditionFailed(__func__,__FILE__,__LINE__,"modulus is too big");

		}

		ModularCrooked (const ModularCrooked<double>& mf) :
			modulus (mf.modulus)
			,up_mod (mf.up_mod)
			,lo_mod (mf.lo_mod)
			,lmodulus (mf.lmodulus)
		{}

		const ModularCrooked &operator= (const ModularCrooked<double> &F)
		{
			modulus = F.modulus;
			up_mod = F.up_mod;
			lo_mod = F.lo_mod;
			lmodulus= F.lmodulus;
			return *this;
		}


		inline integer &cardinality (integer &c) const
		{
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

		inline double &convert (double &x, const Element& y) const
		{
			return x=y;
		}

		inline float &convert (float &x, const Element& y) const
		{
			return x=y;
		}

		std::ostream &write (std::ostream &os) const
		{
			// os << modulus << '(' << lo_mod << ',' << up_mod << ')' << std::endl;
			os << "crooked double mod " << int(modulus) << " @ " ;
			os.precision(2) ;
			os << (double)up_mod/(modulus-1);
			os.precision();
			return os ;
		}

		std::istream &read (std::istream &is)
		{
			is >> modulus;
			if(modulus <= 1)
				throw PreconditionFailed (__FUNCTION__,
							  __LINE__,
							  "modulus must be > 1");
			if(modulus > getMaxModulus())
				throw PreconditionFailed (__FUNCTION__,
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
			// JGD : should'nt it be double tmp ???
			is >> tmp;
			init(x,tmp);
			return is;
		}


		inline Element &init (Element &x, const integer &y) const  {
			// 			return x = (Element)mpz_fdiv_ui(y.get_mpz(),lmodulus );
			return x = (Element)(y%lmodulus);
		}

		inline Element& init(Element& x, const double y =0) const
		{

			x = drem (y, modulus);
			if (x < lo_mod) return x += modulus;
			if (x > up_mod) x -= modulus;
			return x;
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
			if ( x < lo_mod ) return x += modulus;
			if ( x > up_mod ) x -= modulus;
			return x;
		}

		inline Element &sub (Element &x,
				     const Element &y,
				     const Element &z) const
		{
			x = y - z;
			if (x < lo_mod ) return x += modulus;
			if (x > up_mod ) x -= modulus;
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
			if ( tx < lo_mod ) return x = tx + modulus;
			if ( tx > up_mod ) return x = tx - modulus;
			return x = (double) tx;
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
			if ( x < lo_mod ) return x += modulus;
			if ( x > up_mod ) x -= modulus;
			return x;
		}

		inline Element &subin (Element &x, const Element &y) const
		{
			x -= y;
			if ( x < lo_mod ) return x += modulus;
			if ( x > up_mod ) x -= modulus;
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
			double max_double = (double) (1ULL<<DBL_MANT_DIG) - modulus ;
			double p = std::max(up_mod,-lo_mod) ;
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


		static inline double getMaxModulus()
		{ return 67108864.0; } // 2^26

	};

#define SQR(A) \
	((A)*(A))

	template <>
	class FieldAXPY<ModularCrooked<double> > {
	public:

		typedef double Element;
		typedef ModularCrooked<double> Field;

		FieldAXPY (const Field &F) :
			_F (F), _y(0.) , _bound( (double) ((1ULL << 53) - (int) (SQR(std::max(_F.up_mod,-_F.lo_mod)))))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_F (faxpy._F),
			_y(faxpy._y), _bound(faxpy._bound)
		{}

		FieldAXPY<ModularCrooked<double> > &operator = (const FieldAXPY &faxpy)
		{
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
				return _y = drem (_y, _F.modulus);
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
			_y = drem (_y, _F.modulus);
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
				return _y = drem (_y, _F.modulus);
			else
				return _y;
		}

	private:

		Field _F;
		double _y;
		double _bound;
	};

	template <>
	class DotProductDomain<ModularCrooked<double> > : private virtual VectorDomainBase<ModularCrooked<double> > {
	private:
		double _bound;
		size_t _nmax;

	public:
		typedef double Element;
		DotProductDomain (const ModularCrooked<double> &F) :
			VectorDomainBase<ModularCrooked<double> > (F), _bound( (double) ( (1ULL<<53) - (int) (SQR(std::max(_F.up_mod,-_F.lo_mod)))))
		{
			_nmax= (size_t)floor((double(1<<26)* double(1<<26)*2.)/ (SQR(std::max(_F.up_mod,-_F.lo_mod))));
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
				y = drem(y, _F.modulus);
			}
			else{
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=drem(y, _F.modulus);
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=drem(y, _F.modulus);
				y = drem(t, _F.modulus);
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
				y = drem(y, _F.modulus);
			}
			else {
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=drem(y, _F.modulus);
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= drem(y, _F.modulus);
				y = drem(t, _F.modulus);
			}
			return res = y;
		}
	};

#include <iostream>
	template<class T>
	std::ostream& operator<< (std::ostream & o, const ModularCrooked<T> & F)
	{
		return F.write(o);
	}

}
#endif //__LINBOX_modular_crooked_double_H

