/* Copyright (C) 2010 LinBox
 *
 *
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#ifndef __LINBOX_modular_balanced_int_H
#define __LINBOX_modular_balanced_int_H
// This file has been replaced by modular-balanced-int32.h

/* balanced representation for modular<int> field, [-p/2,p/2], p is odd. 
 */

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <linbox/field/field-traits.h>

#ifndef LINBOX_MAX_INT
#define LINBOX_MAX_INT 2147483647
#endif

#ifndef LINBOX_MAX_MODULUS
#define LINBOX_MAX_MODULUS 1073741824

#endif 


// Namespace in which all LinBox code resides
namespace LinBox 
{ 

	template< class Element >
		class ModularBalanced;
	template< class Element >
		class ModularBalancedRandIter;
	
	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<ModularBalanced<Element> >;

	template <>
	struct ClassifyRing<ModularBalanced<int> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	
	/// \ingroup field
	template <>
		class ModularBalanced<int> : public FieldInterface {
		protected:
		int modulus;
		int halfmodulus;
		int nhalfmodulus;
		double modulusinv;

		public:	       

		friend class FieldAXPY<ModularBalanced<int> >;
                friend class DotProductDomain<ModularBalanced<int> >;
			       
		typedef int Element;
		typedef ModularBalancedRandIter<int> RandIter;

		//default modular field,taking 65521 as default modulus
		ModularBalanced () :modulus(65521) {
			modulusinv = 1/(double)65521;
			halfmodulus = (65521 >> 1);
			nhalfmodulus = -halfmodulus;
		}

		ModularBalanced (int value)  : modulus(value) {
			halfmodulus = (modulus >> 1);
			nhalfmodulus = -halfmodulus;
			modulusinv = 1 / ((double) value); 
			if(value <= 1) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			if(value > LINBOX_MAX_MODULUS) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
			if( ! (value % 2) ) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be odd");

		}

		ModularBalanced(const ModularBalanced<int>& mf) : modulus(mf.modulus),halfmodulus(mf.halfmodulus),nhalfmodulus(mf.nhalfmodulus),modulusinv(mf.modulusinv) { }

		const ModularBalanced &operator=(const ModularBalanced<int> &F) {
			modulus = F.modulus;
			halfmodulus = F.halfmodulus;
			nhalfmodulus = F.nhalfmodulus;
			modulusinv = F.modulusinv;

			return *this;
		}

	
		integer &cardinality (integer &c) const{ 
			return c = modulus;
		}

		integer &characteristic (integer &c) const {
			return c = modulus; 
		}

		integer &convert (integer &x, const Element &y) const { 
			if(y >= 0)
				return x = y;
			else 
				return x = y + modulus;
		}
		
		std::ostream &write (std::ostream &os) const {
			return os << "int mod " << modulus;
		}
		
		std::istream &read (std::istream &is) {
			is >> modulus; 
			halfmodulus = modulus/2;
			nhalfmodulus = -halfmodulus;
			modulusinv = 1 /((double) modulus );
                        if(modulus <= 1) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
                        if(modulus > LINBOX_MAX_MODULUS) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
			if( ! (modulus % 2) ) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be oddd");	

			return is;
		}
		
		std::ostream &write (std::ostream &os, const Element &x) const {
			return os << x;
		}

		std::istream &read (std::istream &is, Element &x) const {
			integer tmp;
			is >> tmp;
			init(x,tmp); 
			return is;
                }
		

		Element &init (Element &x, const integer &y) const  {
			x = y % integer (modulus);
			if (x < nhalfmodulus) x += modulus;
			else if (x > halfmodulus) x -= modulus;
			return x;
		}

		inline Element& init(Element& x, int y =0) const {
			x = y % modulus;

			if ( x < nhalfmodulus ) x += modulus;
			else if (x > halfmodulus ) x -= modulus;

			return x;
		}

		inline Element& init(Element& x, long y) const {
			x = y % modulus;
			if ( x < nhalfmodulus ) x += modulus;
			else if ( x > halfmodulus ) x -= modulus;

			return x;
		}
		
		inline Element& assign(Element& x, const Element& y) const {
			return x = y;
		}
									
		
		inline bool areEqual (const Element &x, const Element &y) const {
			return x == y;
		}

		inline  bool isZero (const Element &x) const {
			return x == 0; 
		}
		
		inline bool isOne (const Element &x) const {
			return x == 1; 
		}

		inline Element &add (Element &x, const Element &y, const Element &z) const {
			x = y + z;
			if ( x > halfmodulus ) x -= modulus;
			else if ( x < nhalfmodulus ) x += modulus;

			return x;
		}
 
		inline Element &sub (Element &x, const Element &y, const Element &z) const {
			x = y - z;
			if (x > halfmodulus) x -= modulus;
			else if (x < nhalfmodulus) x += modulus;
			return x;
		}
		
		inline Element &mul (Element &x, const Element &y, const Element &z) const {
			int q;

			q  = (int) ((((double) y) * ((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (int) (y*z - q*modulus);
			
			if (x > halfmodulus)
				x -= modulus;
			else if (x < nhalfmodulus)
				x += modulus;

			return x;
		}
 
		inline Element &div (Element &x, const Element &y, const Element &z) const {
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}
 
		inline Element &neg (Element &x, const Element &y) const {
			return x = -y;
		}
 
		inline Element &inv (Element &x, const Element &y) const {
			int d, t;			
			XGCD(d, x, t, y, modulus);
			if (d != 1) 
				throw PreconditionFailed(__FUNCTION__,__LINE__,"InvMod: inverse undefined");
			if (x > halfmodulus) 
				x -= modulus;
			else if (x < nhalfmodulus)
				x += modulus;

			return x;		
							      
		}

		inline Element &axpy (Element &r, 
				      const Element &a, 
				      const Element &x, 
				      const Element &y) const {
			int q;
			
			q  = (int) (((((double) a) * ((double) x)) + (double)y) * modulusinv);  // q could be off by (+/-) 1
			r = (int) (a * x + y - q*modulus);
			
			
			if (r > halfmodulus)
				r -= modulus;
			else if (r < nhalfmodulus)
				r += modulus;

			return r;

		}

		inline Element &addin (Element &x, const Element &y) const {
			x += y;
			if ( x > halfmodulus ) x -= modulus;
			else if (x < -halfmodulus) x += modulus;

			return x;
		}
 
		inline Element &subin (Element &x, const Element &y) const {
			x -= y;
			if (x > halfmodulus) 
				x -= modulus;
			else if (x < nhalfmodulus)
				x += modulus;

			return x;
		}
 
		inline Element &mulin (Element &x, const Element &y) const {
			return mul(x,x,y);
		}
 
		inline Element &divin (Element &x, const Element &y) const {
			return div(x,x,y);
		}
 
		inline Element &negin (Element &x) const {
			return x = -x;
		}
 
		inline Element &invin (Element &x) const {
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const {
			int q;
			
			q  = (int) (((((double) a)*((double) x)) + (double)r) * modulusinv);  // q could be off by (+/-) 1
			r = (int) (a * x + r - q*modulus);
			
			
			if (r > halfmodulus)
				r -= modulus;
			else if (r < nhalfmodulus)
				r += modulus;

			return r;
		}

		private:

      		inline static void XGCD(int& d, int& s, int& t, int a, int b) {
			int  u, v, u0, v0, u1, v1, u2, v2, q, r;
			
			int aneg = 0, bneg = 0;
			
			if (a < 0) {
				if (a < -LINBOX_MAX_INT) throw PreconditionFailed(__FUNCTION__,__LINE__,"XGCD: integer overflow");
				a = -a;
				aneg = 1;
			}
			
			if (b < 0) {
				if (b < -LINBOX_MAX_INT) throw PreconditionFailed(__FUNCTION__,__LINE__,"XGCD: integer overflow");
				b = -b;
				bneg = 1;
			}
			
			u1 = 1; v1 = 0;
			u2 = 0; v2 = 1;
			u = a; v = b;
			
			while (v != 0) {
				q = u / v;
				r = u % v;
				u = v;
				v = r;
				u0 = u2;
				v0 = v2;
				u2 =  u1 - q*u2;
				v2 = v1- q*v2;
				u1 = u0;
				v1 = v0;
			}
			
			if (aneg)
				u1 = -u1;
			
			if (bneg)
				v1 = -v1;
			
			d = u;
			s = u1;
			t = v1;
		}
		
	};

	template <>
		class FieldAXPY<ModularBalanced<int> > {	  
		public:
	  
		typedef int Element;
		typedef ModularBalanced<int> Field;
	  
		FieldAXPY (const Field &F) : _F (F),_y(0),_times(0) {
		}


		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F), _y (0),_times(0){}
	  
		FieldAXPY<ModularBalanced<int> > &operator = (const FieldAXPY &faxpy) {
			_F = faxpy._F; 
			_y = faxpy._y; 
			_times = faxpy._times;
			return *this; 
		}
	  
		inline int64& mulacc (const Element &a, const Element &x) {
		        int64 t = (int64) a * (int64)   x;
			if (_times < blocksize) {
				++_times;
				return _y += t;
			}

			else {
				_times = 1;
				normalize();
				return _y += t;
			}
		}

		inline int64& accumulate (const Element &t) {
			if (_times < blocksize) {
				++_times;
				return _y += t;
			}

			else {
				_times = 1;
				normalize();
				return _y += t;
			}
		}

		inline Element& get (Element &y) {

			normalize();

			y = _y;
			
			if (y > _F.halfmodulus)
				y -= _F.modulus;
			else if (y < _F.nhalfmodulus)
				y += _F.modulus;

			return y;
		}

		inline FieldAXPY &assign (const Element y) {
			_y = y; 
			return *this;
		}

		inline void reset() {
			_y = 0;
		}

		private:
	  
		Field _F;
		int64 _y;
		int _times;
		static const int blocksize = 32;

		inline void normalize() {
			_y = (int)_y -(int)(int64)((double) _y * _F.modulusinv) * (int)_F.modulus;
		}
			
	};


	template <>
		class DotProductDomain<ModularBalanced<int> > : private virtual VectorDomainBase<ModularBalanced<int> > {

		private:
		const int blocksize;
		
		public:	  
		typedef int Element;	  
		DotProductDomain (const ModularBalanced<int> &F)
			: VectorDomainBase<ModularBalanced<int> > (F) ,blocksize(32){
		}
				
		protected:
		template <class Vector1, class Vector2>
			inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const {
			
			typename Vector1::const_iterator pv1,pv1e;
			typename Vector2::const_iterator pv2;
		  			
			int64 y = 0;
			int64 t;
			int times = 0;

			pv1 = pv1e = v1.begin();
			pv2 = v2.begin();

			for(int i = 0; i < v1.size() / blocksize ;++i) {
				pv1e = pv1e + blocksize;
				for(;pv1 != pv1e;++pv1,++pv2) {
					t = (((int64) *pv1 ) * ((int64) *pv2 ));
					y += t;
				}
				normalize(y);									
			}

			for(;pv1 != v1.end(); ++pv1, ++pv2) {
				t = (((int64) *pv1 ) * ((int64) *pv2 ));
				y += t;
			}
			
			normalize(y);
			res = y;

			if (res > _F.halfmodulus) res -= _F.modulus;
			else if(res < _F.nhalfmodulus) res += _F.modulus;

			return res;

		}
	  
		template <class Vector1, class Vector2>
			inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const {		  

			typename Vector1::first_type::const_iterator i_idx, i_idxe;
			typename Vector1::second_type::const_iterator i_elt;
		  
			int64 y = 0;
			int64 t;

			i_idx = i_idxe = v1.first.begin();
			i_elt = v1.second.begin();

			for(int i = 0; i < v1.first.size() / blocksize ; ++i) {
				i_idxe = i_idxe + blocksize;
				for(;i_idx!= i_idxe;++i_idx, ++i_elt) {
					t = ( (int64) *i_elt ) * ( (int64) v2[*i_idx] );
					y += t;
				}
				normalize(y);
			}


			for(;i_idx!= v1.first.end();++i_idx, ++i_elt) {
				t = ( (int64) *i_elt ) * ( (int64) v2[*i_idx] );
				y += t;
			}

			normalize(y);

			res = y;
			if (res > _F.halfmodulus) res -= _F.modulus;
			else if(res < _F.nhalfmodulus) res += _F.modulus;	

			return res;
		}

		inline void normalize(int64& _y) const {
			_y = (int)_y -(int)(int64)((double) _y * _F.modulusinv) * (int)_F.modulus;
		}
			
	};
}

#include "linbox/randiter/modular-balanced.h"
#endif //__LINBOX_modular_balanced_int_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
