#ifndef __LINBOX_MODULAR__SHORT_H
#define __LINBOX_MODULAR__SHORT_H


#include "linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"

#ifndef LINBOX_MAX_SHORT
#define LINBOX_MAX_SHORT 32767
#endif

#ifndef LINBOX_MAX_SHORT_MODULUS
#define LINBOX_MAX_SHORT_MODULUS 32767
#endif

// Namespace in which all LinBox code resides
namespace LinBox 
{ 

	template<class Element>
		class Modular;
	
	template<class Element>
		class ModularRandIter;

	template <>
		class Modular<short> : public FieldInterface {
		protected:
		short modulus;
		double modulusinv;
		public:	       
		friend class FieldAXPY<Modular<short> >;
                friend class DotProductDomain<Modular<short> >;
			       
		typedef short Element;
		typedef ModularRandIter<short> RandIter;

		//default modular field,taking 65521 as default modulus
		Modular () :modulus(251){modulusinv=1/(double)251;}

		Modular (int value)  : modulus(value) {
			modulusinv = 1 / ((double) value); 
			if(value<=1) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			if(value>LINBOX_MAX_SHORT_MODULUS) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
		}

		Modular(const Modular<short>& mf) : modulus(mf.modulus),modulusinv(mf.modulusinv){}

		const Modular &operator=(const Modular<short> &F) {
			modulus = F.modulus;
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
			return x = y;
		}
		
		std::ostream &write (std::ostream &os) const {
			return os << "short mod " << modulus;
		}
		
		std::istream &read (std::istream &is) {
			int prime;
			is >> prime; 
			modulus=prime;
			modulusinv = 1 /((double) modulus );
			if(prime<=1) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			if(prime>LINBOX_MAX_SHORT_MODULUS) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
		
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
			if (x < 0) x += modulus;
			return x;
		}

		inline Element& init(Element& x, int y =0) const {
			x = y % modulus;
			if ( x < 0 ) x += modulus;
			return x;
		}

		inline Element& init(Element& x, long y) const {
			x = y % modulus;
			if ( x < 0 ) x += modulus;
			return x;
		}
		
		inline Element& assign(Element& x, const Element& y) const {
			return x=y;
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
			if ( (uint16)x >= (uint16)modulus ) x =( (uint16)x )- modulus;
			return x;
		}
 
		inline Element &sub (Element &x, const Element &y, const Element &z) const {
			x = y - z;
			if (x < 0) x += modulus;
			return x;
		}
		
		inline Element &mul (Element &x, const Element &y, const Element &z) const {
			short q;

			double ab=((double) y)* ((double) z);		
			q  = (short)(ab*modulusinv);  // q could be off by (+/-) 1
			x = (short) (ab - ((double) q )* ((double) modulus));
			
			
			if (x >= modulus)
				x -= modulus;
			else if (x < 0)
				x += modulus;

			return x;
		}
 
		inline Element &div (Element &x, const Element &y, const Element &z) const {
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}
 
		inline Element &neg (Element &x, const Element &y) const {
			if(y==0) return x=0;
			else return x=modulus-y;
		}
 
		inline Element &inv (Element &x, const Element &y) const {
			short d, t;			
			XGCD(d, x, t, y, modulus);
			if (d != 1)
				throw PreconditionFailed(__FUNCTION__,__LINE__,"InvMod: inverse undefined");
			if (x < 0)
				return x += modulus;
			else
				return x;
							      
		}

		inline Element &axpy (Element &r, 
				      const Element &a, 
				      const Element &x, 
				      const Element &y) const {
			short q;

			double ab=((double) a)* ((double) x) + ( double ) y;		
			q  = (short)(ab*modulusinv);  // q could be off by (+/-) 1
			r = (short) (ab - ((double) q )* ((double) modulus));
			
			
			if (r >= modulus)
				r -= modulus;
			else if (r < 0)
				r += modulus;

			return r;

		}

		inline Element &addin (Element &x, const Element &y) const {
			x += y;
			if ( ((uint16) x) >= (uint16)modulus ) x = ((uint16) x)-modulus;
			return x;
		}
 
		inline Element &subin (Element &x, const Element &y) const {
			x -= y;
			if (x < 0) x += modulus;
			return x;
		}
 
		inline Element &mulin (Element &x, const Element &y) const {
			return mul(x,x,y);
		}
 
		inline Element &divin (Element &x, const Element &y) const {
			return div(x,x,y);
		}
 
		inline Element &negin (Element &x) const {
			if (x == 0) return x; 
			else return x = modulus - x; 
		}
 
		inline Element &invin (Element &x) const {
			return inv (x, x);
		}

		inline Element &axpyin (Element &r, const Element &a, const Element &x) const {
			
			short q;

			double ab=((double) a)* ((double) x) + ( double ) r;		
			q  = (short)(ab*modulusinv);  // q could be off by (+/-) 1
			r = (short) (ab - ((double) q )* ((double) modulus));
			
			
			if (r >= modulus)
				r -= modulus;
			else if (r < 0)
				r += modulus;

			return r;
		}

		private:

      		static void XGCD(short& d, short& s, short& t, short a, short b) {
			short  u, v, u0, v0, u1, v1, u2, v2, q, r;
			
			short aneg = 0, bneg = 0;
			
			if (a < 0) {
				if (a < -LINBOX_MAX_SHORT) throw PreconditionFailed(__FUNCTION__,__LINE__,"XGCD: integer overflow");
				a = -a;
				aneg = 1;
			}
			
			if (b < 0) {
				if (b < -LINBOX_MAX_SHORT) throw PreconditionFailed(__FUNCTION__,__LINE__,"XGCD: integer overflow");
				b = -b;
				bneg = 1;
			}
			
			u1=1; v1=0;
			u2=0; v2=1;
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
		class FieldAXPY<Modular<short> > {	  
		public:
	  
		typedef short Element;
		typedef Modular<short> Field;
	  
		FieldAXPY (const Field &F) : _F (F),_y(0) {
			uint32 two_64 = 2;
		  
			for (int i = 0; i < 6; ++i)
				two_64 = (two_64 * two_64) % (uint32)_F.modulus;
		  
			_two_64 = (uint16)two_64;
		}

		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F), _y (0),_two_64(faxpy._two_64) {}
	  
		FieldAXPY<Modular<short> > &operator = (const FieldAXPY &faxpy) {
			_F = faxpy._F; 
			_y = faxpy._y; 
			_two_64=faxpy._two_64;
			return *this; 
		}
	  
		inline void accumulate (const Element &a, const Element &x) {
			uint64 t = ( (uint32) a ) * ( (uint32) x );
			_y += t;		 
		}

		inline Element& get (Element &y) {
			y =_y % (uint64) _F.modulus;
			return y;
		}

		inline FieldAXPY &assign (const Element y) {
			_y = y; 
			return *this;
		}

		private:
	  
		Field _F;
		uint64 _y;
		uint16 _two_64;
	};


	template <>
		class DotProductDomain<Modular<short> > : private virtual VectorDomainBase<Modular<short> > {
		private:
		uint16 _two_64;

		public:	  
		typedef short Element;	  
		DotProductDomain (const Modular<short> &F)
			: VectorDomainBase<Modular<short> > (F) {
			uint32 two_64 = 2;
		  
			for (int i = 0; i < 6; ++i)
				two_64 = (two_64 * two_64) % (uint32)_F.modulus;
		  
			_two_64 = (uint16)two_64;
		  
		}
	  
	  
		protected:
		template <class Vector1, class Vector2>
			inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const {
		  
			typename Vector1::const_iterator i;
			typename Vector2::const_iterator j;
		  
			uint64 y = 0;
			uint64 t;
		  
			for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
				y  += ( (uint32) *i ) * ( (uint32) *j );
			}
		       
			y %= (uint64) _F.modulus;
		  
			return res = y;
			
		}
	  
		template <class Vector1, class Vector2>
			inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const {		  
			typename Vector1::first_type::const_iterator i_idx;
			typename Vector1::second_type::const_iterator i_elt;
		  
			uint64 y = 0;
		  
			for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {
				y += ( (uint32) *i_elt ) * ( (uint32) v2[*i_idx] );
			}
		
			y %= (uint64) _F.modulus;
		  
			return res = y;
		}
	  
	};
  
} 

#include "linbox/randiter/modular.h"
#endif
