/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/modular-double.h
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */




#ifndef __LINBOX_MODULAR_DOUBLE_H
#define __LINBOX_MODULAR_DOUBLE_H


#include "linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"




// Namespace in which all LinBox code resides
namespace LinBox { 
	
	template< class Element >
	class Modular;
	template< class Element >
	class ModularRandIter;
	
	template <>
	class Modular<double> : public FieldInterface {
	protected:
		double modulus;
		double inv_modulus;
		
	public:	       
		friend class FieldAXPY<Modular<double> >;
		friend class DotProductDomain<Modular<double> >;
		    
			       
		typedef double Element;
		typedef ModularRandIter<double> RandIter;


		Modular (int p)  : modulus((double)p), inv_modulus(1./(double)p) {}

		Modular (const integer& p) : modulus((double) p), inv_modulus(1./(double)p) {}

		Modular(const Modular<double>& mf) : modulus(mf.modulus),inv_modulus(mf.inv_modulus){}

		const Modular &operator=(const Modular<double> &F) {
			modulus = F.modulus;
			inv_modulus = F.inv_modulus;
			return *this;
		}

	
		integer &cardinality (integer &c) const{ 
			return c = integer(modulus);
		}

		integer &characteristic (integer &c) const {
			return c = integer(modulus); 
		}

		integer &convert (integer &x, const Element &y) const { 
			return x = integer(y);
		}

		double &convert (double &x, const Element& y) const {
			return x=y;
		}
		
		std::ostream &write (std::ostream &os) const {
			return os << "int mod " << (int)modulus;
		}
		
		std::istream &read (std::istream &is) {
			is >> modulus; 
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

		inline Element& init(Element& x, double y =0) const {		  
			double tmp=y;
			int sign=0;
			if (tmp < 0.0) {
				tmp=-tmp;
				sign=1;
			}	
		  
			if (tmp >= modulus) 
				tmp -= (modulus * floor( tmp*inv_modulus));
						
			if ((sign)&&(tmp))
				tmp= modulus - tmp;
		  
			return  x=tmp;		       
		}

		
		
		inline Element& assign(Element& x, const Element& y) const {
			return x = y;
		}
									
		
		inline bool areEqual (const Element &x, const Element &y) const {
			return x == y;
		}

		inline  bool isZero (const Element &x) const {
			return x == 0.; 
		}
		
		inline bool isOne (const Element &x) const {
			return x == 1.; 
		}

		inline Element &add (Element &x, const Element &y, const Element &z) const {
			x = y + z;
			if ( x >= modulus ) x -= modulus;
			return x;
		}
 
		inline Element &sub (Element &x, const Element &y, const Element &z) const {
			x = y - z;
			if (x < 0) x += modulus;
			return x;
		}
		
		inline Element &mul (Element &x, const Element &y, const Element &z) const {		
			double tmp= y*z;
			x= tmp - floor(tmp*inv_modulus)*modulus;
		  
			return x;
		}
 
		inline Element &div (Element &x, const Element &y, const Element &z) const {
			Element temp;
			inv (temp, z);
			return mul (x, y, temp);
		}
 
		inline Element &neg (Element &x, const Element &y) const {
			if(y == 0) return x=0;
			else return x = modulus-y;
		}
 
		inline Element &inv (Element &x, const Element &y) const {
			// The extended Euclidean algoritm 
			int x_int, y_int, q, tx, ty, temp;
			x_int = (int) modulus;
			y_int = (int) y;
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
		  
			if (tx < 0) tx += (int)modulus;
		  
			// now x_int = gcd (modulus,residue)
			return x = (double)tx;
		  
		  
		}

		inline Element &axpy (Element &r, 
				      const Element &a, 
				      const Element &x, 
				      const Element &y) const {
			double tmp= a*x+y;
			return r= tmp- floor(tmp*inv_modulus)*modulus; 

		}

		inline Element &addin (Element &x, const Element &y) const {
			x += y;
			if (  x >= modulus ) x -= modulus;
			return x;
		}
 
		inline Element &subin (Element &x, const Element &y) const {
			x -= y;
			if (x < 0.) x += modulus;
			return x;
		}
 
		inline Element &mulin (Element &x, const Element &y) const {
			return mul(x,x,y);
		}
 
		inline Element &divin (Element &x, const Element &y) const {
			return div(x,x,y);
		}
 
		inline Element &negin (Element &x) const {
			if (x == 0.) return x; 
			else return x = modulus - x; 
		}
		
		inline Element &invin (Element &x) const {
			return inv (x, x);
		}
		
		inline Element &axpyin (Element &r, const Element &a, const Element &x) const {
			double tmp=r+a*x;
			return r= tmp- floor(tmp*inv_modulus)*modulus; 
		}
		
	};

	template <>
	class FieldAXPY<Modular<double> > {	  
	public:
	  
		typedef double Element;
		typedef Modular<double> Field;
	  
		FieldAXPY (const Field &F) : _F (F) , _invmod(1./_F.modulus), _y(0.) , _bound( (double) (1>>53 - (int) (_F.modulus*_F.modulus))) {}
	  
		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F), _invmod(faxpy._invmod) , _y(faxpy._y), _bound(faxpy._bound) {}
	  
		FieldAXPY<Modular<double> > &operator = (const FieldAXPY &faxpy) {
			_F = faxpy._F; 
			_invmod= faxpy._invmod;
			_y= faxpy._y;
			_bound= faxpy._bound;
			return *this; 
		}
	  
		inline void accumulate (const Element &a, const Element &x) {
			Element tmp= a*x;		  
			_y += tmp;
	    
			if (_y > _bound) {
				_y-= floor(_y*_invmod)*_F.modulus; 
			}
		}
	  
		inline Element& get (Element &y) {
			_y-= floor(_y*_invmod)*_F.modulus;
			return y=_y ;
		}
	  
		inline FieldAXPY &assign (const Element y) {
			_y = y; 
			return *this;
		}
	  
	private:
	  
		Field _F;
		double _invmod;
		double _y;
		double _bound;		
	};
	
	
	template <>
	class DotProductDomain<Modular<double> > : private virtual VectorDomainBase<Modular<double> > {
	private:
		double _bound;
		double _invmod;
	  
	public:	  
		typedef double Element;	  
		DotProductDomain (const Modular<double> &F)
			: VectorDomainBase<Modular<double> > (F), _bound( (double) (1>>53 - (int) (_F.modulus*_F.modulus))), _invmod(1./_F.modulus) {}
	  
	  
	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const {
	    
			typename Vector1::const_iterator i;
			typename Vector2::const_iterator j;
	    
			double y = 0;
			double t;
	    
			for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
				t = *i * *j ;
				y += t;
	      
				if (y > _bound)
					y-= floor(y*_invmod)*_F.modulus;
	      
			}
	    
			if (y > _F.modulus)
				y-= floor(y*_invmod)*_F.modulus;
	    
			return res = y;
		}

	      

	  
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const {		  
			typename Vector1::first_type::const_iterator i_idx;
			typename Vector1::second_type::const_iterator i_elt;
	    
			double y = 0;
			double t;
	    
			for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {
				t =  *i_elt  *  v2[*i_idx] ;
				y += t;
	      
				if (y > _bound)
					y-= floor(y*_invmod)*_F.modulus;
			}
	      	      
			if (y > _F.modulus)
				y-= floor(y*_invmod)*_F.modulus;
	    
			return res = y;
		}

	};


}
 

#include "linbox/randiter/modular.h"

#endif
