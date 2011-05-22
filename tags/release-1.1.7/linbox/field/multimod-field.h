/* linbox/field/multimod-field.h
 * Copyright (C) 2005 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */




#ifndef __LINBOX_multimod_field_H
#define __LINBOX_multimod_field_H


#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/field/field-traits.h"
#include "linbox/util/field-axpy.h"
#include "linbox/util/debug.h"
#include <math.h>
#include <vector>




// Namespace in which all LinBox code resides
namespace LinBox { 
	
	class MultiModDouble;
	class MultiModRandIter;

	template <class Ring>
	struct ClassifyRing;
	
	template <>
	struct ClassifyRing<MultiModDouble>{
		typedef RingCategories::ModularTag categoryTag;
	};
	

	
	class MultiModDouble : public FieldInterface {
		
	protected:
		
		std::vector<Modular<double> >              _fields;		
		size_t                                       _size;
		std::vector<integer>                 _crt_constant;
		std::vector<double >                  _crt_inverse;
		integer                                _crt_modulo; 
	

	public:	       
		
		friend class FieldAXPY<MultiModDouble>;
		friend class DotProductDomain<MultiModDouble>;
		friend class MultiModRandIter;

		typedef std::vector<double>              Element;
		typedef MultiModRandIter                RandIter;
		
		MultiModDouble (): _size(0) {}

		MultiModDouble (const std::vector<integer> &primes) : _fields(primes.size()), _size(primes.size()),
								      _crt_constant(primes.size()), _crt_inverse(primes.size()) 
		{
			_crt_modulo=1;
			for (size_t i=0; i<_size; ++i){
				_fields[i]   = Modular<double> (primes[i]);
				_crt_modulo *= primes[i];
			}
			double tmp;
			for (size_t i=0; i<_size; ++i){
				_crt_constant[i]= _crt_modulo/(integer)primes[i];	
				_fields[i].init(tmp, _crt_constant[i]);
				_fields[i].inv(_crt_inverse[i],tmp);
			}

		}

		
		MultiModDouble (const std::vector<double> &primes) : _fields(primes.size()), _size(primes.size()),
								     _crt_constant(primes.size()), _crt_inverse(primes.size()) 
		{
			_crt_modulo=1;
			for (size_t i=0; i<_size; ++i){
				_fields[i]   = Modular<double> (primes[i]);
				_crt_modulo *= primes[i];
			}
			double tmp;
			for (size_t i=0; i<_size; ++i){
				_crt_constant[i]= _crt_modulo/(integer)primes[i];	
				_fields[i].init(tmp, _crt_constant[i]);
				_fields[i].inv(_crt_inverse[i],tmp);
			}		
		}


		MultiModDouble(const MultiModDouble& F) : _fields(F._fields), _size(F._size),
							  _crt_constant(F._crt_constant), _crt_inverse(F._crt_inverse),
							  _crt_modulo(F._crt_modulo) {}

		const MultiModDouble &operator=(const MultiModDouble &F) {
			_fields       = F._fields;
			_size         = F._size;
			_crt_constant = F._crt_constant;
			_crt_modulo   = F._crt_modulo;
			_crt_inverse  = F._crt_inverse;
			return *this;
		}
		
	
		size_t size() const 
		{return this->_size;}
	
		const Modular<double>& getBase(size_t i) const 
		{ return this->_fields[i]; }

		double getModulo(size_t i) const 
		{ return this->_fields[i].modulus;}
			

		const integer& getCRTmodulo() const {return _crt_modulo;}	       

		integer &cardinality (integer &c) const{ 
			c=1;
			for (size_t i=0; i<_size; ++i){
				c*=_fields[i].modulus;
			}				
			return c;
		}

		integer &characteristic (integer &c) const {
			return c=integer(0);
		}

		integer &convert (integer &x, const Element &y) const { 
			x=0;
			double tmp;
			for (size_t i=0;i<_size; ++i){
				_fields[i].mul(tmp, y[i], _crt_inverse[i]);
				integer res= tmp;
				x= x + ( res*_crt_constant[i]);
				if (x > _crt_modulo)
					x-= _crt_modulo;
			}
			return x;
		}

		//double &convert (double &x, const Element& y) const {
		//	return x=y;
		//}
		
		std::ostream &write (std::ostream &os) const {
		       os << "multimod double (";
		       for (size_t i=0;i<_size-1;++i)
			       os<<integer(_fields[i].modulus)<<",";
		       os<<integer(_fields[_size-1].modulus)<<")\n";
		       return os;
		}
		
		std::istream &read (std::istream &is) {
			/*
			is >> modulus; 
			if(modulus <= 1) 
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			if(modulus > 94906265) 
				throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
			*/
			return is;
		}
		
		std::ostream &write (std::ostream &os, const Element &x) const {
			os<<"(";
			for (size_t i=0;i<x.size()-1;++i)
				os<<x[i]<<",";
			os<<x[x.size()-1]<<")\n";
			return os;
		}

		std::istream &read (std::istream &is, Element &x) const {
			integer tmp;
			is >> tmp;
			init(x,tmp); 
			return is;
		}
		

		Element &init (Element &x, const integer &y) const  {
			x.resize(_size);
			for (size_t i=0;i<_size; ++i)
				_fields[i].init(x[i], y);	
			return x;
		}

		inline Element& init(Element& x, double y =0) const {		  
			x.resize(_size);
			for (size_t i=0;i<_size; ++i)
				_fields[i].init(x[i], y);
			return x;
		}

		
		
		inline Element& assign(Element& x, const Element& y) const {
			return x = y;
		}
		
		
		inline bool areEqual (const Element &x, const Element &y) const {
			return x == y;
		}

		inline  bool isZero (const Element &x) const {
			return x == std::vector<double> (_size, 0.);
		}
		
		inline bool isOne (const Element &x) const {
			return x == std::vector<double>(_size, 1.); 
		}

		inline Element &add (Element &x, const Element &y, const Element &z) const {
			for (size_t i=0;i<_size;++i) {
				_fields[i].add(x[i], y[i], z[i]);
			}
			return x;
		}
 
		inline Element &sub (Element &x, const Element &y, const Element &z) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].sub(x[i], y[i], z[i]);
			}			
			return x;
		}
		
		inline Element &mul (Element &x, const Element &y, const Element &z) const {		
			for (size_t i=0;i<_size;++i){
				_fields[i].mul(x[i], y[i], z[i]);
			}
			return x;
		}
 
		inline Element &div (Element &x, const Element &y, const Element &z) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].div(x[i], y[i], z[i]);
			}
			return x;
		}
 
		inline Element &neg (Element &x, const Element &y) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].neg(x[i], y[i]);
			}
			return x;
		}
 
		inline Element &inv (Element &x, const Element &y) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].inv(x[i], y[i]);
			}
			return x;		  		  
		}

		inline Element &axpy (Element &r, 
				      const Element &a, 
				      const Element &x, 
				      const Element &y) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].axpy(r[i], a[i], x[i], y[i]);
			}				\
			return r;
		}

		inline Element &addin (Element &x, const Element &y) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].addin(x[i], y[i]);
			}
			return x;
		}
 
		inline Element &subin (Element &x, const Element &y) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].subin(x[i], y[i]);
			}	
			return x;
		}
 
		inline Element &mulin (Element &x, const Element &y) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].mulin(x[i], y[i]);
			}
			return x;
		}
 
		inline Element &divin (Element &x, const Element &y) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].divin(x[i], y[i]);
			}
			return x;
		}
 
		inline Element &negin (Element &x) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].negin(x[i]);
			}
			return x;
		}
		
		inline Element &invin (Element &x) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].invin(x[i]);
			}
			return x;
		}
		
		inline Element &axpyin (Element &r, const Element &a, const Element &x) const {
			for (size_t i=0;i<_size;++i){
				_fields[i].axpyin(r[i], a[i], x[i]);
			}
			return r;
		}
		
		static inline double getMaxModulus()
		{ return 94906265.0; } // floor( 2^26.5 )
		
	};// end of class MultiModField

	class MultiModRandIter {
	
	public:
		
		MultiModRandIter() {}
		
		MultiModRandIter(const MultiModDouble &F,
				 const integer   &size=0,
				 const integer   &seed=0) : _field(F), _size(size), _seed(seed), _randiter(F._size)
		{			
			for (size_t i=0;i< F._size;++i)
				_randiter[i] = new Modular<double>::RandIter(F._fields[i],size,seed) ;
		}
		
		~MultiModRandIter() {for  (size_t i=0;i< _randiter.size();++i) delete _randiter[i]; }
			
		const MultiModRandIter& operator= (const MultiModRandIter &R) {
			_seed  = R._seed;
			_size  = R._size;
			_field = R._field;
			_randiter = R._randiter;				
			return *this;
		}

		std::vector<double>& random(std::vector<double> &x) const { 			
			for (size_t i=0;i<x.size();++i)
				_randiter[i]->random(x[i]);
			return x;
		}						
		
	protected:
		MultiModDouble        _field;
		integer                _size;	
		integer                _seed;
		std::vector<Modular<double>::RandIter*> _randiter;
	};// end of class MultiModRandIter

/*	
	template <>
	class DotProductDomain<MultiModDouble > : private virtual VectorDomainBase<MultiModDouble> {
	private:
		//std::vector<double> _bound;
		std::vector<size_t>  _nmax;
		//double _invmod;
	  
	public:	  
		typedef std::vector<double> Element;
		
		DotProductDomain (const MultiModDouble &F)
			: VectorDomainBase<MultiModDouble > (F) //, _invmod(1./_F.modulus) 
		{
				for (size_t i=0; i<F.size();++i){
					//_bound[i]=  (double) (1<<53 - (int) (_F.getModulo(i)*_F.getModulo(i))))
					_nmax[i] =  (size_t)floor((double(1<<26)* double(1<<26)*2.)/ (_F.getModulo(i) * _F.getModulo(i)));
				}
			}
	  
	protected:
		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const {
			
			for (size_t k=0;k<_F.size();++k){
				double y = 0.;
				double t = 0.;
				if (v1.size() < _nmax[k]) {
					for (size_t i = 0; i< v1.size();++i)
						y += v1[i][k] * v2[i][k] ;				
					y = fmod(y, _F.getModulo(k));
				}
				else{			
					size_t i=0;
					for (;i< v1.size()- _nmax[k] ;i=i+_nmax[k]){
						for (size_t j=i;j<i+_nmax[k];++j)
							y += v1[j][k] * v2[j][k];
						t+=fmod(y, _F.getModulo(k));
						y=0.;							
					}
					for (;i < v1.size();++i)
						y += v1[i][k] * v2[i][k];
					t+=fmod(y, _F.getModulo(k));
					y = fmod(t, _F.getModulo(k));
				}
				res[k]=y;
			}
			return res;
		}

		template <class Vector1, class Vector2>
		inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const {		  
			
			for (size_t k=0;k<_F.size();++k){
				double y = 0.;
				double t =0.;
								
				if (v1.first.size() < _nmax[k]) {
					for (size_t i=0;i<v1.first.size();++i)
						y+= v1.second[i] * v2[v1.first[i]];
					y = fmod(y, _F.getModulo(k));
				}
				else {			
					size_t i=0;
					for (;i< v1.first.size()- _nmax[k] ;i=i+_nmax[k]){
						for (size_t j=i;j<i+_nmax[k];++j)
							y += v1.second[j] * v2[v1.first[j]];
						t+=fmod(y, _F.getModulo(k));
						y=0.;							
					}
					for (;i < v1.first.size();++i)
						y += v1.second[i] * v2[v1.first[i]];
					t+= fmod(y, _F.getModulo(k));
					y = fmod(t, _F.getModulo(k));
				}
			}
			return res;
		}
	};
*/
}



#endif // __LINBOX_multimod_field_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
