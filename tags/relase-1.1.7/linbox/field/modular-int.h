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


#warning "use of modular-int.h is now deprecated.  Use modular-int32.h."
#include <linbox/field/modular-int32.h>

#if 0

#ifndef __LINBOX_modular_int_H
#define __LINBOX_modular_int_H

/* Use of this class is depreciated.  Use Modular<int32> in modular-int32.h */

#include "linbox/linbox-config.h"
#include "linbox/integer.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/field/field-interface.h"
#include "linbox/util/debug.h"
#include <math.h>
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
		class Modular;
	template< class Element >
		class ModularRandIter;
	template<class Field>
		class DotProductDomain;
	template<class Field>
		class FieldAXPY;
	template<class Field>
		class MVProductDomain;
	
	template <class Ring>
	struct ClassifyRing; 
	
	template <class Element>
	struct ClassifyRing<Modular<Element> >;

	template <>
	struct ClassifyRing<Modular<int> >{
		typedef RingCategories::ModularTag categoryTag;
	};
	
	/** \brief Specialization of Modular to int element type with efficient dot product.
	 * 
	 * Efficient element operations for dot product, mul, axpy, by using floating point
	 * inverse of modulus (borrowed from NTL) and some use of non-normalized intermediate values.
	 *
	 * For some uses this is the most efficient field for primes in the range from half word
	 * to 2^30.
	 *
	 * Requires: Modulus < 2^30.
	 * Intended use: 2^15 < prime modulus < 2^30.
	\ingroup field
	 */
	template <>
		class Modular<int> : public FieldInterface {

		protected:

		int modulus;

		double modulusinv;
		
		int _two64;

		public:	       

		friend class FieldAXPY<Modular<int> >;
                friend class DotProductDomain<Modular<int> >;
		friend class MVProductDomain<Modular<int> >;
	       
		typedef int Element;
		typedef ModularRandIter<int> RandIter;
                
                public:	       
                    const Element zero,one;
                    Element mone;

		//default modular field,taking 65521 as default modulus
		Modular () :modulus(65521),zero(0),one(1),mone(65520) {
			modulusinv=1/(double)65521;
			
			_two64 = (int) ((uint64) (-1) % (uint64) 65521);
			_two64 += 1;
			if (_two64 >= 65521) _two64 -= 65521;
		}

		Modular (int value)  : modulus(value),zero(0),one(1),mone(value-1) {
			modulusinv = 1 / ((double) value); 
			if(value<=1) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
			if(value>LINBOX_MAX_MODULUS) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
			_two64 = (int) ((uint64) (-1) % (uint64) value);
			_two64 += 1;
			if (_two64 >= value) _two64 -= value;
		}

		Modular(const Modular<int>& mf) : modulus(mf.modulus),modulusinv(mf.modulusinv),_two64(mf._two64),zero(0),one(1),mone(mf.mone){}

		Modular &operator=(const Modular<int> &F) {
			modulus = F.modulus;
			modulusinv = F.modulusinv;
			_two64 = F._two64;
                        mone = F.mone;
			return *this;
		}

	
		inline integer &cardinality (integer &c) const{ 
			return c = modulus;
		}

		inline integer &characteristic (integer &c) const {
			return c = modulus; 
		}

		inline integer &convert (integer &x, const Element &y) const { 
			return x = y;
		}

		inline double & convert (double &x, const Element &y) const { 
			return x = (double) y;
		}

		inline float & convert (float &x, const Element &y) const { 
			return x = (float) y;
		}
		
		inline std::ostream &write (std::ostream &os) const {
			return os << "int mod " << modulus;
		}
		
		inline std::istream &read (std::istream &is) {
			is >> modulus; 
			modulusinv = 1 /((double) modulus );
                        if(modulus <= 1) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be > 1");
                        if(modulus > LINBOX_MAX_MODULUS) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus is too big");
			_two64 = (int) ((uint64) (-1) % (uint64) modulus);
			_two64 += 1;
			if (_two64 >= modulus) _two64 -= modulus;
			
			return is;
		}
		
		inline std::ostream &write (std::ostream &os, const Element &x) const {
			return os << x;
		}

		inline std::istream &read (std::istream &is, Element &x) const {
			integer tmp;
			is >> tmp;
			init(x,tmp); 
			return is;
                }
		

		template<class Element1>
		inline Element &init (Element & x, const Element1 &y) const {
			x = y % modulus;
			if (x < 0) x += modulus;
			return x;
		}

		inline Element &init (Element &x, const double &y) const  { 
		  double z = fmod(y, (double)modulus);
		  if (z < 0) z += (double)modulus;
		  //z += 0.5; //C Pernet: sounds nasty and not necessary
		  return x = static_cast<long>(z); //rounds towards 0
		}

		inline Element &init (Element &x, const float &y) const  { 
		  return init(x , (double) y);
		}

		inline Element &init (Element &x, const integer &y) const  {
			x = y % modulus;
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
			if ( x >= modulus ) x -= modulus;
			return x;
		}
 
		inline Element &sub (Element &x, const Element &y, const Element &z) const {
			x = y - z;
			if (x < 0) x += modulus;
			return x;
		}
		
		inline Element &mul (Element &x, const Element &y, const Element &z) const {
			int q;

			q  = (int) ((((double) y)*((double) z)) * modulusinv);  // q could be off by (+/-) 1
			x = (int) (y*z - q*modulus);
			
			
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
			if(y == 0) return x=0;
			else return x = modulus-y;
		}
 
		inline Element &inv (Element &x, const Element &y) const {
			int d, t;			
			XGCD(d, x, t, y, modulus);
			if (d != 1) {
				throw PreconditionFailed(__FUNCTION__,__LINE__,"InvMod: Input is not invertible ");
			}
			if (x < 0)
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
			
			
			if (r >= modulus)
				r -= modulus;
			else if (r < 0)
				r += modulus;

			return r;

		}

		inline Element &addin (Element &x, const Element &y) const {
			x += y;
			if (  x >= modulus ) x -= modulus;
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
			int q;
			
			q  = (int) (((((double) a) * ((double) x)) + (double) r) * modulusinv);  // q could be off by (+/-) 1
			r = (int) (a * x + r - q*modulus);
			
			
			if (r >= modulus)
				r -= modulus;
			else if (r < 0)
				r += modulus;

			return r;
		}

		private:

      		static void XGCD(int& d, int& s, int& t, int a, int b) {
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
		class FieldAXPY<Modular<int> > {	  
		public:
	  
		typedef int Element;
		typedef Modular<int> Field;
	  
		FieldAXPY (const Field &F) : _F (F),_y(0) {}
		

		FieldAXPY (const FieldAXPY &faxpy) : _F (faxpy._F), _y (0) {}
	  
		FieldAXPY<Modular<int> > &operator = (const FieldAXPY &faxpy) {
			_F = faxpy._F; 
			_y = faxpy._y; 		       
			return *this; 
		}
	  
		inline uint64& mulacc (const Element &a, const Element &x) {
			uint64 t = (uint64) a * (uint64) x;
			_y += t;
			if (_y < t)
				return _y += _F._two64;
                        else
                            return _y;
		}

		inline uint64& accumulate (const Element &t) {
			_y += t;
			if (_y < (uint64)t)
				return _y += _F._two64;
                        else
                            return _y;
		}

		inline Element& get (Element &y) {
			y =_y % (uint64) _F.modulus;
			return y;
		}

		inline FieldAXPY &assign (const Element y) {
			_y = y; 
			return *this;
		}

		inline void reset() {
			_y = 0;
		}

	  	protected:
		Field _F;
		uint64 _y;		
	};


	template <>
		class DotProductDomain<Modular<int> > : private virtual VectorDomainBase<Modular<int> > {	       

		public:	  
		typedef int Element;	  
		DotProductDomain (const Modular<int> &F)
			: VectorDomainBase<Modular<int> > (F) {}
	  
	  
		protected:
		template <class Vector1, class Vector2>
			inline Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const {
		  
			typename Vector1::const_iterator i;
			typename Vector2::const_iterator j;
		  
			uint64 y = 0;
			uint64 t;
		  
			for (i = v1.begin (), j = v2.begin (); i < v1.end (); ++i, ++j) {
				t = ( (uint64) *i ) * ( (uint64) *j );
				y += t;
			  
				if (y < t)
					y += _F._two64;
			}
		  
			y %= (uint64) _F.modulus; 
			return res = y;

		}
	  
		template <class Vector1, class Vector2>
			inline Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const {		  
			typename Vector1::first_type::const_iterator i_idx;
			typename Vector1::second_type::const_iterator i_elt;
		  
			uint64 y = 0;
			uint64 t;
		  
			for (i_idx = v1.first.begin (), i_elt = v1.second.begin (); i_idx != v1.first.end (); ++i_idx, ++i_elt) {
				t = ( (uint64) *i_elt ) * ( (uint64) v2[*i_idx] );
				y += t;
			  
				if (y < t)
					y += _F._two64;
			}
		  

			y %= (uint64) _F.modulus;
		  
			return res = y;
		}
	};
	  
	// Specialization of MVProductDomain for int32 modular field	

	template <>
		class MVProductDomain<Modular<int32> >
		{
		public:

			typedef int32 Element;

		protected:
			template <class Vector1, class Matrix, class Vector2>
				inline Vector1 &mulColDense
				(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v) const
				{
					return mulColDenseSpecialized
						(VD, w, A, v, VectorTraits<typename Matrix::Column>::VectorCategory ());
				}

		private:
			template <class Vector1, class Matrix, class Vector2>
				Vector1 &mulColDenseSpecialized
				(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
				 VectorCategories::DenseVectorTag) const;
			template <class Vector1, class Matrix, class Vector2>
				Vector1 &mulColDenseSpecialized
				(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
				 VectorCategories::SparseSequenceVectorTag) const;
			template <class Vector1, class Matrix, class Vector2>
				Vector1 &mulColDenseSpecialized
				(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
				 VectorCategories::SparseAssociativeVectorTag) const;
			template <class Vector1, class Matrix, class Vector2>
				Vector1 &mulColDenseSpecialized
				(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
				 VectorCategories::SparseParallelVectorTag) const;

			mutable std::vector<uint64> _tmp;
		};

	template <class Vector1, class Matrix, class Vector2>
		Vector1 &MVProductDomain<Modular<int32> >::mulColDenseSpecialized
		(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::DenseVectorTag) const {
		
		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());
		
		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<uint64>::iterator l;

		uint64 t;

		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());
		
		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);
		
		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint64) *k) * ((uint64) *j);

				*l += t;
				
				if (*l < t)
					*l += VD.field ()._two64;
			}
		}
		
		typename Vector1::iterator w_j;
		
		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;
		
		return w;
	}
	
	template <class Vector1, class Matrix, class Vector2>
		Vector1 &MVProductDomain<Modular<int32> >::mulColDenseSpecialized
		(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseSequenceVectorTag) const
		{
			linbox_check (A.coldim () == v.size ());
			linbox_check (A.rowdim () == w.size ());
			
			typename Matrix::ConstColIterator i = A.colBegin ();
			typename Vector2::const_iterator j;
			typename Matrix::Column::const_iterator k;
			std::vector<uint64>::iterator l;
			
			uint64 t;
			
			if (_tmp.size () < w.size ())
				_tmp.resize (w.size ());
			
			std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);
			
			for (j = v.begin (); j != v.end (); ++j, ++i) {
				for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
					t = ((uint64) k->second) * ((uint64) *j);

					_tmp[k->first] += t;
					
					if (_tmp[k->first] < t)
						_tmp[k->first] += VD.field ()._two64;
				}
			}
			
			typename Vector1::iterator w_j;
			
			for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
				*w_j = *l % VD.field ().modulus;
			
			return w;
		}
	
	template <class Vector1, class Matrix, class Vector2>
		Vector1 &MVProductDomain<Modular<int32> >::mulColDenseSpecialized
		(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseAssociativeVectorTag) const {

		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());
		
		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::const_iterator k;
		std::vector<uint64>::iterator l;
		
		uint64 t;
		
		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());
		
		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);
		
		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k = i->begin (), l = _tmp.begin (); k != i->end (); ++k, ++l) {
				t = ((uint64) k->second) * ((uint64) *j);
				
				_tmp[k->first] += t;
				
				if (_tmp[k->first] < t)
					_tmp[k->first] += VD.field ()._two64;
			}
		}
		
		typename Vector1::iterator w_j;
		
		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;
		
		return w;
	}

	template <class Vector1, class Matrix, class Vector2>
		Vector1 &MVProductDomain<Modular<int32> >::mulColDenseSpecialized
		(const VectorDomain<Modular<int32> > &VD, Vector1 &w, const Matrix &A, const Vector2 &v,
		 VectorCategories::SparseParallelVectorTag) const {
		
		linbox_check (A.coldim () == v.size ());
		linbox_check (A.rowdim () == w.size ());
		
		typename Matrix::ConstColIterator i = A.colBegin ();
		typename Vector2::const_iterator j;
		typename Matrix::Column::first_type::const_iterator k_idx;
		typename Matrix::Column::second_type::const_iterator k_elt;
		std::vector<uint64>::iterator l;
		
		uint64 t;
		
		if (_tmp.size () < w.size ())
			_tmp.resize (w.size ());
		
		std::fill (_tmp.begin (), _tmp.begin () + w.size (), 0);
		
		for (j = v.begin (); j != v.end (); ++j, ++i) {
			for (k_idx = i->first.begin (), k_elt = i->second.begin (), l = _tmp.begin ();
			     k_idx != i->first.end ();
			     ++k_idx, ++k_elt, ++l)
				{
					t = ((uint64) *k_elt) * ((uint64) *j);

					_tmp[*k_idx] += t;

					if (_tmp[*k_idx] < t)
						_tmp[*k_idx] += VD.field ()._two64;
				}
		}

		typename Vector1::iterator w_j;

		for (w_j = w.begin (), l = _tmp.begin (); w_j != w.end (); ++w_j, ++l)
			*w_j = *l % VD.field ().modulus;

		return w;
	}
  	  

} 

#include "linbox/randiter/modular.h"
#endif //__LINBOX_modular_int_H

#endif 

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
