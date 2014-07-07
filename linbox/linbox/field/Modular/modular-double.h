/* linbox/field/modular-double.h
 * Copyright (C) 2003 Pascal Giorgi
 *
 * Written by Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * ------------------------------------
 *
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the  GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file field/Modular/modular-double.h
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
#include "linbox/randiter/nonzero.h"
#include "linbox/randiter/modular.h"
#include "linbox/util/write-mm.h"

#include <fflas-ffpack/field/modular-double.h>


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
		      typedef FFPACK::Modular<double> Father_t;
		      friend class FieldAXPY<Modular<Element> >;
		      friend class DotProductDomain<Modular<Element> >;
		      friend class MultiModDouble;

	      public:

		      typedef ModularRandIter<Element> RandIter;

		      static ClassifyRing<Modular<Element> >::categoryTag getCategory()
		      {
			      return ClassifyRing<Modular<Element> >::categoryTag();
		      }

		      Modular (const integer& p, int e=1) :
			      Father_t((unsigned long) p)
		      {
			      linbox_check(e==1);
#ifdef DEBUG
			      if(modulus <= 1)
				      throw PreconditionFailed(LB_FILE_LOC,"modulus must be > 1");
			      if(modulus > getMaxModulus())
				      throw PreconditionFailed(LB_FILE_LOC,"modulus is too big");
#endif
		      }

		      Modular () : Father_t() {};

		      using Father_t ::cardinality ;
		      integer &cardinality (integer &c) const
		      {
			      return c = integer(modulus);
		      }

		      using Father_t ::characteristic;
		      integer &characteristic (integer &c) const
		      {
			      return c = integer(modulus);
		      }

		      using Father_t ::convert;
		      integer &convert (integer &x, const Element &y) const
		      {
			      return x = integer(y);
		      }


		      using Father_t ::inv;
		      Element &inv (Element &x, const Integer&y) const
		      {
			      init(x,y);
			      return invin(x);
		      }

		      //!@bug use FFPACK operator
		      const Modular<double> &operator=(const Modular<double> &F)
		      {
			      if ( this == &F)
				      return *this;
			      modulus  = F.modulus;
			      lmodulus = F.lmodulus;

			      F.assign(const_cast<Element&>(one),F.one);
			      F.assign(const_cast<Element&>(zero),F.zero);
			      F.assign(const_cast<Element&>(mOne),F.mOne);
			      return *this;
		      }


		      using Father_t ::init;
		      Element &init (Element &x, const integer &y) const
		      {
			      x = (Element)(y%lmodulus);
			      if (x<0) x+= modulus ;
			      linbox_check(x < lmodulus);
			      linbox_check(!(x < 0));
			      return x  ;
		      }

		      Element &init (Element &x) const
		      {
			      return x = 0 ;
		      }

		       bool isMinusOne (const Element &x) const
		      {
			      return (x == mOne);
		      }

		      /** Max number of operations before reducing
		       * @param r if \c r=0, we consider how many \c += are performable.
		       * if \c r=1, then we look for the maximum \c axpy operations doable.
		       * @return \p 0 if the field is too big, a positive number otherwise, \p -1 if infinity
		       * on general fields, it is \p 1.
		       */
		      unsigned long AccBound(const Element r) const
		      {
			      // Element One, Zero ; init(One,1UL) ; init(Zero,0UL);
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
			      }
			      else
				      throw LinboxError("Bad input, expecting 0 or 1");
			      return 0;
		      }

		/*- Print field as a constructor call.
		 * @return output stream to which field is written.
		 * @param  os  output stream to which field is written.
		 * @param  F  optional name to give the field in the description.  IF F is the null string, only the class name is written.
		 * Example: For element type double and modulus 101,
		 * write(os) produces      "Modular< double > ( 101 )"  on os,
 		 * write(os, "F") produces "Modular< double > F( 101 )" on os, and
 		 * write(os, "") produces  "Modular< double >"          on os.
		 */
		std::ostream &write (std::ostream &os) const
		{
		  integer p = cardinality();
		  return os << "Modular<" << eltype( Element() ) << " >( " << p << " )";
		}

		std::ostream &write (std::ostream &os, std::string F) const
		{
		  os << "Modular<" << eltype( Element() ) << " > "; // class name
		  if (F != "") {
		    integer p = cardinality();
		    os << F << "( " << p << " )"; // show constuctor args
		  }
		  return os;
		}

        std::ostream &write (std::ostream & os, const Element & x) const
	{
            return Father_t::write(os,x);
        }

		Element& next(Element &x) const
		{
			return addin(x,one);
		}

    };

} // LinBox

// FieldAXPY/DotProductDomain
namespace LinBox
{

	template <>
	class FieldAXPY<Modular<double> > {
	public:

		typedef double Element;
		typedef double Abnormal;
		typedef Modular<double> Field;

		FieldAXPY (const Field &F) :
			_field (&F) , //_invmod(1./field().modulus),
			_y(0.) , _bound( (double) ((1ULL << 53) - (unsigned long int) (field().modulus*field().modulus)))
		{}

		FieldAXPY (const FieldAXPY &faxpy) :
			_field (faxpy._field),// _invmod(faxpy._invmod) ,
			_y(faxpy._y), _bound(faxpy._bound)
		{}

#if 0
		FieldAXPY<Modular<double> > &operator = (const FieldAXPY &faxpy)
		{
			_field = faxpy._field;
			//_invmod= faxpy._invmod;
			_y= faxpy._y;
			_bound= faxpy._bound;
			return *this;
		}
#endif

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
				return _y = fmod (_y, field().modulus);
			else
				return _y;
		}

		 Element& subumulate (const Element &tmp)
		{
			_y -= tmp;
			if (_y < 0)
				return _y += field().modulus;
			else
				return _y;
		}

		 Element& get (Element &y)
		{
			_y = fmod (_y, field().modulus);
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
				return _y = fmod (_y, field().modulus);
			else
				return _y;
		}

		inline const Field & field() const { return *_field; }

	protected:

		const Field *_field;
		//double _invmod;
		double _y;
		double _bound;
	};

	template <>
	class DotProductDomain<Modular<double> > : public virtual VectorDomainBase<Modular<double> > {
	private:
		// double _bound; // BB : not used
		size_t _nmax;
		//double _invmod;

	public:
		DotProductDomain () { /*std::cerr << "DPD-Md def cstor" << std::endl;*/ }
		typedef double Element;
		DotProductDomain (const Modular<double> &F) :
			VectorDomainBase<Modular<double> > (F)
			// , _bound( (double) ( (1ULL<<53) - (unsigned long int) (F.modulus*F.modulus)))
			, _nmax(0)//, _invmod(1./field().modulus)
		{
			_nmax= (size_t)floor((double(1<<26)* double(1<<26)*2.)/ (F.modulus * F.modulus));
			_nmax = (_nmax>0?_nmax:1);
		}

		using VectorDomainBase<Modular<double> >::field;

	protected:
		template <class Vector1, class Vector2>
		 Element &dotSpecializedDD (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;
			if (v1.size() < _nmax) {
				for (size_t i = 0; i< v1.size();++i)
					y += v1[i] * v2[i] ;
				y = fmod(y, field().modulus);
			}
			else{
				double t = 0.;
				size_t i=0;
				for (;i< v1.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1[j] * v2[j];
					t+=fmod(y, field().modulus);
					y=0.;
				}
				for (;i < v1.size();++i)
					y += v1[i] * v2[i];
				t+=fmod(y, field().modulus);
				y = fmod(t, field().modulus);
			}
			return res = y;
		}


		template <class Vector1, class Vector2>
		 Element &dotSpecializedDSP (Element &res, const Vector1 &v1, const Vector2 &v2) const
		{

			double y = 0.;

			if (v1.first.size() < _nmax) {
				for (size_t i=0;i<v1.first.size();++i)
					y+= v1.second[i] * v2[v1.first[i]];
				y = fmod(y, field().modulus);
			} else {
				double t = 0.;
				size_t i=0;
				for (;i< v1.first.size()- _nmax ;i=i+_nmax){
					for (size_t j=i;j<i+_nmax;++j)
						y += v1.second[j] * v2[v1.first[j]];
					t+=fmod(y, field().modulus);
					y=0.;
				}
				for (;i < v1.first.size();++i)
					y += v1.second[i] * v2[v1.first[i]];
				t+= fmod(y, field().modulus);
				y = fmod(t, field().modulus);
			}
			return res = y;
		}
	};
}


#endif //__LINBOX_modular_double_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
