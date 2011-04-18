/* linbox/field/givaro-rational.h
 * Copyright (C) 2004 Gilles Villard
 *
 * Created  Gilles Villard < Gilles.Villard@ens-lyon.fr> 
 * ven oct  8 13:30:05 CEST 2004
 * see COPYING for licence information
 */


#ifndef __LINBOX_givaro_rational_H
#define __LINBOX_givaro_rational_H


#include "linbox/integer.h"
#include "linbox/field/unparametric.h"
#include <linbox/field/field-traits.h>


//------------------------------------
// Files of Givaro library

#include "givaro/givrational.h"
//------------------------------------


namespace LinBox 
{ 

	template <class Ring>
    struct ClassifyRing;

	class GivaroRational; 
	
	template<> 
	struct ClassifyRing<GivaroRational> {
		typedef RingCategories::RationalTag categoryTag;
	}; 

 class GivaroRational : public UnparametricField<Rational>
  {
 
  public:

    /** Element type.
     *  
     */
    typedef  Rational Element;
    

    Element& init(Element& x , const integer& y) const
	  { return x=Rational(y);}

    Element& init(Element& x , const integer& n,const integer& d) const
	  { return x=Rational(n,d);}

    template<class XX>
    Element& init(Element& x , const XX& y) const
	  { return x=Rational(y);}

    integer& convert(integer& i, const Element& r) const 
          { return i=r.nume(); }
      

    Element& assign(Element& x , const Rational& y) const
          { return x=y;}

    // x = numerator of y
    integer& get_num (integer& x, const Element& y)  const
	  { return x = y.nume(); }

    // x = denominator of y
    integer& get_den (integer& x, const Element& y) const 
          { return x = y.deno(); }


  }; // class 
 


} // namespace LinBox


// Specialization of homomorphism for basefield
#include "linbox/field/hom.h"
namespace LinBox 
{
	template <class _Target>
	class Hom<GivaroRational, _Target> {

	public:
		typedef GivaroRational Source;
		typedef _Target Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;
	
            Hom(const Source& S, const Target& T) : _source (S), _target(T){ }
            Elt& image(Elt& t, const SrcElt& s) {
                _source. get_num (num, s);
                _source. get_den (den, s);
                if (den == 1) {
                    return _target.init(t,num);
                } else if (num == 1) {
                    _target.init(t,den);
                    return _target.invin(t);
                } else {
                    _target. init (tmp, den);
                    _target. init (t, num);
                    return _target. divin (t, tmp);
                }
// 			_target. init (t, den);
// 			return _target. invin (t);
            }
            SrcElt& preimage(SrcElt& s, const Elt& t) {
                _target. convert (num, t);
                _source. init (s, num);
                return s;
            }
            const Source& source() { return _source;}
            const Target& target() { return _target;}
            
	protected:
		integer num, den;
		Elt tmp;
		Source _source;
		Target _target;
	}; // end Hom 

	template <>
	class Hom<GivaroRational,GivaroRational> {

	public:
		typedef GivaroRational Source;
		typedef Source Target;
		typedef Source::Element SrcElt;
		typedef Target::Element Elt;
	
		Hom(const Source& S, const Target& T) : _source (S), _target(T){}
		Elt& image(Elt& t, const SrcElt& s) {
			_target.assign(t, s);
			return t;
		}
		SrcElt& preimage(SrcElt& s, const Elt& t) {
			_source.assign(s, t);
			return s;
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		Source _source;
		Target _target;
	}; // end Hom 
}

#endif // __LINBOX_givaro_rational_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
