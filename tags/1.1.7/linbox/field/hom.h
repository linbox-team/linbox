/* linbox/field/hom.h
 * Copyright(C) LinBox
 * Written by David Saunders
 * See COPYING for license information.
 */

#ifndef __LINBOX_hom_H
#define __LINBOX_hom_H

#include "linbox/linbox-config.h"
#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include <linbox/util/error.h>

#ifdef __LINBOX_HAVE_NTL
#include <linbox/field/ntl-ZZ.h>
#endif //__LINBOX_HAVE_NTL

namespace LinBox 
{
	/// Error object for attempt to establish a Hom that cannot exist.
	class NoHomError {};

	/**  \brief map element of source ring(field) to target ring
	\ingroup field

	 * An instance of Hom is a homomorphism from a ring of type Source
	 * to a ring (usually field) of type Target.  The intended use is that
	 * it will be a natural mapping.  For instance: 
	 * 
	 * Hom<Unparametric<Integers>, Modular<integer> >(Z, Zp) nat; // is the mod p mapping.
	 *
	 * Hom<<NTL_ZZp, Modular<integer> >(Zp, Mp) nat; 
	 *
	 * // is an isomorphism, provided the Zp and Mp have same characteristic.
	 * Hom<Unparametric<NTL_ZZp, Unparameteric<NTL_ZZpEx> >(Z3, Z27) nat; 
	 * // is the embedding of the prime field in the extension.
	 */

	template< class Source, class Target > 
	class Hom 
	{   public:
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;

		//Hom(){}
		/**
		 * Construct a homomorphism from a specific source ring S and target 
		 * field T with Hom(S, T).  
		 * Specializations define all actual homomorphisms.
		 */
		Hom(const Source& S, const Target& T) : _source(S), _target(T){ }

		/** 
		 * image(t, s) implements the homomorphism, assigning the 
		 * t the value of the image of s under the mapping.
		 *
		 * The default behaviour goes through integers.
		 */
		Elt& image(Elt& t, const SrcElt& s) {
                    return _target.init(t, _source.convert(tmp,s));
                }

		/** If possible, preimage(s,t) assigns a value to s such that 
		 * the image of s is t.  Otherwise behaviour is unspecified.
		 * An error may be thrown, a conventional value may be set, or
		 * an arb value set.
		 *
		 * The default behaviour goes through integers.
		 */
		SrcElt& preimage(SrcElt& s, const Elt& t) {
                    return _source.init(s, _target.convert(tmp,t));
                }
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	private:
            integer tmp;
		Source _source;
		Target _target;
	}; // end Hom 



	// specialization for equal domain TYPES
	// WARNING this FORBIDS same type homomorphism
	template <class Source>
	class Hom<Source, Source> {

	public:
		typedef Source Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;
	
		Hom(const Source& S, const Target& T) : _source (S){}
		Elt& image(Elt& t, const SrcElt& s) {
			_source. assign (t, s);
			return t;
		}
		SrcElt& preimage(SrcElt& s, const Elt& t) {
			_source. assign (s, t);
			return s;
		}
		const Source& source() { return _source;}
		const Target& target() { return _source;}

	protected:
		Source _source;
	}; // end Hom 


}

#ifdef __LINBOX_field_modular_H
// including specialization to modular
//#include "linbox/field/modular.h"
/// Specialization to Modular<uint16> --> Modular<uint_32>.
// Just a trial.  delete this when better examples exist.
namespace LinBox 
{
	template<> inline Hom<Modular<uint16>, Modular<uint32> >::
	Hom(const Modular<uint16>& S, const Modular<uint32>& T ): _source(S),_target(T)
	{
		integer ps, pt;
		if (S.characteristic(ps) != T.characteristic(pt)) throw NoHomError();
	}

	template<> inline Modular<uint32>::Element& Hom<Modular<uint16>, Modular<uint32> >::
	image(Modular<uint32>::Element& t, const Modular<uint16>::Element& s) { return t = s; }

	// assumes t normalized.
	template<> inline Modular<uint16>::Element& Hom<Modular<uint16>, Modular<uint32> >::
	preimage(Modular<uint16>::Element& s, const Modular<uint32>::Element& t) { return s = t; }

}// namespace LinBox
#endif //__LINBOX_field_modular_H

#ifdef __LINBOX_field_unparametric_H
namespace LinBox
{

	template<class _Target>
	class Hom<UnparametricField<integer>, _Target> {

	public:
		typedef UnparametricField<integer> Source;
		typedef _Target Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;
	
		Hom(const Source& S, const Target& T) : _source (S), _target (T) {}
		inline Elt& image(Elt& t, const SrcElt& s) {
			_target. init (t, s);
			return t;
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t) {
			_target. convert (s, t);
			return s;
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		Source _source;
		Target _target;
	}; // end Hom 

	template<>
	class Hom<UnparametricField<integer>, UnparametricField<integer> > {

	public:
		typedef UnparametricField<integer> Source;
		typedef UnparametricField<integer> Target;
		typedef Source::Element SrcElt;
		typedef Target::Element Elt;
	
		Hom(const Source& S, const Target& T) : _source (S), _target (T) {}
		inline Elt& image(Elt& t, const SrcElt& s) {
			t = s;
			return t;
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t) {
			s = t;
			return s;
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		Source _source;
		Target _target;
	}; // end Hom 

	template<class _Target>
	class Hom<PID_integer, _Target> {

	public:
		typedef PID_integer Source;
		typedef _Target Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;
	
		Hom(const Source& S, const Target& T) : _source (S), _target (T) {}
		inline Elt& image(Elt& t, const SrcElt& s) {
			_target. init (t, s);
			return t;
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t) {
			_target. convert (s, t);
			return s;
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		Source _source;
		Target _target;
	}; // end Hom 


	template<>
	class Hom<PID_integer, PID_integer> {

	public:
		typedef PID_integer Source;
		typedef Source Target;
		typedef Source::Element SrcElt;
		typedef Target::Element Elt;
	
		Hom(const Source& S, const Target& T) : _source (S), _target (T) {}
		inline Elt& image(Elt& t, const SrcElt& s) {
			_target. assign (t, s);
			return t;
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t) {
			_target. assign (s, t);
			return s;
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		Source _source;
		Target _target;
	}; // end Hom 


#if 0
	  #ifdef __FIELD_MODULAR_H
	  // Dan Roche mapping from UnparametricField to Modular - for integer
	  // computations that use mod one or more primes and possibly chinese
	  // remaindering.
	  template<class _INT1, class _INT2>
	  class Hom<UnparametricField<_INT1 >,Modular<_INT2 > > {

	  public:
	  typedef UnparametricField<_INT1 > Source;
	  typedef Modular<_INT2 > Target;
	  typedef _INT1 SrcElt;
	  typedef _INT2 Elt;

	  Hom(const Source& S, const Target& T) :_source(S), _target(T) {}

	  inline Elt& image(Elt& t, const SrcElt& s) {
	  integer temp;
	  return _target.init(t,_source.convert(temp,s));
	  }
	  inline SrcElt& preimage(SrcElt& s, const Elt& t) {
	  integer temp;
	  return _source.init(s,_source.convert(temp,t));
	  }
	  const Source& source() { return _source; }
	  const Target& target() { return _target; }

	  protected:
	  Source _source;
	  Target _target;
	  }; // end Hom

	  #endif // __FIELD_MODULAR_H
#endif

} // namespace LinBox
#endif //__LINBOX_field_unparametric_H

#ifdef __LINBOX_HAVE_NTL
namespace LinBox 
{

	template<class _Target > 
	class Hom <NTL_ZZ , _Target> {

	public:
		typedef NTL_ZZ Source;
		typedef _Target Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;
	
		Hom(const Source& S, const Target& T) : _source(S), _target(T){}
		inline Elt& image(Elt& t, const SrcElt& s) {
			return _target. init (t, _source. convert (tmp, s));
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t) {
			return _source. init (s, _target. convert (tmp, t) );
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		integer tmp;
		Source _source;
		Target _target;
	}; // end Hom 
	
	template<>
	class Hom <NTL_ZZ , NTL_ZZ> {
	public:
		typedef NTL_ZZ Source;
		typedef NTL_ZZ Target;
		typedef Source::Element SrcElt;
		typedef Target::Element Elt;
	
		Hom(const Source& S, const Target& T) : _source(S), _target(T){}
		inline Elt& image(Elt& t, const SrcElt& s) {
			return _target.assign(t, s);
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t) {
			return _source.assign(s, t);
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		Source _source;
		Target _target;
	}; // end Hom 
	
}
#endif //__LINBOX_HAVE_NTL

#ifdef __LINBOX_field_gmp_rational_H
namespace LinBox 
{
	template <class _Target>
	class Hom<GMPRationalField, _Target> {

	public:
		typedef GMPRationalField Source;
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
	class Hom<GMPRationalField, GMPRationalField> {

	public:
		typedef GMPRationalField Source;
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
#endif //__LINBOX_field_gmp_rational_H

#endif // __LINBOX_hom_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
