/* linbox/field/hom.h
 * Copyright(C) LinBox
 * Written by David Saunders
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

#ifndef __LINBOX_hom_H
#define __LINBOX_hom_H

#include "linbox/linbox-config.h"
#include "linbox/field/modular.h"
#include "linbox/field/PID-integer.h"
#include "linbox/util/error.h"

#ifdef __LINBOX_HAVE_NTL
#include "linbox/field/NTL/ntl-ZZ.h"
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
	class Hom {
	public:
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;

		//Hom(){}
		/**
		 * Construct a homomorphism from a specific source ring S and target
		 * field T with Hom(S, T).
		 * Specializations define all actual homomorphisms.
		 */
		Hom(const Source& S, const Target& T) :
			_source(S), _target(T)
		{ }

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

		Hom(const Source& S, const Target& T) :
			_source (S)
		{}
		Elt& image(Elt& t, const SrcElt& s)
		{
			_source. assign (t, s);
			return t;
		}
		SrcElt& preimage(SrcElt& s, const Elt& t)
		{
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
/// Specialization to Modular<uint16_t> --> Modular<uint_32>.
// Just a trial.  delete this when better examples exist.
namespace LinBox
{
	template<> inline Hom<Modular<uint16_t>, Modular<uint32_t> >::
	Hom(const Modular<uint16_t>& S, const Modular<uint32_t>& T ) :
		_source(S),_target(T)
	{
		integer ps, pt;
		if (S.characteristic(ps) != T.characteristic(pt)) throw NoHomError();
	}

	template<> inline Modular<uint32_t>::Element& Hom<Modular<uint16_t>, Modular<uint32_t> >::
	image(Modular<uint32_t>::Element& t, const Modular<uint16_t>::Element& s)
	{
	       	return t = s;
	}

	// assumes t normalized.
	template<> inline Modular<uint16_t>::Element& Hom<Modular<uint16_t>, Modular<uint32_t> >::
	preimage(Modular<uint16_t>::Element& s, const Modular<uint32_t>::Element& t)
	{
		linbox_check(t < UINT8_MAX) ;
	       	return s = (uint16_t) t;
	}

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

		Hom(const Source& S, const Target& T) :
			_source (S), _target (T)
		{}
		inline Elt& image(Elt& t, const SrcElt& s)
		{
			_target. init (t, s);
			return t;
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t)
		{
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

		Hom(const Source& S, const Target& T) :
			_source (S), _target (T)
		{}
		inline Elt& image(Elt& t, const SrcElt& s)
		{
			t = s;
			return t;
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t)
		{
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

		Hom(const Source& S, const Target& T) :
			_source (S), _target (T)
		{}
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

		Hom(const Source& S, const Target& T) :
			_source (S), _target (T)
		{}
		inline Elt& image(Elt& t, const SrcElt& s)
		{
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

		Hom(const Source& S, const Target& T) :
			_source(S), _target(T)
		{}

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

		Hom(const Source& S, const Target& T) :
			_source(S), _target(T)
		{}
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

		Hom(const Source& S, const Target& T) :
			_source(S), _target(T)
		{}
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

		Hom(const Source& S, const Target& T) :
			_source (S), _target(T)
		{ }
		Elt& image(Elt& t, const SrcElt& s) {
			_source. get_num (num, s);
			_source. get_den (den, s);
			if (den == 1) {
				return _target.init(t,num);
			}
			else if (num == 1) {
				_target.init(t,den);
				return _target.invin(t);
			}
			else {
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

		Hom(const Source& S, const Target& T) :
			_source (S), _target(T)
		{}
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

#ifdef __LINBOX_field_gf2_H
namespace LinBox
{

	template <>
	class Hom<GF2,GF2> {

	public:
		typedef GF2 Target;
		typedef GF2 Source;
		typedef Source::Element SrcElt;
		typedef Target::Element Elt;

		Hom(const Source& S, const Target& ) :
			_source (S)
		{}
		Elt& image(Elt& t, const SrcElt& s)
		{
			return _source.assign (t, s);
		}
		SrcElt& preimage(SrcElt& s, const Elt& t)
		{
			return _source.assign (s, t);
		}
		const Source& source()
		{
			return _source;
		}
		const Target& target()
		{
			return _source;
		}

	protected:
		Source _source;
	};

	template<class Target >
	class Hom<GF2, Target > {
	public:
		typedef GF2 Source;
		typedef typename GF2::Element SrcElt;
		typedef typename Target::Element Elt;

		Hom(const Source& S, const Target& T) :
			_source(S), _target(T)
		{ }
		Elt& image(Elt& t, const SrcElt& s)
		{
			return _source.convert(t,s);
		}
		SrcElt& preimage(SrcElt& s, const Elt& t)
		{
			return _target.convert(s,t);
		}
		stdBitReference preimage(stdBitReference s, const Elt& t) const
		{
			int ts;
			return s = _target.convert(ts, t);
		}

		const Source& source()
		{
			return _source;
		}
		const Target& target()
		{
			return _target;
		}

	private:
		Source _source;
		Target _target;
	}; // end Hom


}

#include "linbox/field/Givaro/givaro-extension.h"
namespace LinBox
{
	template<>
	class Hom < GF2, GivaroExtension<GF2> > {
		typedef GF2 Source;
		typedef GivaroExtension<GF2> Target;
	public:
		typedef Source::Element SrcElt;
		typedef Target::Element Elt;

		//Hom(){}
		/**
		 * Construct a homomorphism from a specific source ring S and target
		 * field T with Hom(S, T).  The default behaviour is error.
		 * Specializations define all actual homomorphisms.
		 */
		Hom(const Source& S, const Target& T) :
			_source(S), _target(T)
		{}

		/**
		 * image(t, s) implements the homomorphism, assigning the
		 * t the value of the image of s under the mapping.
		 *
		 * The default behaviour is a no-op.
		 */
		Elt& image(Elt& t, const SrcElt& s) const
		{
			return _target.assign(t,s);
		}

		/** If possible, preimage(s,t) assigns a value to s such that
		 * the image of s is t.  Otherwise behaviour is unspecified.
		 * An error may be thrown, a conventional value may be set, or
		 * an arb value set.
		 *
		 * The default behaviour is a no-op.
		 */
		SrcElt& preimage(SrcElt& s, const Elt& t) const
		{
			return _target.convert(s, t);
		}
		stdBitReference preimage(stdBitReference s, const Elt& t) const
		{
			bool ts;
			return s = _target.convert(ts, t);
		}

		const Source& source() const
		{
			return _source;
		}
		const Target& target() const
		{
			return _target;
		}

	private:
		Source _source;
		Target _target;
	}; // end Hom
}
#endif // __LINBOX_field_gf2_H

#ifdef __LINBOX_field_givaro_gfq_H
namespace LinBox
{
	template<>
	class Hom <GivaroGfq,GivaroGfq> {
	public:
		typedef GivaroGfq Source;
		typedef GivaroGfq Target;

		typedef Source::Element SrcElt;
		typedef Target::Element Elt;

		//Hom(){}
		/**
		 * Construct a homomorphism from a specific source ring S and target
		 * field T with Hom(S, T).
		 * Specializations define all actual homomorphisms.
		 */
		Hom(const Source& S, const Target& T) :
			_source(S), _target(T)
		{ }

		/**
		 * image(t, s) implements the homomorphism, assigning the
		 * t the value of the image of s under the mapping.
		 *
		 * The default behaviour goes through integers.
		 */
		Elt& image(Elt& t, const SrcElt& s)
		{
			return _target.init(t, _source.convert(tmp,s));
		}

		/** If possible, preimage(s,t) assigns a value to s such that
		 * the image of s is t.  Otherwise behaviour is unspecified.
		 * An error may be thrown, a conventional value may be set, or
		 * an arb value set.
		 *
		 * The default behaviour goes through integers.
		 */
		SrcElt& preimage(SrcElt& s, const Elt& t)
		{
			return _source.init(s, _target.convert(tmp,t));
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	private:
		integer tmp;
		Source _source;
		Target _target;
	}; // end Hom
}
#endif



#endif // __LINBOX_hom_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

