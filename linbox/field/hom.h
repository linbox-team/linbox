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
#include "linbox/util/error.h"
#include "linbox/integer.h"
#include <givaro/givcaster.h>
#include <givaro/qfield.h>
#include <givaro/modular.h>
#include <givaro/zring.h>

#ifdef __LINBOX_HAVE_NTL
#include "linbox/ring/ntl/ntl-zz.h"
#endif //__LINBOX_HAVE_NTL

#include "givaro/givrational.h"

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
	 * Hom<Unparametric<Integers>, Givaro::Modular<integer> >(Z, Zp) nat; // is the mod p mapping.
	 *
	 * Hom<<NTL_ZZp, Givaro::Modular<integer> >(Zp, Mp) nat;
	 *
	 * // is an isomorphism, provided the Zp and Mp have same characteristic.
	 * Hom<Unparametric<NTL_ZZp, Unparameteric<NTL_ZZpEx> >(Z3, Z27) nat;
	 * // is the embedding of the prime field in the extension.
	 */

	template< class Source, class Target, class Enabled = void > // Enabled is just for std::enable_if 
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
		const Source& _source;
		const Target& _target;
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
			_source.assign (t, s);
			return t;
		}

		SrcElt& preimage(SrcElt& s, const Elt& t)
		{
			_source.assign (s, t);
			return s;
		}

		const Source& source() { return _source;}
		const Target& target() { return _source;}

	protected:
		const Source& _source;
	};


}

// // Specialization to modular
// namespace LinBox
// {
//     // GENERIC SPECIALISATION FOR MODULAR HOMORMOPHISM
//     // The enable_if is to protect undetermination with Hom<Source, Source>
// 	template<typename T1, typename C1, typename T2, typename C2>
// 	class Hom<Givaro::Modular<T1, C1>, Givaro::Modular<T2, C2>,
//         typename std::enable_if<!std::is_same<Givaro::Modular<T1, C1>, Givaro::Modular<T2, C2>>::value>::type>
//     {
// 	public:
// 	    using Element1 = typename Givaro::Modular<T1, C1>::Element;
// 	    using Element2 = typename Givaro::Modular<T2, C2>::Element;
	
// 		Hom(const Givaro::Modular<T1, C1>& S, const Givaro::Modular<T2, C2>& T)
// 		{
// 			integer ps, pt;
// 			if (S.characteristic(ps) != T.characteristic(pt)) throw NoHomError();
// 		}
		
// 		inline Element2& image(Element2& t, const Element1& s) const
// 		{
// 			return Givaro::Caster<Element2,Element1>(t,s);
// 		}
		
// 		// assumes t normalized.
// 		inline Element1& preimage(Element1& s, const Element2& t) const
// 		{
// 			return Givaro::Caster<Element1,Element2>(s,t);
// 		}
// 	};
// }// namespace LinBox

namespace LinBox
{

    template<class _Source>
	class Hom<_Source, Givaro::ZRing<Integer> > {

	public:
        typedef Givaro::ZRing<Integer> Target;
		typedef _Source Source;
		typedef typename Source::Element SrcElt;
		typedef Integer Elt;

		Hom(const Source& S, const Target& T) :
			_source (S), _target (T)
		{}
		inline Integer& image(Integer& t, const SrcElt& s)
		{
			return Givaro::Caster(t,s);
		}
		inline SrcElt& preimage(SrcElt& s, const Integer& t)
		{
                        return Givaro::Caster(s,t);
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		const Source& _source;
		const Target& _target;
	}; // end Hom

        
	template<class _Target>
	class Hom<Givaro::ZRing<Integer>, _Target> {

	public:
		typedef Givaro::ZRing<Integer> Source;
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
		const Source& _source;
		const Target& _target;
	}; // end Hom

	template<>
	class Hom<Givaro::ZRing<Integer>, Givaro::ZRing<Integer> > {

	public:
		typedef Givaro::ZRing<Integer> Source;
		typedef Source Target;
		typedef Integer SrcElt;
		typedef Integer Elt;

		Hom(const Source& S, const Target& T) :
			_source (S), _target (T)
		{}
		inline Elt& image(Elt& t, const SrcElt& s)
        {   
			return t=s;
		}
		inline SrcElt& preimage(SrcElt& s, const Elt& t)
		{
			return s=t;
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		const Source& _source;
		const Target& _target;
	}; // end Hom

	template <class _Target>
	class Hom<Givaro::QField<Givaro::Rational>, _Target> {

	public:
		typedef Givaro::QField<Givaro::Rational> Source;
		typedef _Target Target;
		typedef typename Source::Element SrcElt;
		typedef typename Target::Element Elt;

		Hom(const Source& S, const Target& T) :
			_source (S), _target(T)
		{ }
		Elt& image(Elt& t, const SrcElt& s) {
			if (_source.isOne(s.deno())) {
				return _target.init(t,s.nume());
			}
			else if (_source.isOne(s.nume())) {
				_target.init(t,s.deno());
				return _target.invin(t);
			}
			else {
				_target. init (tmp, s.deno());
				_target. init (t, s.nume());
				return _target. divin (t, tmp);
			}
		}
            // Warning, by default convert only the numerator
		SrcElt& preimage(SrcElt& s, const Elt& t) {
			_target. convert (s.nume(), t);
			_source. init (s, s.nume());
			return s;
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	protected:
		Elt tmp;
		Source _source;
		Target _target;
	}; // end Hom


	template <>
	class Hom<Givaro::QField<Givaro::Rational>, Givaro::QField<Givaro::Rational>> {

	public:
		typedef Givaro::QField<Givaro::Rational> Domain;
		typedef Domain::Element Elt;

		Hom(const Domain& S, const Domain& T) : _d(S) { }
		Elt& image(Elt& t, const Elt& s) const { return t=s; }
            // Warning, by default convert only the numerator
		Elt& preimage(Elt& s, const Elt& t) const { return s=t; }
		const Domain& source() const { return _d;}
		const Domain& target() const { return _d;}

	protected:
		const Domain& _d;
	}; // end Hom


} // namespace LinBox


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
		Integer tmp;
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

#include <givaro/extension.h>

namespace LinBox
{
	template<>
	class Hom < GF2, Givaro::Extension<GF2> > {
		typedef GF2 Source;
		typedef Givaro::Extension<GF2> Target;
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
	class Hom <Givaro::GFq,Givaro::GFq> {
	public:
		typedef Givaro::GFq Source;
		typedef Givaro::GFq Target;

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
		 * The default behaviour goes through Integers.
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
		 * The default behaviour goes through Integers.
		 */
		SrcElt& preimage(SrcElt& s, const Elt& t)
		{
			return _source.init(s, _target.convert(tmp,t));
		}
		const Source& source() { return _source;}
		const Target& target() { return _target;}

	private:
		Integer tmp;
		Source _source;
		Target _target;
	}; // end Hom
}
#endif



#endif // __LINBOX_hom_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
