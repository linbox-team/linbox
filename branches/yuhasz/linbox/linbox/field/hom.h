/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/hom.h
 * Written by David Saunders
 * See COPYING for license information.
 */

#ifndef __HOM_H
#define __HOM_H

namespace LinBox {
/** 
 * An instance of Hom is a homomorphism from a ring of type Source
 * to a ring (usually field) of type Target.  The intended use is that
 * it will be a natural mapping.  For instance: 
 * 
 * Hom<Integers, Modular<integer> >(Z, Zp) nat; // is the mod p mapping.
 *
 * Hom<Unparametric<NTL_ZZp, Modular<integer> >(Zp, Mp) nat; 
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

	Hom(){}
	/**
	 * Construct a homomorphism from a specific source ring S and target 
	 * field T with Hom(S, T).  The default behaviour is error.  
	 * Specializations define all actual homomorphisms.
	 */
	Hom(Source& S, Target& T){ throw NoHomError(); }

	/** 
	 * image(t, s) implements the homomorphism, assigning the 
	 * t the value of the image of s under the mapping.
	 *
	 * The default behaviour is a no-op.
	 */
	Elt& image(Elt& t, const SrcElt& s) {}

	/** If possible, preimage(s,t) assigns a value to s such that 
	 * the image of s is t.  Otherwise behaviour is unspecified.
	 * An error may be thrown, a conventional value may be set, or
	 * an arb value set.
	 *
	 * The default behaviour is a no-op.
	 */
	SrcElt& preimage(SrcElt& s, const Elt& t) {}
	const Source& source() { return _source;}
	const Target& target() { return _target;}

	private:
	Source _source;
	Target _target;
};

/// Error object for attempt to establish a Hom that cannot exist.
class NoHomError{};

#include "linbox/field/modular.h"
/// Specialization to Modular<uint16> --> Modular<uint_32>.
// Just a trial.  delete this when better examples exist.
template<> Hom<Modular<uint16>, Modular<uint32> >::
	Hom(Modular<uint16>& S, Modular<uint32>& T ): _source(S),_target(T)
	{
			integer ps, pt;
			if (S.characteristic(ps) != T.characteristic(pt)) throw NoHomError();
	}

template<> Modular<uint32>::Element& Hom<Modular<uint16>, Modular<uint32> >::
image(Modular<uint32>::Element& t, const Modular<uint16>::Element& s) { return t = s; }

// assumes t normalized.
template<> Modular<uint16>::Element& Hom<Modular<uint16>, Modular<uint32> >::
preimage(Modular<uint16>::Element& s, const Modular<uint32>::Element& t) { return s = t; }

}// namespace LinBox
#endif // __HOM_H
