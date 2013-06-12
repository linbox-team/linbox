/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/image-field.h
 * Written by David Saunders
 * See COPYING for license information.
 */

#ifndef __IMAGE_FIELD_H
#define __IMAGE_FIELD_H

#include "linbox/field/hom.h"
namespace LinBox
{

/**
 * ImageFields are fields which are targets of a ring homomorphism from a source ring.
 * Usually the source ring is ZZ or another field.
 * The ImageField can be viewed as a module over the source ring.  Hence we could think
 * of the source as the scalar domain.  Thus a mul of a source element times a
 * image element can be thought of as a scalar multiplication and we call the op smul().
 * Often this can be done faster than taking the image of the scalar and multiplying
 * in the image field.
 *
 * ImageFields provide smul() and saxpy() as well as the image() and preimage() maps of the
 * homomorphism.
 *
 * The generic default is to build the ImageField functions directly from a Hom object,
 * but the intention is to provide efficient specializations for specific cases.
 *
 * Issues abound.  How do we convey that the scalars are small so a implementation
 * of a modular image can be efficient in that case?
 */
template<class Source, class Target>
class ImageField: public Target
{   protected:
	typedef typename Target::Element Elt;
	typedef typename Source::Element SrcElt;
	Hom<Source, Target> _nat;
    public:
	ImageField(){}
	ImageField(Hom<Source, Target> nat):Target(nat.target()),_nat(nat) {}

	ImageField(Source& S, Target& T) : Target(T),_nat(Hom<Source,Target>(S, T)) {}

	Elt & image(Elt& t, const SrcElt& s)
	{ return _nat.image(t, s); }

	SrcElt& preimage(SrcElt& s, const Elt& t)
	{ return _nat.preimage(s, t); }

	/// y <-- a*x
	Elt& smul(Elt& y, const SrcElt& a, const Elt& x)
	{ Elt w;mul(y,_nat.image(w,a),x);}

	/// x <-- a*x
	Elt& smulin(Elt& x, const SrcElt& a)
	{ Elt w;   return this->mulin(x, _nat.image(w, a)); }

	/// z <-- a*x + y
	Elt& saxpy(Elt& z, const SrcElt& a, const Elt& x, const Elt& y)
	{ Elt w;   return this->axpy(z, _nat.image(w, a), x, y); }

	/// y <-- a*x + y
	Elt& saxpyin(Elt& y, const SrcElt& a, const Elt& x)
	{ Elt w;   return this->axpyin(y, _nat.image(w, a), x); }

	/*
	/// dot product of source vector and image vector
	Elt& dot(vector of src, vector of target)
	{ allocate vector of target; map vector of src to it; target dot }

	/// accumulator of source * target product (i.e accum results of smul's).
	accumulator...
	*/
}; // ImageField<Source, Target>


// Trial specialization.  Try one member. one could similarly do other member functions.
// Delete this when better examples are available later.
//
template<> uint32_t& ImageField<Modular<uint16_t>, Modular<uint32_t> >::
smul(uint32_t& y, const uint16_t& a, const uint32_t & x)
	{ return mul(y, static_cast<uint32_t>(a), x); }

}// namespace LinBox
#endif // __IMAGE_FIELD_H
