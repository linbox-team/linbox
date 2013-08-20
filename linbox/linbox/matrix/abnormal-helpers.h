/* Copyright (C) 2013 LinBox
 * Written by AJS <stachnik@udel.edu>
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file   linbox/matrix/abnormal-helpers.h
 * @ingroup linbox/matrix
 * @brief
 */

#ifndef __LINBOX_ABNORMAL_HELPERS_H
#define __LINBOX_ABNORMAL_HELPERS_H

#include <stdlib.h>
#include <fstream>

namespace LinBox
{

template <class Field>
class AbnormalHelper {
public:
        typedef typename Field::Element Element;
	typedef typename Field::Element Abnormal;

	AbnormalHelper () {}

	AbnormalHelper(const Field& F) : field_(&F) {}

	void init(const Field& field) { field_=&field; }

	inline Abnormal& mulacc(Abnormal& x, const Element& y, const Element& z) const {
		field_->mulacc(x,y,z);
		return x;
	}

	inline Element normalize(Abnormal& elt) const {
		Element d;
		field_->init(d,elt);
		return d;
	}

protected:
	Field* field_;
};

template <>
class AbnormalHelper<Modular<double> > {
public:
	typedef Modular<double> Field;
        typedef double Element;
	typedef double Abnormal;

	AbnormalHelper () {}

	AbnormalHelper(const Field& field) {init(field);}

	void init(const Field& field) {
		modulus_=field.characteristic();
		unsigned long long maxDouble = 1ULL<<52;
		bound_=(double)maxDouble;
		field_=&field;
	}

	inline Abnormal& mulacc(Abnormal& x, const Element& y, const Element& z) const
        {
		Element d;
		d=y*z;//at most (2**26-1)**2
		x=x+d;//at most (2**26-1)**2+(2**52-1)<2**53
                maybeNormalize(x);
                return x;
        }

	inline Element& normalize(Abnormal& elt) const
	{
		return elt=fmod(elt,modulus_);
	}

protected:
	inline Abnormal& maybeNormalize(Abnormal& elt) const
	{
		if (elt>=bound_) {
			elt=fmod(elt,modulus_);
		}
		return elt;
	}

	double modulus_;

	double bound_;

	const Field* field_;
};

template <>
class AbnormalHelper<Modular<uint64_t> > {
public:
	typedef Modular<uint64_t> Field;
        typedef uint64_t Element;
	typedef uint64_t Abnormal;

	AbnormalHelper () {}

	AbnormalHelper(const Field& field) {init(field);}

	void init(const Field& field) {
		modulus_=field.characteristic();
		bound_=(uint64_t)(1ULL<<62);
		field_=&field;
	}

	inline Abnormal& mulacc(Abnormal& x, const Element& y, const Element& z) const
        {
		Abnormal d;
		d=y*z;
		x=x+d;
                maybeNormalize(x);
                return x;
        }

	inline Element normalize(Abnormal& elt) const
	{
		return elt=elt%modulus_;
	}

protected:
	inline Abnormal& maybeNormalize(Abnormal& elt) const
	{
		if (elt>=bound_) {
                        elt=elt%modulus_;
		}
		return elt;
	}

	Abnormal modulus_;

	Abnormal bound_;

	const Field* field_;
};

}

#endif // __LINBOX_ABNORMAL_HELPERS_H

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
