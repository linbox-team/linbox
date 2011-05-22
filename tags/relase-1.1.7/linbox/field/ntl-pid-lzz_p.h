/* Copyright (C) 2010 LinBox
 * written by bds
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


#ifndef __LINBOX_ntl_pid_zz_p_H
#define __LINBOX_ntl_pid_zz_p_H

#include "linbox/field/ntl-lzz_p.h"
#include "linbox/util/debug.h"
#include "linbox/linbox-config.h"
#include <NTL/ZZ.h>
#include <linbox/field/field-traits.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{
	template <class Ring>
	struct ClassifyRing;

	class NTL_PID_zz_p;

	template<>
	struct ClassifyRing<NTL_PID_zz_p> {
		typedef RingCategories::ModularTag categoryTag;
	};

    /** \brief extend Wrapper of zz_p from NTL.  Add PID functions
	\ingroup field
     */
    struct NTL_PID_zz_p: public NTL_zz_p
    {
	protected: long _modulus;
	public:
	NTL_PID_zz_p(long pp, int exp = 1) 
    	: NTL_zz_p(pp), _modulus(pp) {
		if( exp != 1 ) throw PreconditionFailed(__FUNCTION__,__LINE__,"exponent must be 1");
	}

        Element& gcd(Element& g, const Element& a, const Element& b) const
	{   g = NTL::GCD(NTL::rep(a), NTL::rep(b));  
	    g = NTL::GCD(NTL::rep(g), _modulus);
	    return g;
	}

        Element& gcdin(Element& a, const Element& b) const
	{   return gcd(a,a, b);  }

	bool isUnit(const Element& a) const
	{   return 1 == NTL::GCD(NTL::rep(a), _modulus);  }

	Element& div(Element& c, const Element& a, const Element& b) const
	{   return c = NTL::rep(a)/NTL::GCD(NTL::rep(a),NTL::rep(b));   }
	Element& divin(Element& a, const Element& b) const
	{   return div(a, a, b);   }

	static inline double getMaxModulus() { return (double)NTL_SP_BOUND; }
    };
	    
} // namespace LinBox

#endif // __LINBOX_ntl_pid_zz_p_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
