/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* linbox/field/ntl-pid-zz_p.h
 * written by bds
 *
 */

#ifndef __NTL_PID_zz_p_H
#define __NTL_PID_zz_p_H

#include "linbox/field/ntl-zz_p.h"
#include "linbox/util/debug.h"
#include "linbox-config.h"
#include <NTL/ZZ.h>

// Namespace in which all LinBox library code resides
namespace LinBox
{
  
    /** @memo extend Wrapper of zz_p from NTL.  Add PID functions
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
    };
	    
} // namespace LinBox

#endif // __NTL_PID_zz_p_H
