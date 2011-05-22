/* linbox/field/ntl-z_pE.h
 * Copyright (C) 2004  Pascal Giorgi
 *
 * Written by  Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * Modified by W. J. Turner <wjturner@acm.org>
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


#ifndef __LINBOX_ntl_lzz_pe_H
#define __LINBOX_ntl_lzz_pe_H


#include <linbox/field/unparametric.h>
#include <linbox/randiter/unparametric.h>
#include <NTL/lzz_pXFactoring.h>
#include <NTL/lzz_pE.h>
#include <time.h>
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>


namespace LinBox 
{


	template <class Ring>
	struct ClassifyRing;

	template<>
	struct ClassifyRing<UnparametricRandIter<NTL::zz_pE> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	template<>
	class UnparametricRandIter<NTL::zz_pE>
	{
	public:
		typedef NTL::zz_pE Element;
		UnparametricRandIter<NTL::zz_pE>(const UnparametricField<NTL::zz_pE>& F =UnparametricField<NTL::zz_pE>(), 
						 const size_t& size = 0,
						 const size_t& seed = 0
						 )
			: _size(size), _seed(seed)
		{
			if(_seed == 0)
				NTL::SetSeed(NTL::to_ZZ(time(0)));
			else
				NTL::SetSeed(NTL::to_ZZ(_seed));
		}

		UnparametricRandIter<NTL::zz_pE>(const UnparametricRandIter<NTL::zz_pE>& R)
			: _size(R._size), _seed(R._seed) 
	
		{
			if(_seed == 0)
				NTL::SetSeed(NTL::to_ZZ(time(0)));
			else
				NTL::SetSeed(NTL::to_ZZ(_seed));
		}
      
		Element& random (Element& x) const
		{
			NTL::random(x);
			return x;
		}

	protected:
		size_t _size;
		size_t _seed;
	};



	/*
	 * Define a parameterized class to easily handle UnparametricField<NTL::zz_pE> field
	 */
	
	/// \brief for large cardinality, small prime.  \ingroup field
 	class NTL_zz_pE : public UnparametricField<NTL::zz_pE>
	{
	public:
		NTL_zz_pE (const integer &p, const integer &k) {
			
			NTL::zz_p::init( (long) p);
			NTL::zz_pX irredPoly = NTL::BuildIrred_zz_pX ((long) k);
			NTL::zz_pE::init(irredPoly);
		}

		Element& random (Element& x) const
		{
			NTL::random(x);
			return x;
		}
		
	}; // end o class NTL_zz_pE
	

	/*
	 * Specialization of UnparametricField<> for NTL::zz_pE type
	 */
	template<>
	NTL::zz_pE& UnparametricField<NTL::zz_pE>::init (NTL::zz_pE &x, const integer &y) const
	{
		x=NTL::to_zz_pE(static_cast<long>(y));
		return x;
	}
	template<>
	NTL::zz_pE& UnparametricField<NTL::zz_pE>::init (NTL::zz_pE &x, const double &y) const
	{
		x=NTL::to_zz_pE(static_cast<long>(y));
		return x;
	}
	
	template<>
	integer& UnparametricField<NTL::zz_pE>::convert (integer& x, const NTL::zz_pE &y) const	{
		NTL::zz_pX poly = rep(y);
		integer base = static_cast<integer>(NTL::zz_p::modulus());		
		long i;		
		x = 0;
		for(i = deg(poly); i >= 0; --i) {
			x *= base;
			x +=  NTL::to_long(rep(coeff(poly, i)));
		}
		return x;
	}
	


	template<>
	bool UnparametricField<NTL::zz_pE>::isZero (const NTL::zz_pE& a) const
	{
		return NTL::IsZero(a);
	}
  
	template<>
	bool UnparametricField<NTL::zz_pE>::isOne (const NTL::zz_pE& a) const
	{
		return NTL::IsOne(a);
	}
  
  
	template<>
	integer& UnparametricField<NTL::zz_pE>::characteristic (integer &c) const
	{
		return c = static_cast<integer>(NTL::zz_p::modulus()); 
	}
  
	template<>
	integer& UnparametricField<NTL::zz_pE>::cardinality(integer& c) const
	{
		NTL::ZZ card = NTL::zz_pE::cardinality();
		long b = NumBytes(card);
		unsigned char* byteArray;
		byteArray = new unsigned char[(size_t)b ];
		BytesFromZZ(byteArray, card, b);
      
		integer base(256);
		c= integer(0);
      
		for(long i = b - 1; i >= 0; --i) {
			c *= base;
			c += integer(byteArray[i]);
		}
		delete [] byteArray;		

		return c;
	}


	template<>
	NTL::zz_pE& UnparametricField<NTL::zz_pE>::inv(NTL::zz_pE& x, const NTL::zz_pE& y) const
	{
		x=NTL::to_zz_pE(1)/y;
		return x;
	}
	template<>
	NTL::zz_pE& UnparametricField<NTL::zz_pE>::invin(NTL::zz_pE& x) const
	{
		x=NTL::to_zz_pE(1)/x;
		return x;
	}


	template<>
	std::istream& UnparametricField<NTL::zz_pE>::read(std::istream& is, NTL::zz_pE& x) const
	{
		long tmp;
		is>>tmp;
		x=NTL::to_zz_pE(tmp);
		return is;
	}
  

}

#endif //__LINBOX_ntl_lzz_pe_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
