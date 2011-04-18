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


#ifndef __LINBOX_ntl_gf2e_H
#define __LINBOX_ntl_gf2e_H


#include <linbox/util/debug.h>
#include <linbox/field/unparametric.h>
#include <linbox/randiter/unparametric.h>
#include <linbox/util/debug.h>
#include <NTL/GF2XFactoring.h>
#include <NTL/GF2E.h>
#include <time.h>
#include "linbox/linbox-config.h"
#include <linbox/field/field-traits.h>

namespace LinBox 
{

	template <class Ring>
	struct ClassifyRing; 

	template <class Element>
	struct ClassifyRing<UnparametricRandIter<Element> >;
	
	template <>
	struct ClassifyRing<UnparametricRandIter<NTL::GF2E> > {
		typedef RingCategories::ModularTag categoryTag;
	};


	/// \ingroup field
	template<>
	class UnparametricRandIter<NTL::GF2E>
	{
	public:
		typedef NTL::GF2E Element;
		UnparametricRandIter<NTL::GF2E>(const UnparametricField<NTL::GF2E>& F =UnparametricField<NTL::GF2E>(), 
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

		UnparametricRandIter<NTL::GF2E>(const UnparametricRandIter<NTL::GF2E>& R)
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
	 * Define a parameterized class to easily handle UnparametricField<NTL::GF2E> field
	 */
	
 	class NTL_GF2E : public UnparametricField<NTL::GF2E>
	{
	public:
		NTL_GF2E (const integer &p, const integer &k) {	
		  if(p != 2) throw PreconditionFailed(__FUNCTION__,__LINE__,"modulus must be 2");
		  NTL::GF2X irredPoly = NTL::BuildSparseIrred_GF2X((long) k);
		  NTL::GF2E::init(irredPoly);
		}
		
	}; // end o class NTL_GF2E
	

	/*
	 * Specialization of UnparametricField<> for NTL::GF2E type
	 */
	template<>
	NTL::GF2E& UnparametricField<NTL::GF2E>::init (NTL::GF2E &x, const integer &y) const
	{
		x=NTL::to_GF2E(static_cast<long>(y));
		return x;
	}
	template<>
	NTL::GF2E& UnparametricField<NTL::GF2E>::init (NTL::GF2E &x, const double &y) const
	{
		x=NTL::to_GF2E(static_cast<long>(y));
		return x;
	}
	

	template<>
	integer& UnparametricField<NTL::GF2E>::convert (integer& x, const NTL::GF2E &y) const	{
		NTL::GF2X poly = rep(y);
		
		long i;		
		x = 0;
		for(i = deg(poly); i >= 0; --i) {
		  x <<= 1;
		  x += rep(coeff(poly, i));
		}
		return x;
	}
	


	template<>
	bool UnparametricField<NTL::GF2E>::isZero (const NTL::GF2E& a) const
	{
		return NTL::IsZero(a);
	}
  
	template<>
	bool UnparametricField<NTL::GF2E>::isOne (const NTL::GF2E& a) const
	{
		return NTL::IsOne(a);
	}
  
  
	template<>
	integer& UnparametricField<NTL::GF2E>::characteristic (integer &c) const
	{
		return c = 2; 
	}
  
	template<>
	integer& UnparametricField<NTL::GF2E>::cardinality(integer& c) const
	  {
	    c=1;
	    c<<= NTL::GF2E::degree();
	    return c;
	  }


	template<>
	NTL::GF2E& UnparametricField<NTL::GF2E>::inv(NTL::GF2E& x, const NTL::GF2E& y) const
	{
		x=NTL::to_GF2E(1)/y;
		return x;
	}
	template<>
	NTL::GF2E& UnparametricField<NTL::GF2E>::invin(NTL::GF2E& x) const
	{
		x=NTL::to_GF2E(1)/x;
		return x;
	}


	template<>
	std::istream& UnparametricField<NTL::GF2E>::read(std::istream& is, NTL::GF2E& x) const
	{
		long tmp;
		is>>tmp;
		x=NTL::to_GF2E(tmp);
		return is;
	}
  

}

#endif //__LINBOX_ntl_gf2e_H

/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
