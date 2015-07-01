/* linbox/field/givaro-field.h
 * Copyright (C) 2009 JGD
 *
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

#ifndef __LINBOX_givaro_field_gf2_H
#define __LINBOX_givaro_field_gf2_H

#include "linbox/field/gf2.h"
#include "linbox/field/Givaro/givaro-field.h"

// Specialization of GivaroField for GF2
namespace LinBox
{

	/**
	  \brief give LinBox fields an allure of Givaro Fields
	  \ingroup field

	 *  This class adds the necessary requirements allowing
	 *  the construction of an extension of a LinBox field.
	 */
	template<>
	struct GivaroField<LinBox::GF2> : public LinBox::GF2 {
		typedef LinBox::GF2 BaseField;
		typedef BaseField::Element TT;
		typedef Signed_Trait<TT>::unsigned_type UTT;
		typedef TT Rep;
		typedef GivaroField<BaseField> Self_t;
		typedef Rep Element;
		typedef UTT Residu_t;

		Element zero, one;
		GivaroField(const BaseField& bf) :
			BaseField(bf)
		{
			this->init(zero,0UL);
			this->init(one, 1UL);
		}


		// -- amxy: r <- c - a * b mod p
		Rep& amxy (Rep& r, const Rep a, const Rep b, const Rep c) const
		{
			Rep tmp;
			this->mul(tmp, a, b);
			this->assign(r,c);
			return this->subin(r,tmp);
		}
		stdBitReference amxy (stdBitReference r, const Rep a, const Rep b, const Rep c) const
		{
			Rep tmp;
			this->mul(tmp, a, b);
			this->assign(r,c);
			return this->subin(r,tmp);
		}


		// -- maxpy: r <- y - a * x
		Rep& maxpy (Rep& r, const Rep a, const Rep x, const Rep y) const
		{
			Rep tmp; this->mul(tmp, a, x);
			return this->sub(r,y,tmp);
		}
		stdBitReference maxpy (stdBitReference r, const Rep a, const Rep x, const Rep y) const
		{
			Rep tmp; this->mul(tmp, a, x);
			return this->sub(r,y,tmp);
		}
		// -- axmyin: r <-  a * x - r
		Rep& axmyin (Rep& r, const Rep a, const Rep x) const
		{
			maxpyin(r,a,x);
			return negin(r);
		}
		stdBitReference axmyin (stdBitReference r, const Rep a, const Rep x) const
		{
			maxpyin(r,a,x);
			return negin(r);
		}
		// -- maxpyin: r <- r - a * x
		Rep& maxpyin (Rep& r, const Rep a, const Rep x) const
		{
			Rep tmp; this->mul(tmp, a, x);
			return this->subin(r,tmp);
		}
		stdBitReference maxpyin (stdBitReference r, const Rep a, const Rep x) const
		{
			Rep tmp; this->mul(tmp, a, x);
			return this->subin(r,tmp);
		}



		bool areNEqual ( const Rep a, const Rep b) const
		{
			return ! this->areEqual(a,b);
		}

		// Access to the modulus, characteristic, size, exponent
		UTT residu() const
		{
			integer c;
			BaseField::characteristic(c);
			return UTT(c);
		}

		UTT characteristic() const
		{
			integer c; BaseField::characteristic(c); return UTT(c);
		}
		UTT cardinality() const
		{
			integer c; BaseField::cardinality(c); return UTT(c);
		}
		UTT exponent() const
		{
			return 1;
		}
		UTT size() const
		{
			integer c;
			BaseField::cardinality(c);
			return UTT(c);
		}


		// ----- random generators
		template<class RandIter> Rep& random(RandIter& g, Rep& r) const
		{
			return r = g() ;
		}
		template<class RandIter> Rep& random(RandIter& g, Rep& r, long s) const
		{
			return r = g() ;
		}
		template<class RandIter> Rep& random(RandIter& g, Rep& r, const Rep& b) const
		{
			return r = g() ;
		}
		template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const
		{
			return r = g() ;
		}
		template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, long s) const
		{
			return r = g() ;
		}
		template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, const Rep& b) const
		{
			return r = g() ;
		}

		template<class RandIter> stdBitReference random(RandIter& g, stdBitReference r) const
		{
			return r = g() ;
		}
		template<class RandIter> stdBitReference random(RandIter& g, stdBitReference r, long s) const
		{
			return r = g() ;
		}
		template<class RandIter> stdBitReference random(RandIter& g, stdBitReference r, const stdBitReference b) const
		{
			return r = g() ;
		}
		template<class RandIter> stdBitReference nonzerorandom(RandIter& g, stdBitReference r) const
		{
			return r = g() ;
		}
		template<class RandIter> stdBitReference nonzerorandom(RandIter& g, stdBitReference r, long s) const
		{
			return r = g() ;
		}
		template<class RandIter> stdBitReference nonzerorandom(RandIter& g, stdBitReference r, const Rep& b) const
		{
			return r = g() ;
		}
		template<class RandIter> stdBitReference nonzerorandom(RandIter& g, stdBitReference r, const stdBitReference b) const
		{
			return r = g() ;
		}

	};

}

#endif // __LINBOX_givaro_field_gf2_H


// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

