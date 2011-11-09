/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
/* linbox/field/givaro-field.h
 * Copyright (C) 2009 JGD
 *
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#ifndef __LINBOX_givaro_field_H
#define __LINBOX_givaro_field_H

#include <givaro/givconfig.h>
#include "linbox/integer.h"
#include "linbox/field/field-traits.h"
#include "linbox/field/field-interface.h"

namespace LinBox
{

	/**
	  \brief give LinBox fields an allure of Givaro Fields
	  \ingroup field

	 *  This class adds the necessary requirements allowing
	 *  the construction of an extension of a LinBox field
	 *  or a givaro polynomial of a LinBox field ...
	 */
	template< class BaseField >
	struct GivaroField : public BaseField {
		typedef typename BaseField::Element TT;
		typedef typename Signed_Trait<TT>::unsigned_type UTT;
		typedef TT Rep;
		typedef GivaroField<BaseField> Self_t;
		typedef Rep Element;
		typedef UTT Residu_t;

		Element zero, one, mOne;
		GivaroField() :
			BaseField()
		{
			this->init(zero,0UL);
			this->init(one, 1UL);
			init(mOne);
			neg(mOne,one);
		}

		GivaroField(const BaseField& bf) :
			BaseField(bf)
		{
			this->init(zero,0UL);
			this->init(one, 1UL);
			init(mOne);
			neg(mOne,one);
			// neg(init(mOne),one);
		}


		// -- axmy: r <- a * x - y
		Rep& axmy (Rep& r, const Rep a, const Rep x, const Rep y) const
		{
			Rep tmp; this->mul(tmp, a, x);
			return this->sub(r,tmp,y);
		}

		// -- maxpy: r <- y - a * x
		Rep& maxpy (Rep& r, const Rep a, const Rep x, const Rep y) const
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

		// -- maxpyin: r <- r - a * x
		Rep& maxpyin (Rep& r, const Rep a, const Rep x) const
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
			integer c;
			BaseField::characteristic(c);
			return UTT(c);
		}

		integer& characteristic(integer& i) const
		{
			return BaseField::characteristic(i);
		}

		UTT cardinality() const
		{
			integer c;
			BaseField::cardinality(c);
			return UTT(c);
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
			return this->init(r,g()) ;
		}

		template<class RandIter> Rep& random(RandIter& g, Rep& r, long s) const
		{
			return this->init(r,g()) ;
		}

		template<class RandIter> Rep& random(RandIter& g, Rep& r, const Rep& b) const
		{
			return this->init(r,g()) ;
		}

		template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r) const
		{
			do { this->init(r,g()); } while( this->isZero(r) );
			return r;
		}

		template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, long s) const
		{
			do { this->init(r,g()); } while( this->isZero(r) );
			return r;
		}

		template<class RandIter> Rep& nonzerorandom(RandIter& g, Rep& r, const Rep& b) const
		{
			do { this->init(r,g()); } while( this->isZero(r) ); return r;
		}

	};



} // end namespace LinBox

// Specialization of GivaroField for GF2
#ifdef __LINBOX_field_gf2_H
namespace LinBox
{

	/**
	  \brief give LinBox fields an allure of Givaro Fields
	  \ingroup field

	 *  This class adds the necessary requirements allowing
	 *  the construction of an extension of a LinBox field.
	 */
	template<>
	struct GivaroField<LinBox::GF2> : public LinBox::GF2
	{
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
		std::_Bit_reference amxy (std::_Bit_reference r, const Rep a, const Rep b, const Rep c) const
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
		std::_Bit_reference maxpy (std::_Bit_reference r, const Rep a, const Rep x, const Rep y) const
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
		std::_Bit_reference axmyin (std::_Bit_reference r, const Rep a, const Rep x) const
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
		std::_Bit_reference maxpyin (std::_Bit_reference r, const Rep a, const Rep x) const
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
		{ integer c; BaseField::cardinality(c); return UTT(c);
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

		template<class RandIter> std::_Bit_reference random(RandIter& g, std::_Bit_reference r) const
		{
			return r = g() ;
		}
		template<class RandIter> std::_Bit_reference random(RandIter& g, std::_Bit_reference r, long s) const
		{
			return r = g() ;
		}
		template<class RandIter> std::_Bit_reference random(RandIter& g, std::_Bit_reference r, const std::_Bit_reference b) const
		{
			return r = g() ;
		}
		template<class RandIter> std::_Bit_reference nonzerorandom(RandIter& g, std::_Bit_reference r) const
		{
			return r = g() ;
		}
		template<class RandIter> std::_Bit_reference nonzerorandom(RandIter& g, std::_Bit_reference r, long s) const
		{
			return r = g() ;
		}
		template<class RandIter> std::_Bit_reference nonzerorandom(RandIter& g, std::_Bit_reference r, const Rep& b) const
		{
			return r = g() ;
		}
		template<class RandIter> std::_Bit_reference nonzerorandom(RandIter& g, std::_Bit_reference r, const std::_Bit_reference b) const
		{
			return r = g() ;
		}

	};

}

#endif

#endif // __LINBOX_givaro_field_H

