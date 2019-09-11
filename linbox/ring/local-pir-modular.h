/* Copyright (C) 2010 LinBox
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


#ifndef __LINBOX_local_pir_modular_H
#define __LINBOX_local_pir_modular_H

#include <string>
#include <givaro/modular.h>
#include <givaro/givpower.h>

#include "linbox/field/field-traits.h"

// Namespace in which all LinBox code resides
namespace LinBox
{


	template <typename intType>
	class LocalPIRModular;

	template <typename intType>
	struct ClassifyRing<LocalPIRModular<intType> >  {
		typedef RingCategories::ModularTag categoryTag;
	};

        /// \ingroup ring
	template <typename intType>
	class LocalPIRModular : public Givaro::Modular<intType>  {

	public:

		using Parent_t = Givaro::Modular<intType>;

		using typename Givaro::Modular<intType>::Element;
            //typedef typename Givaro::Modular<intType>::Element Element;

		using typename Givaro::Modular<intType>::RandIter;

            //No default cstor

		LocalPIRModular (intType value, uint32_t exp = 1) :
                Givaro::Modular<intType>(Givaro::power(value, exp)),
			_exponent(exp)
            {}

        ~LocalPIRModular() noexcept {};

        using Parent_t:: zero;
        using Parent_t:: one;
        using Parent_t:: mOne;

        using Parent_t:: minElement;
        using Parent_t:: maxElement;

        using Parent_t:: characteristic;
        using Parent_t:: cardinality;

        using Parent_t:: maxCardinality;
        using Parent_t:: minCardinality;

            // ----- Checkers
        using Parent_t:: isZero;
        using Parent_t:: isOne ;
        using Parent_t:: isMOne;
        using Parent_t:: areEqual;
        using Parent_t:: length;

            // ----- Ring-wise operators
        using Parent_t:: operator==;
        using Parent_t:: operator!=;
        using Parent_t:: operator=;

            // ----- Initialisation
        using Parent_t:: init ;

        using Parent_t:: assign ;

        using Parent_t:: convert;

        using Parent_t:: reduce ;

        using Parent_t:: mul;
        using Parent_t:: add;
        using Parent_t:: sub;
        using Parent_t:: neg;
        using Parent_t:: inv;

        using Parent_t:: mulin;
        using Parent_t:: addin;
        using Parent_t:: subin;
        using Parent_t:: negin;
        using Parent_t:: invin;

        using Parent_t:: axpy  ;
        using Parent_t:: axpyin;

        using Parent_t:: axmy;
        using Parent_t:: axmyin;

        using Parent_t:: maxpy;
        using Parent_t:: maxpyin;

            // ----- Random generators
        using typename Parent_t:: NonZeroRandIter;
        using Parent_t:: random;
        using Parent_t:: nonzerorandom;

            // --- IO methods
        using Parent_t:: read ;
        using Parent_t:: write;

        std::istream& read (std::istream& is)
            { std::string s; return Parent_t::read(is >> s)>> s >> Parent_t::residu() >> s >> _exponent; }
        std::ostream& write(std::ostream& os) const
            { return Parent_t::write(os<<"Local- ") << "irred: " << Parent_t::residu() << ", exponent: " << _exponent; }

        Element& gcdin (Element& a, const Element& b) const {
            Givaro::gcd(a, a, b);
            Givaro::gcd(a, a, Parent_t::residu());
            return reduce(a);
        }

        Element& gcd(Element& g, const Element& a, const Element& b) const {
            return Givaro::gcd(g,a,b);
        }

        Element& xgcd(Element& g, Element& s, Element& t, const Element& a, const Element& b) const {
            return Givaro::gcd(g,s,t,a,b);
        }


        bool isUnit(const Element& a) const {
            Element g;
            Givaro::gcd(g, a, Parent_t::residu());
            return isOne(g);
        }

        bool isDivisor(const Element& a, const Element& b) const {
            Integer g;
            if (this->isZero(a)) return false;
            else if (this->isZero(b)) return true;
            else {
                Givaro::gcd(g, Integer(a), Integer(Parent_t::residu()));
                return this->isZero(Integer(b) % g);
            }
        }
        Element& div(Element& r, const Element& a, const Element& b) const {
            Integer g, ia(a), ib(b);
            Givaro::gcd(g, ia, ib);
            ia /= g;
            ib /= g;
            Element iv;
            this->inv(iv, ib);
            return this->mul(r, ia, iv);
        }

        Element& divin(Element& r, const Element& b) const {
            Element a(r); // @enhancement: could use inplace variants
            return div(r,a,b);
        }

		Element& normal (Element& a, const Element& b) const
		{

			if (this->isZero(b))
                return a = this->zero;
			else
				return this->gcd(a, b, Parent_t::residu() );
		}

		Element& normalIn (Element& a) const
		{
			if (this->isZero(a))
                return a;
            else
				return this->gcdin(a, Parent_t::residu() );
		}


	protected:
		uint32_t _exponent;
	};

}


#endif //__LINBOX_local_pir_modular_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
