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
			_irred(value), 
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
        using Parent_t:: div;
        using Parent_t:: add;
        using Parent_t:: sub;
        using Parent_t:: neg;
        using Parent_t:: inv;

        using Parent_t:: mulin;
        using Parent_t:: divin;
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
	{ std::string s; return Parent_t::read(is >> s)>> s >> _irred >> s >> _exponent; }
        std::ostream& write(std::ostream& os) const
	{ return Parent_t::write(os<<"Local- ") << "irred: " << _irred << ", exponent: " << _exponent; }

	Element& gcdin (Element& a, const Element& b) const
                {
			Givaro::gcd(a, a, b);
                        Givaro::gcd(a, a, Parent_t::residu()); 
			return reduce(a);
                }


	bool isUnit(const Element& a) const
                {
                        Element g;
                        Givaro::gcd(g, a, _irred);
                        return isOne(g);
                }

        
	protected:
		intType _irred;
		uint32_t _exponent;
	};

}


#endif //__LINBOX_local_pir_modular_H


// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
