/* Copyright (C) 2007  LinBox
 * Written by JG Dumas
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

#ifndef __LINBOX_rational_full_multip_cra_H
#define __LINBOX_rational_full_multip_cra_H

#include "givaro/zring.h"
#include "linbox/algorithms/cra-builder-full-multip.h"

namespace LinBox
{

	template<class Domain_Type>
	struct RationalCRABuilderFullMultip : public virtual CRABuilderFullMultip<Domain_Type> {
		typedef Domain_Type				Domain;
		typedef CRABuilderFullMultip<Domain> 			Father_t;
		typedef typename Father_t::DomainElement 	DomainElement;
		typedef RationalCRABuilderFullMultip<Domain>		Self_t;
		Givaro::ZRing<Integer> _ZZ;
	public:

		RationalCRABuilderFullMultip(const double log2Bound = 0.0) :
			Father_t(log2Bound)
		{}


        template <class Vect>
		Vect& result (Vect &num, Integer& den)
		{
            Father_t::result(num, false);
            den = 1;
            const auto& mod = Father_t::getModulus();
            Integer s, nd;
            _ZZ.sqrt(s, mod);
            for (auto num_it = num.begin(); num_it != num.end(); ++num_it) {
                iterativeratrecon(*num_it, nd, den, mod, s);

                if (nd > 1) {
                    for (auto t02 = num.begin(); t02 != num_it; ++t02)
                        *t02 *= nd;
                    den *= nd;
                }
            }
            return num;
        }

	protected:
		Integer& iterativeratrecon(Integer& u1, Integer& new_den, const Integer& old_den, const Integer& m1, const Integer& s)
		{
			Integer a;
			_ZZ.RationalReconstruction(a, new_den, u1*=old_den, m1, s);
			return u1=a;
		}
	};
}

#endif //__LINBOX_rational_full_multip_cra_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
