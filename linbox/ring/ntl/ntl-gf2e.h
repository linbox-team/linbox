/* linbox/field/ntl-z_pE.h
 * Copyright (C) 2004  Pascal Giorgi
 * Copyright (C) 2011 LinBox
 *
 * Written by  Pascal Giorgi <pascal.giorgi@ens-lyon.fr>
 *
 * Modified by W. J. Turner <wjturner@acm.org>
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

/*! @file field/ntl/ntl-GF2E.h
 * @ingroup field
 * @ingroup NTL
 * @brief NO DOC
 */

#ifndef __LINBOX_field_ntl_gf2e_H
#define __LINBOX_field_ntl_gf2e_H

#ifndef __LINBOX_HAVE_NTL
#error "you need NTL here"
#endif

#include <NTL/GF2XFactoring.h>
#include <NTL/GF2E.h>
#include <time.h>

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include <givaro/zring.h>
#include "linbox/field/field-traits.h"


#include "linbox/integer.h"
#include "ntl-zz.h"

namespace Givaro
{
	template<>
        //	NTL::GF2E& Givaro::ZRing<NTL::GF2E>::init (NTL::GF2E &x, const Integer &y) const
	NTL::GF2E& Caster(NTL::GF2E &x, const Integer &y)
	{
		x=NTL::to_GF2E(static_cast<int64_t>(y));
		return x;
	}
	template<>
        //	NTL::GF2E& Givaro::ZRing<NTL::GF2E>::init (NTL::GF2E &x, const double &y) const
	NTL::GF2E& Caster(NTL::GF2E &x, const double &y)
	{
		x=NTL::to_GF2E(static_cast<long>(y));
		return x;
	}


	template<>
        //	Integer& Givaro::ZRing<NTL::GF2E>::convert (Integer& x, const NTL::GF2E &y) const	{
	Integer& Caster(Integer& x, const NTL::GF2E &y)
	{
		NTL::GF2X poly = rep(y);

		long i;
		x = 0;
		for(i = deg(poly); i >= 0; --i) {
			x <<= 1;
			x += static_cast<int64_t>(rep(coeff(poly, i)));
		}
		return x;
	}

} // namespace Givaro




// Namespace in which all LinBox library code resides
namespace LinBox
{
        /// \ingroup field

class NTL_GF2E_Initialiser {
public :
    NTL_GF2E_Initialiser( const Integer & p, size_t k = 1) {
        linbox_check(p == 2);
			// if(p != 2) throw PreconditionFailed(LB_FILE_LOC,"modulus must be 2");
        NTL::GF2X irredPoly = NTL::BuildSparseIrred_GF2X((long) k);
        NTL::GF2E::init(irredPoly);

    }

		// template <class ElementInt>
		// NTL_GF2E_Initialiser(const ElementInt& d) {
        // NTL::ZZ_p::init (NTL::to_ZZ(d));
		// }

		// NTL_GF2E_Initialiser (const NTL::ZZ& d) {
        // NTL::ZZ_p::init(d);
		// }

};


        /*
         * Define a parameterized class to easily handle Givaro::ZRing<NTL::GF2E> field
         */

	class NTL_GF2E :  public NTL_GF2E_Initialiser, public Givaro::UnparametricOperations<NTL::GF2E> {
	public:
		typedef NTL::GF2E Element ;
		typedef Givaro::UnparametricOperations<Element> Father_t ;
		typedef UnparametricRandIter<Element> RandIter;

		const Element zero,one,mOne ;

		NTL_GF2E (const integer &p, const int32_t &k) :
                NTL_GF2E_Initialiser(p,k),Father_t ()
			,zero( NTL::to_GF2E(0)),one( NTL::to_GF2E(1)),mOne(-one)
            { }

		bool isZero (const Element& a) const
            {
                return NTL::IsZero(a);
            }


		bool isOne (const Element& a) const
            {
                return NTL::IsOne(a);
            }

		bool isUnit (const Element& a) const
            {
                return (deg(rep(a))==0 &&
                        NTL::IsOne(ConstTerm(rep(a))));
            }

		bool isMOne (const Element& x) const
            {
                Element y ; neg(y,x);
                return isOne(y);
            }

		integer& characteristic (integer &c) const
            {
                return c = 2;
            }

		integer& cardinality(integer& c) const
            {
                c=1;
                c<<= static_cast<int64_t>(Element::degree());
                return c;
            }

		Element& inv(Element& x, const Element& y) const
            {
                x=NTL::to_GF2E(1)/y;
                return x;
            }

		Element& invin(Element& x) const
            {
                x=NTL::to_GF2E(1)/x;
                return x;
            }

		std::istream& read(std::istream& is, Element& x) const
            {
                long tmp= 0;
                is>>tmp;
                x=NTL::to_GF2E(tmp);
                return is;
            }
	}; // end o class NTL_GF2E

	template <class Ring>
	struct ClassifyRing;

	template <class Element>
	struct ClassifyRing<UnparametricRandIter<Element> >;

	template <>
	struct ClassifyRing<UnparametricRandIter<NTL::GF2E> > {
		typedef RingCategories::ModularTag categoryTag;
	};

	template<>
	class UnparametricRandIter<NTL::GF2E> {
	public:
		typedef NTL::GF2E Element;
        typedef size_t Residu_t;

		UnparametricRandIter(const NTL_GF2E & F,
                                     const uint64_t seed = 0,
                                     const Residu_t& size = 0) :
                _size(size), _seed(seed)
            {
                if(_seed == 0)
                    NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(std::time(nullptr))));
                else
                    NTL::SetSeed(Caster<NTL::ZZ,uint64_t>(_seed));
            }

		UnparametricRandIter(const UnparametricRandIter<NTL::GF2E>& R) :
                _size(R._size), _seed(R._seed)

            {
                if(_seed == 0)
                    NTL::SetSeed(NTL::to_ZZ(static_cast<long unsigned int>(std::time(nullptr))));
                else
                    NTL::SetSeed(Caster<NTL::ZZ,uint64_t>(_seed));
            }

		Element& random (Element& x) const
            {
                NTL::random(x);
                return x;
            }

	protected:
		size_t _size;
		uint64_t _seed;
	};


}

#endif //__LINBOX_ntl_gf2e_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
