/* linbox/polynomial/givaro-polynomial.h
 * Copyright(C) The LinBox group
 * Written
 *  by Clement Pernet <clement.pernet@univ-grenoble-alpes.fr>
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

/*! @file polynomial/dense-polynomial.h
 * @ingroup polynomial
 * @brief no doc.
 *
 */

#ifndef __LINBOX_polynomial_H
#define __LINBOX_polynomial_H

#include "givaro/givpoly1.h"
#include "givaro/givpoly1dense.h"
#include "linbox/field/hom.h"
namespace LinBox {
    
    template <class BaseRing, class Storage_Tag>
    class PolynomialRing;
    
	/*! Dense Polynomial representation using Givaro.
     * @ingroup polynomial
     * A \p dense GivaroPolynomial is an Element of Givaro::Poly1Dom and a reference 
     * to a PolynomialDomain which handles most operations over such polynomials.
	 */
    template<class Field>
    class DensePolynomial : public Givaro::Poly1Dom<Field, Givaro::Dense>::Element {

    public:

        typedef DensePolynomial<Field> Self_t;
        typedef Givaro::Poly1Dom<Field, Givaro::Dense> Domain_t;
        typedef typename Domain_t::Element Storage_t;

        DensePolynomial () : Domain_t::Element(), _field(NULL) {}

        DensePolynomial (const Field& F) : _field (&F) {}

        template<typename... Args>
        DensePolynomial(const Field& F, Args... args) :
                Domain_t::Element(args...),
                _field(&F) {}

        template <class _OtherPoly >
        DensePolynomial (const _OtherPoly& P, const Field& F) :
                Domain_t::Element(P.size()),
                _field(&F)
            {
                typename _OtherPoly::template rebind<Field>()(*this,P);
            }

        const Field& field () const {return *_field;};

        template<typename _Tp1>
        struct rebind {
            typedef DensePolynomial<_Tp1> other;

            void operator() (other & P2, const Self_t& P1) {
                typedef typename Self_t::const_iterator ConstSelfIterator ;
                typedef typename other::iterator OtherIterator ;
                OtherIterator    P2_i = P2.begin();
                ConstSelfIterator P1_i = P1.begin();
                Hom<Field, _Tp1> hom (P1.field(), P2.field()) ;
                for ( ; P1_i != P1.end(); ++P1_i, ++P2_i)
                    hom.image (*P2_i, *P1_i);
			}
        };

    protected:
	// DensePolynomials can be constructed with no Field
	// Therefore PolynomialRing will set the Field pointer from the "init" member
        template <class BaseRing, class Storage_Tag> friend class PolynomialRing;
        const Field* _field;
    };

} // namespace LinBox

#endif // __LINBOX_polynomial_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
