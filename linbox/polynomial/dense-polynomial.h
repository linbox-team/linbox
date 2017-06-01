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


namespace LinBox {
    
    template <class BaseRing, class Storage_Tag>
    class PolynomialRing;
    
	/*! Dense Polynomial representation using Givaro.
	 * @ingroup polynomial
	 * A \p dense GivaroPolynomial is an Element of Givaro::Poly1Dom and a reference 
         * to a PolynomialDomain which handles most operations over such polynomials.
	 */
    template<class Field>
    class DensePolynomial : public Givaro::Poly1FactorDom<Field, Givaro::Dense>::Element {

    public:

        typedef DensePolynomial<Field> Self_t;
        typedef Givaro::Poly1FactorDom<Field, Givaro::Dense> Domain_t;

        DensePolynomial (const Field& F) : _field (F) {}
        DensePolynomial (const Field& F, const size_t s) : Domain_t::Element(s), _field (F) {}
        DensePolynomial (const typename Domain_t::Element& P, const Field& F) :
                Domain_t::Element(P),
                _field(F)
            {}

        DensePolynomial& operator=(const DensePolynomial & P)  {
            return *this = P;
        }

        const Domain_t& getField () {return _field;};

        template<typename _Tp1>
        struct rebind {
            typedef DensePolynomial<_Tp1> other;

            void operator() (other & P2,
                             const Self_t& P1,
                             const _Tp1& F)
                {
                    typename Self_t::const_iterator it1 = P1.begin();
                    typename other::iterator it2 = P2.begin();
                    for (; it1 != P1.end(); ++it1, ++it2)
                        F.init (*it2, *it1);
                }
        };

    protected:

        const Field& _field;
    };

} // namespace LinBox

#endif // __LINBOX_polynomial_H
