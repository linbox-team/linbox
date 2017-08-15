
/* tests/__LINBOX_smith_form_kannan_bachem.h
 * Copyright (C) 2014 Gavin Harrison,
 *
 * Written by Gavin Harrison <gmh33@drexel.edu>,
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
 
#include <iostream>
#include "linbox/matrix/dense-matrix.h"

// Givaro Polynomial
#include <givaro/givpoly1.h>
#include "linbox/ring/givaro-poly.h"

// Polynomial Matrix
#include "linbox/matrix/polynomial-matrix.h"
#include "linbox/algorithms/polynomial-matrix/polynomial-matrix-domain.h"

#ifndef __LINBOX_poly_smith_form_iliopoulos_domain_H
#define __LINBOX_poly_smith_form_iliopoulos_domain_H

namespace LinBox
{
	template<class Field>
	class PolynomialSmithFormIliopoulosDomain
	{
	public:
		typedef typename Field::Element Element;
		typedef GivaroPoly<Field> PolyRing;
		typedef typename PolyRing::Element Polynomial;

	private:
		Field _F;
		PolyRing _PD;
		PolynomialMatrixMulDomain<Field> _PMD;
		
		size_t NOT_FOUND = -1;

	public:
		PolynomialSmithFormIliopoulosDomain(const Field &F) : _F(F), _PD(F, "x"), _PMD(F) {
		}

	private:
		
		template<typename PMatrix>
		void writePolynomialToMatrix(PMatrix &M, size_t r, size_t c, const Polynomial &p) const {
			Givaro::Degree d = _PD.degree(p).value();
			
			for (size_t k = 0; k <= M.degree(); k++) {
				if (d >= k) {
					typename PolyRing::Scalar_t e;
					_PD.getEntry(e, Givaro::Degree(k), p);
					_F.assign(M.ref(r, c, k), e);
				}
			}
		}
		
		// Ensures that if a=b then s=u=v=1 and t=0 to avoid an infinite loop
		void dxgcd(Polynomial &s, Polynomial &t, Polynomial &u, Polynomial &v, const Polynomial &a,
			const Polynomial &b) const {
			if (_PD.isDivisor(b, a)) {
				Polynomial quo;
				_PD.quo(quo, b, a);
				
				_PD.assign(s, _PD.one);
				_PD.assign(t, _PD.zero);
				_PD.assign(u, _PD.one);
				_PD.assign(v, quo);
				return;
			}

			Polynomial g;
			_PD.dxgcd(g,s,t,u,v,a,b);
		}
		
		template<typename PMatrix>
		void swapRows(PMatrix& A, size_t a, size_t b) const {
			Element tmp;
			
			for (size_t j = 0; j < A.coldim(); j++) {
				for (size_t k = 0; k < A.degree(); k++) {
					_F.assign(tmp, A.ref(a, j, k));
					_F.assign(A.ref(a, j, k), A.ref(b, j, k));
					_F.assign(A.ref(b, j, k), tmp);
				}
			}
		}
		
		template<typename PMatrix>
		void swapCols(PMatrix& A, size_t a, size_t b) const {
			Element tmp;
			
			for (size_t i = 0; i < A.coldim(); i++) {
				for (size_t k = 0; k < A.degree(); k++) {
					_F.assign(tmp, A.ref(i, a, k));
					_F.assign(A.ref(i, a, k), A.ref(i, b, k));
					_F.assign(A.ref(i, b, k), tmp);
				}
			}
		}
		
		template<typename Vector>
		bool isZero(const Vector &P) const {
			size_t d = P.size();
			
			for (size_t k = 0; k < d; k++) {
				if (!_F.isZero(P[k])) {
					return false;
				}
			}
			
			return true;
		}
		
		void normalize(Polynomial &a, const Polynomial &b, const Polynomial &d) const {
			_PD.gcd(a, b, d);
			
			if (_PD.degree(b).value() == _PD.degree(d).value()) {
				_PD.assign(a, _PD.zero);
			}
		}
		
		// return index of element that divides all other nonzero elements in list
		// or index of zero element if no divisors are found.
		size_t findDivisorOrZero(const std::vector<Polynomial> &lst, const Polynomial &d) const {
			Polynomial divisor;
			size_t divisor_idx = NOT_FOUND;
			size_t zero_idx = NOT_FOUND;
			
			for (size_t i = 0; i < lst.size(); i++) {
				Polynomial tmp;
				normalize(tmp, lst[i], d);
				bool isZero = _PD.isZero(tmp);
				
				if (isZero) {
					zero_idx = zero_idx != NOT_FOUND ? zero_idx : i;
					continue;
				}
				
				if (_PD.isUnit(tmp)) {
					return i;
				}
				
				if (divisor_idx == NOT_FOUND && !isZero) {
					_PD.assign(divisor, tmp);
					divisor_idx = i;
					continue;
				}
				
				if (_PD.isDivisor(tmp, divisor)) {
					continue;
				}
				
				if (_PD.isDivisor(divisor, tmp)) {
					_PD.assign(divisor, tmp);
					divisor_idx = i;
					continue;
				}
				
				for (size_t j = i+1; j < lst.size(); j++) {
					Polynomial tmp2;
					normalize(tmp2, lst[j], d);
					if (_PD.isZero(tmp2)) {
						return j;
					}
				}
				return NOT_FOUND;
			}
			return divisor_idx;
		}
		
		void eliminateCol(PMatrix &M, size_t p, Polynomial &d) {
			std::vector<Polynomial> lst;
			for (size_t i = p; i < M.rowdim(); i++) {
				Polynomial tmp;
				lst.push_back(_PD.init(tmp, M(i, p)));
			}
			size_t p_row = p + findDivisorOrZero(lst, d);
			
			if (p_row == NOT_FOUND) {
				zeroOutPivotRow(M, p, d);
			} else if (p_row != p) {
				swapRows(M, p, p_row);
			}
		}
		
	public:
		template<typename PMatrix>
		void solve(PMatrix &M, Polynomial &d) const {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			std::vector<Polynomial> lst;
			for (size_t i = 0; i < M.rowdim(); i++) {
				Polynomial tmp;
				lst.push_back(_PD.init(tmp, M(i, 0)));
			}
			std::cout << "Idx: " << findDivisorOrZero(lst, d) << std::endl;
		}
	};
}

#endif // __LINBOX_smith_form_iliopoulos_domain_H
