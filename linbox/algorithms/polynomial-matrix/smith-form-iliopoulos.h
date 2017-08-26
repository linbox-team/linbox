
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
				} else {
					_F.assign(M.ref(r, c, k), _F.zero);
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
		
		void xgcd(Polynomial &g, std::vector<Polynomial> &ts, const std::vector<Polynomial> &as, const Polynomial &d) {
			ts.resize(as.size());
			_PD.assign(ts[0], _F.one);
			_PD.assign(g, as[0]);
			
			for (size_t i = 1; i < as.size(); i++) {
				Polynomial r, t1, t2;
				_PD.gcd(r, t1, t2, g, as[i]);
				
				for (size_t j = 0; j < i; j++) {
					_PD.mulin(ts[j], t1);
				}
				
				_PD.assign(ts[i], t2);
				_PD.assign(g, r);
			}
			
			Polynomial r, t1, t2;
			_PD.gcd(r, t1, t2, g, d);
			
			for (size_t i = 0; i < ts.size(); i++) {
				_PD.mulin(ts[i], t1);
			}
			
			_PD.assign(g, r);
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
		
		template<typename PMatrix>
		void reduceMatrix(PMatrix &M, const Polynomial &d) const {
			for (size_t r = 0; r < M.rowdim(); r++) {
				for (size_t c = 0; c < M.coldim(); c++) {
					Polynomial a, b;
					_PD.init(a, M(r, c));
					_PD.rem(b, a, d);
					writePolynomialToMatrix(M, r, c, b);
				}
			}
			
			M.setsize(M.real_degree() + 1);
		}
		
		int64_t degree(const Polynomial &p) const {
			return _PD.degree(p).value();
		}
		
		size_t max(const std::vector<Polynomial> &ps) const {
			int64_t d = 0;
			
			for (size_t i = 0; i < ps.size(); i++) {
				int64_t tmp = degree(ps[i]);
				
				d = d < tmp ? tmp : d;
			}
			
			return d < 0 ? 0 : (size_t) d;
		}
		
		// Column Elimination Methods
		
		template<typename PMatrix>
		size_t findColPivotOrZero(PMatrix &M, size_t p) {
			size_t zeroIdx = NOT_FOUND;
			
			for (size_t i = p; i < M.rowdim(); i++) {
				Polynomial tmp;
				_PD.init(tmp, M(i, p));
				
				if (_PD.isZero(tmp)) {
					zeroIdx = i;
				} else if (degree(tmp) == 0) {
					return i;
				}
			}
			
			return zeroIdx;
		}
		
		template<typename PMatrix>
		void makeLeftElim1(PMatrix &L, size_t pivotRow, size_t otherRow, const Polynomial &s, const Polynomial &t,
			const Polynomial &u, const Polynomial &v) const {
			
			size_t dim = L.rowdim();
			for (size_t i = 0; i < dim; i++) {
				_F.assign(L.ref(i, i, 0), _F.one);
			}
			
			writePolynomialToMatrix(L, pivotRow, pivotRow, s);
			writePolynomialToMatrix(L, pivotRow, otherRow, t);
			
			Polynomial nv;
			_PD.neg(nv, v);
			writePolynomialToMatrix(L, otherRow, pivotRow, nv);
			writePolynomialToMatrix(L, otherRow, otherRow, u);
		}
		
		template<typename PMatrix>
		void zeroPivotCol(PMatrix &M, size_t p) {
			Polynomial pivot;
			Polynomial other;
			
			_PD.init(pivot, M(p, p));
			_PD.init(other, M(p+1, p));
			
			Polynomial s, t, u, v;
			dxgcd(s, t, u, v, pivot, other);
			
			std::vector<Polynomial> tmp;
			size_t size = max(tmp = {s,t,u,v}) + 1;
			PMatrix L(_F, M.rowdim(), M.rowdim(), size);
			makeLeftElim1(L, p, p+1, s, t, u, v);
			
			PMatrix Z(_F, M.coldim(), M.rowdim(), 1);
			_PMD.mul(Z, L, M);
			
			size_t d = Z.real_degree();
			M.setsize(d + 1);
			M.copy(Z, 0, d);
		}
		
		template<typename PMatrix>
		void makeColPivot(PMatrix &M, size_t p, Polynomial &d) {
			std::vector<Polynomial> elms;
			for (size_t i = p+1; i < M.rowdim(); i++) {
				Polynomial tmp;
				_PD.init(tmp, M(i, p));
				elms.push_back(tmp);
			}
			
			Polynomial g;
			std::vector<Polynomial> ts;
			xgcd(g, ts, elms, d);
			
			size_t size = max(ts) + 1;
			PMatrix L(_F, M.rowdim(), M.rowdim(), size);
			
			size_t dim = L.rowdim();
			for (size_t i = 0; i < dim; i++) {
				_F.assign(L.ref(i, i, 0), _F.one);
			}
			
			for (size_t i = 0; i < ts.size(); i++) {
				writePolynomialToMatrix(L, p, p + i + 1, ts[i]);
			}
			
			PMatrix Z(_F, M.rowdim(), M.coldim(), 1);
			_PMD.mul(Z, L, M);
			
			size_t degree = Z.real_degree();
			degree = ((int64_t) degree) < this->degree(g) ? this->degree(g) : degree;
			M.setsize(degree + 1);
			M.copy(Z, 0, degree);
			
			writePolynomialToMatrix(M, p, p, g);
		}
		
		template<typename PMatrix>
		void eliminateCol(PMatrix &M, size_t p, Polynomial &d) {
			size_t zeroOrUnit = findColPivotOrZero(M, p);
			
			if (zeroOrUnit == NOT_FOUND) {
				zeroPivotCol(M, p);
			} else if (zeroOrUnit != p) {
				swapRows(M, p, zeroOrUnit);
			}
			
			if (isZero(M(p, p))) {
				makeColPivot(M, p, d);
			}
			
			// Construct elimination matrix
			Polynomial pivot;
			_PD.init(pivot, M(p, p));
			
			std::vector<Polynomial> qs;
			for (size_t i = p+1; i < M.rowdim(); i++) {
				Polynomial other;
				_PD.init(other, M(i, p));
				
				Polynomial q;
				_PD.quo(q, other, pivot);
				_PD.negin(q);
				
				qs.push_back(q);
			}
			
			size_t size = max(qs) + 1;
			PMatrix L(_F, M.rowdim(), M.rowdim(), size);
			for (size_t i = 0; i < L.rowdim(); i++) {
				_F.assign(L.ref(i, i, 0), _F.one);
			}
			
			for (size_t i = 0; i < qs.size(); i++) {
				writePolynomialToMatrix(L, p + i + 1, p, qs[i]);
			}
			
			PMatrix Z(_F, M.rowdim(), M.coldim(), 1);
			_PMD.mul(Z, L, M);
			
			size_t degree = Z.real_degree();
			M.setsize(degree + 1);
			M.copy(Z, 0, degree);
			
			reduceMatrix(M, d);
		}
		
		// Row Elimination Methods
		
		template<typename PMatrix>
		size_t findRowPivotOrZero(PMatrix &M, size_t p) {
			size_t zeroIdx = NOT_FOUND;
			
			for (size_t i = p; i < M.coldim(); i++) {
				Polynomial tmp;
				_PD.init(tmp, M(p, i));
				
				if (_PD.isZero(tmp)) {
					zeroIdx = i;
				} else if (degree(tmp) == 0) {
					return i;
				}
			}
			
			return zeroIdx;
		}
		
		template<typename PMatrix>
		void makeRightElim1(PMatrix &L, size_t pivotRow, size_t otherRow, const Polynomial &s, const Polynomial &t,
			const Polynomial &u, const Polynomial &v) const {
			
			size_t dim = L.rowdim();
			for (size_t i = 0; i < dim; i++) {
				_F.assign(L.ref(i, i, 0), _F.one);
			}
			
			writePolynomialToMatrix(L, pivotRow, pivotRow, s);
			writePolynomialToMatrix(L, otherRow, pivotRow, t);
			
			Polynomial nv;
			_PD.neg(nv, v);
			writePolynomialToMatrix(L, pivotRow, otherRow, nv);
			writePolynomialToMatrix(L, otherRow, otherRow, u);
		}
		
		template<typename PMatrix>
		void zeroPivotRow(PMatrix &M, size_t p) {
			Polynomial pivot;
			Polynomial other;
			
			_PD.init(pivot, M(p, p));
			_PD.init(other, M(p, p+1));
			
			Polynomial s, t, u, v;
			dxgcd(s, t, u, v, pivot, other);
			
			std::vector<Polynomial> tmp;
			size_t size = max(tmp = {s,t,u,v}) + 1;
			PMatrix R(_F, M.coldim(), M.coldim(), size);
			makeRightElim1(R, p, p+1, s, t, u, v);
			
			PMatrix Z(_F, M.rowdim(), M.coldim(), 1);
			_PMD.mul(Z, M, R);
			
			size_t d = Z.real_degree();
			M.setsize(d + 1);
			M.copy(Z, 0, d);
		}
		
		template<typename PMatrix>
		void makeRowPivot(PMatrix &M, size_t p, Polynomial &d) {
			std::vector<Polynomial> elms;
			for (size_t i = p+1; i < M.coldim(); i++) {
				Polynomial tmp;
				_PD.init(tmp, M(p, i));
				elms.push_back(tmp);
			}
			
			Polynomial g;
			std::vector<Polynomial> ts;
			xgcd(g, ts, elms, d);
			
			size_t size = max(ts) + 1;
			PMatrix R(_F, M.coldim(), M.coldim(), size);
			
			size_t dim = R.rowdim();
			for (size_t i = 0; i < dim; i++) {
				_F.assign(R.ref(i, i, 0), _F.one);
			}
			
			for (size_t i = 0; i < ts.size(); i++) {
				writePolynomialToMatrix(R, p + i + 1, p, ts[i]);
			}
			
			PMatrix Z(_F, M.rowdim(), M.coldim(), 1);
			_PMD.mul(Z, M, R);
			
			size_t degree = Z.real_degree();
			degree = ((int64_t) degree) < this->degree(g) ? this->degree(g) : degree;
			M.setsize(degree + 1);
			M.copy(Z, 0, degree);
			
			writePolynomialToMatrix(M, p, p, g);
		}
		
		template<typename PMatrix>
		void eliminateRow(PMatrix &M, size_t p, Polynomial &d) {
			size_t zeroOrUnit = findRowPivotOrZero(M, p);
			
			if (zeroOrUnit == NOT_FOUND) {
				zeroPivotRow(M, p);
			} else if (zeroOrUnit != p) {
				swapCols(M, p, zeroOrUnit);
			}
			
			if (isZero(M(p, p))) {
				makeRowPivot(M, p, d);
			}
			
			// Construct elimination matrix
			Polynomial pivot;
			_PD.init(pivot, M(p, p));
			
			std::vector<Polynomial> qs;
			for (size_t i = p+1; i < M.coldim(); i++) {
				Polynomial other;
				_PD.init(other, M(p, i));
								
				Polynomial q;
				_PD.quo(q, other, pivot);
				_PD.negin(q);
				
				qs.push_back(q);
			}
			
			size_t size = max(qs) + 1;
			PMatrix R(_F, M.coldim(), M.coldim(), size);
			for (size_t i = 0; i < R.rowdim(); i++) {
				_F.assign(R.ref(i, i, 0), _F.one);
			}
			
			for (size_t i = 0; i < qs.size(); i++) {
				writePolynomialToMatrix(R, p, p + i + 1, qs[i]);
			}
			
			PMatrix Z(_F, M.rowdim(), M.coldim(), 1);
			_PMD.mul(Z, M, R);
			
			size_t degree = Z.real_degree();
			M.setsize(degree + 1);
			M.copy(Z, 0, degree);
			
			reduceMatrix(M, d);
		}
		
		template<typename PMatrix>
		void zeroRow(PMatrix &M, size_t p) {
			for (size_t i = p+1; i < M.coldim(); i++) {
				writePolynomialToMatrix(M, p, i, _PD.zero);
			}
		}
		
		// Move nonzero column to pivot
		
		template<typename PMatrix>
		bool findNonZeroCol(PMatrix &M, size_t p) {
			for (size_t c = p; c < M.coldim(); c++) {
				for (size_t r = p; r < M.rowdim(); r++) {
					if (!isZero(M(r, c))) {
						if (c != p) {
							swapCols(M, p, c);
						}
						return true;
					}
				}
			}
			return false;
		}
		
		template<typename PMatrix>
		bool isDiagonalized(PMatrix &M, size_t p) {
			for (size_t r = p+1; r < M.rowdim(); r++) {
				if (!isZero(M(r, p))) {
					return false;
				}
			}
			
			for (size_t c = p+1; c < M.coldim(); c++) {
				if (!isZero(M(p, c))) {
					return false;
				}
			}
			
			return true;
		}
		
		template<typename PMatrix>
		bool isUnit(PMatrix &M, size_t p, Polynomial &d) {
			Polynomial pivot;
			_PD.init(pivot, M(p, p));
			
			Polynomial g;
			_PD.gcd(g, pivot, d);
			
			return degree(g) == 0;
		}
		
		void fixDiagonal(std::vector<Polynomial> &v) {
			for (size_t i = 0; i < v.size() - 1; i++) {
				Polynomial h;
				_PD.assign(h, v[i]);
				
				_PD.gcd(v[i], v[i+1], h);
				
				_PD.mulin(v[i+1], h);
				_PD.divin(v[i+1], v[i]);
			}
		}
		
	public:
		template<typename PMatrix>
		void solve(std::vector<Polynomial> &v, PMatrix &M, Polynomial &d) {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			for (size_t p = 0; p < dim - 1; p++) {
				if (!findNonZeroCol(M, p)) {
					break;
				}
				
				while (!isDiagonalized(M, p)) {
					eliminateCol(M, p, d);
					
					if (isUnit(M, p, d)) {
						zeroRow(M, p);
					} else {
						eliminateRow(M, p, d);
					}
				}
			}
			
			v.resize(dim);
			for (size_t i = 0; i < dim; i++) {
				Polynomial e;
				_PD.init(e, M(i, i));
				normalize(v[i], e, d);
				
			}
			
			fixDiagonal(v);
		}
	};
}

#endif // __LINBOX_smith_form_iliopoulos_domain_H
