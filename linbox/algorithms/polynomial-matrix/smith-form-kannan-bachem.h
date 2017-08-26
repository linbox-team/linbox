
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

#ifndef __LINBOX_poly_smith_form_kannan_bachem_domain_H
#define __LINBOX_poly_smith_form_kannan_bachem_domain_H

namespace LinBox
{
	template<class Field>
	class PolynomialSmithFormKannanBachemDomain
	{
	public:
		typedef typename Field::Element Element;
		typedef GivaroPoly<Field> PolyRing;
		typedef typename PolyRing::Element Polynomial;

	private:
		Field _F;
		PolyRing _PD;
		PolynomialMatrixMulDomain<Field> _PMD;

	public:
		PolynomialSmithFormKannanBachemDomain(const Field &F) : _F(F), _PD(F, "x"), _PMD(F) {
		}

	private:
		
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
		
		template<typename PMatrix>
		bool findPivot(PMatrix &M, size_t p) const {
			size_t r = M.rowdim();
			size_t c = M.coldim();
			
			for (size_t i = p; i < r; i++) {
				for (size_t j = p; j < c; j++) {
					if (!isZero(M(i, j))) {
						if (i != p) {
							swapRows(M, i, p);
						}
						
						if (j != p) {
							swapCols(M, j, p);
						}
						
						return true;
					}
				}
			}
			
			return false;
		}
		
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
		
		template<typename PMatrix>
		void makeLeftElim(PMatrix &L, size_t pivotRow, size_t otherRow, const Polynomial &s, const Polynomial &t,
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
		
		int64_t degree(const Polynomial &p) const {
			return _PD.degree(p).value();
		}
		
		size_t max(const Polynomial &a, const Polynomial &b, const Polynomial &c, const Polynomial &d) const {
			int64_t ad, bd, cd, dd;
			ad = degree(a);
			bd = degree(b);
			cd = degree(c);
			dd = degree(d);
			
			int64_t r = ad < bd ? bd : ad;
			r = r < cd ? cd : r;
			r = r < dd ? dd : r;
			return r < 0 ? 0 : (size_t) r;
		}
		
		// use row pivotRow to eliminate row otherRow in the pivotRow-th column
		template<typename PMatrix>
		void eliminateCol(PMatrix &M, size_t pivotRow, size_t otherRow) const {
			Polynomial other;
			_PD.init(other, M(otherRow, pivotRow));
			
			if (_PD.isZero(other)) {
				return;
			}
			
			Polynomial pivot;
			_PD.init(pivot, M(pivotRow, pivotRow));
			
			Polynomial s, t, u, v;
			dxgcd(s, t, u, v, pivot, other);
			
			size_t size = max(s,t,u,v) + 1;
			PMatrix L(_F, M.rowdim(), M.rowdim(), size);
			makeLeftElim(L, pivotRow, otherRow, s, t, u, v);
			
			PMatrix Z(_F, M.coldim(), M.rowdim(), 1);
			_PMD.mul(Z, L, M);
			
			size_t d = Z.real_degree();
			M.setsize(d + 1);
			M.copy(Z, 0, d);
		}
		
		template<typename PMatrix>
		void eliminateCol(PMatrix &M, size_t pivotRow) const {
			for (size_t otherRow = pivotRow + 1; otherRow < M.rowdim(); otherRow++) {
				eliminateCol(M, pivotRow, otherRow);
			}
		}
		
		template<typename PMatrix>
		void makeRightElim(PMatrix &R, size_t pivotCol, size_t otherCol, const Polynomial &s, const Polynomial &t,
			const Polynomial &u, const Polynomial &v) const {
			
			size_t dim = R.rowdim();
			for (size_t i = 0; i < dim; i++) {
				_F.assign(R.ref(i, i, 0), _F.one);
			}
			
			writePolynomialToMatrix(R, pivotCol, pivotCol, s);
			writePolynomialToMatrix(R, otherCol, pivotCol, t);
			
			Polynomial nv;
			_PD.neg(nv, v);
			writePolynomialToMatrix(R, pivotCol, otherCol, nv);
			writePolynomialToMatrix(R, otherCol, otherCol, u);
		}
		
		// use col pivotCol to eliminate col otherCol in the pivotCol-th row
		template<typename PMatrix>
		void eliminateRow(PMatrix &M, size_t pivotCol, size_t otherCol) const {
			Polynomial other;
			_PD.init(other, M(pivotCol, otherCol));
			
			if (_PD.isZero(other)) {
				return;
			}
			
			Polynomial pivot;
			_PD.init(pivot, M(pivotCol, pivotCol));
			
			Polynomial s, t, u, v;
			dxgcd(s, t, u, v, pivot, other);
			
			size_t size = max(s,t,u,v) + 1;
			PMatrix R(_F, M.coldim(), M.coldim(), size);
			makeRightElim(R, pivotCol, otherCol, s, t, u, v);
			
			PMatrix Z(_F, M.coldim(), M.rowdim(), 1);
			_PMD.mul(Z, M, R);
			
			size_t d = Z.real_degree();
			M.setsize(d + 1);
			M.copy(Z, 0, d);
		}
		
		template<typename PMatrix>
		void eliminateRow(PMatrix &M, size_t pivotCol) const {
			for (size_t otherCol = pivotCol + 1; otherCol < M.coldim(); otherCol++) {
				eliminateRow(M, pivotCol, otherCol);
			}
		}

		// True if pivot row/col is zero for all indexes greater other than pivot
		template<typename PMatrix>
		bool isDiagonalized(const PMatrix &M, size_t pivot) const {
			for (size_t i = pivot + 1; i < M.coldim(); i++) {
				Polynomial elm;
				_PD.init(elm, M(pivot, i));
				
				if (!_PD.isZero(elm)) {
					return false;
				}
			}
			
			for (size_t i = pivot + 1; i < M.rowdim(); i++) {
				Polynomial elm;
				_PD.init(elm, M(i, pivot));
				
				if (!_PD.isZero(elm)) {
					return false;
				}
			}
			
			return true;
		}
		
		template<typename PMatrix>
		void fixDiagonal(PMatrix &M) const {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			bool done;
			do {
				done = true;
				for (size_t i = 0; i < dim - 1; i++) {
					Polynomial other;
					_PD.init(other, M(i + 1, i + 1));
					
					if (_PD.isZero(other)) {
						continue;
					}
					
					Polynomial pivot;
					_PD.init(pivot, M(i, i));
					
					Polynomial g;
					_PD.gcd(g, pivot, other);
					
					if (_PD.areAssociates(g, pivot)) {
						continue;
					}
					
					writePolynomialToMatrix(M, i, i, g);
					done = false;
				}
			} while (!done);
		}
		
		template<typename PMatrix>
		void makeReduce(PMatrix &L, size_t row, size_t col, const Polynomial &quo) const {
			size_t dim = L.rowdim() < L.coldim() ? L.rowdim() : L.coldim();
			
			for (size_t i = 0; i < dim; i++) {
				_F.assign(L.ref(i, i, 0), _F.one);
			}
			
			Polynomial nquo;
			_PD.neg(nquo, quo);
			writePolynomialToMatrix(L, row, col, nquo);
		}
		
		/* 
		 * Given a upper triangular matrix, reduces the off diagonal elements to be of degree less than
		 * the element on the diagonal below it, i.e., M[1,3] = M[1,3] - (M[1,3] / M[3,3]) * M[3,3]
		 * 
		 * Uses Chou/Collins ordering:
		 * * 4 5 6
		 * 0 * 2 3
		 * 0 0 * 1
		 * 0 0 0 *    
		 */
		template<typename PMatrix>
		void reduceOffDiagonal(PMatrix &M) const {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			for (size_t row = dim-1; 0 < row; row--) {
				for (size_t col = row+1; col < M.coldim(); col++) {
					Polynomial pivot, other;
					_PD.init(pivot, M(col, col));
					_PD.init(other, M(row, col));
					
					if (degree(other) < degree(pivot)) {
						continue;
					}
					
					Polynomial q;
					_PD.quo(q, other, pivot);
					
					PMatrix L(_F, M.rowdim(), M.rowdim(), degree(q) + 1);
					makeReduce(L, row, col, q);
					
					PMatrix Z(_F, M.rowdim(), M.coldim(), 1);
					_PMD.mul(Z, L, M);
					
					size_t d = Z.real_degree();
					M.setsize(d + 1);
					M.copy(Z, 0, d);
				}
			}
		}
		
		// Makes the bottom right n-by-n into a hermite matrix.
		template<typename PMatrix>
		void hermite(PMatrix &M, size_t n) const {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			for (size_t pivot = n; pivot < dim; pivot++) {
				if (!findPivot(M, pivot)) {
					break;
				}
				
				eliminateCol(M, pivot);
			}
			
			reduceOffDiagonal(M);
		}
		
	public:
		
		template<typename PMatrix>
		void solve(PMatrix &M) const {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			for (size_t pivot = 0; pivot < dim; pivot++) {
				if (!findPivot(M, pivot)) {
					break;
				}
				
				while (!isDiagonalized(M, pivot)) {
					eliminateRow(M, pivot);
					hermite(M, pivot);
				}
			}
			
			fixDiagonal(M);
		}
		
		template<typename PMatrix>
		void solveTextbook(PMatrix &M) const {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			for (size_t pivot = 0; pivot < dim; pivot++) {
				if(!findPivot(M, pivot)) {
					break;
				}
				
				while (!isDiagonalized(M, pivot)) {
					eliminateCol(M, pivot);
					eliminateRow(M, pivot);
				}
			}
			
			fixDiagonal(M);
		}
	};
}

#endif // __LINBOX_smith_form_kannan_bachem_domain_H
