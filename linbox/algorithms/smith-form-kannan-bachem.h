
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
#include "linbox/ring/givaro-poly.h"
#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"

#ifndef __LINBOX_smith_form_kannan_bachem_domain_H
#define __LINBOX_smith_form_kannan_bachem_domain_H

namespace LinBox
{
	template<class Field>
	class SmithFormKannanBachemDomain
	{
	public:
		typedef typename Field::Element Element;
		
		typedef MatrixDomain<Field> MatrixDom;
		
		typedef typename MatrixDom::OwnMatrix OwnMatrix;
		typedef typename MatrixDom::Matrix SubMatrix;

	private:
		Field _F;
		MatrixDom _MD;

	public:
		SmithFormKannanBachemDomain(const Field &F) : _F(F), _MD(_F) {}
		SmithFormKannanBachemDomain(const SmithFormKannanBachemDomain &D) : _F(D._F), _MD(D._MD) {}

	private:
		template<typename Matrix>
		void swapRows(Matrix &M, size_t r1, size_t r2) const {
			SubMatrix row1(M, r1, 0, 1, M.coldim());
			SubMatrix row2(M, r2, 0, 1, M.coldim());
			
			row1.swap(row2);
		}
		
		template<typename Matrix>
		void swapCols(Matrix &M, size_t c1, size_t c2) const {
			SubMatrix col1(M, 0, c1, M.rowdim(), 1);
			SubMatrix col2(M, 0, c2, M.rowdim(), 1);
			
			col1.swap(col2);
		}
		
		void dxgcd(Element &s, Element &t, Element &u, Element &v, const Element &a, const Element &b) const {
			if (_F.isDivisor(b, a)) {
				_F.assign(s, _F.one);
				_F.assign(t, _F.zero);
				_F.assign(u, _F.one);
				_F.quo(v, b, a);
				return;
			}
			
			Element g;
			_F.dxgcd(g, s, t, u, v, a, b);
		}
		
		template<class Matrix>
		bool findPivot(Matrix &A) {
			for (size_t i = 0; i < A.rowdim(); i++) {
				for (size_t j = 0; j < A.coldim(); j++) {
					Element tmp;
					A.getEntry(tmp, i, j);
					
					if (!_F.isZero(tmp)) {
						if (i > 0) {
							swapRows(A, 0, i);
						}
						
						if (j > 0) {
							swapCols(A, 0, j);
						}
						
						return true;
					}
				}
			}
			
			return false;
		}
		
		template<class Matrix>
		void eliminateRow1(Matrix &A, size_t idx) {
			Element pivot, other;
			
			A.getEntry(other, 0, idx);
			
			if (_F.isZero(other)) {
				return;
			}
			
			A.getEntry(pivot, 0, 0);
			
			Element s, t, u, v;
			dxgcd(s,t,u,v,pivot,other);
			_F.negin(v);
			
			SubMatrix pivotCol(A, 0, 0, A.rowdim(), 1);
			SubMatrix otherCol(A, 0, idx, A.rowdim(), 1);
			
			OwnMatrix pivotCopy(pivotCol);
			
			_MD.mulin(pivotCol, s);
			_MD.saxpyin(pivotCol, t, otherCol);
			
			_MD.mulin(otherCol, u);
			_MD.saxpyin(otherCol, v, pivotCopy);
		}
		
		template<class Matrix>
		void eliminateRow(Matrix &A) {
			for (size_t i = 1; i < A.coldim(); i++) {
				eliminateRow1(A, i);
			}
		}
		
		template<class Matrix>
		void zeroRow(Matrix &A) {
			for (size_t i = 1; i < A.coldim(); i++) {
				A.setEntry(0, i, _F.zero);
			}
		}
		
		template<class Matrix>
		void eliminateCol1(Matrix &A, size_t idx) {
			Element pivot, other;
			
			A.getEntry(other, idx, 0);
			
			if (_F.isZero(other)) {
				return;
			}
			
			A.getEntry(pivot, 0, 0);
			
			Element s, t, u, v;
			dxgcd(s,t,u,v,pivot,other);
			_F.negin(v);
			
			SubMatrix pivotRow(A, 0, 0, 1, A.coldim());
			SubMatrix otherRow(A, idx, 0, 1, A.coldim());
			
			OwnMatrix pivotCopy(pivotRow);
			
			_MD.mulin(pivotRow, s);
			_MD.saxpyin(pivotRow, t, otherRow);
			
			_MD.mulin(otherRow, u);
			_MD.saxpyin(otherRow, v, pivotCopy);
		}
		
		template<class Matrix>
		void eliminateCol(Matrix &A) {
			for (size_t i = 1; i < A.rowdim(); i++) {
				eliminateCol1(A, i);
			}
		}
		
		template<class Matrix>
		bool isDiagonalized(Matrix &A) {
			Element tmp;
			for (size_t i = 1; i < A.rowdim(); i++) {
				if (!_F.isZero(A.getEntry(tmp, i, 0))) {
					return false;
				}
			}
			
			for (size_t i = 1; i < A.coldim(); i++) {
				if (!_F.isZero(A.getEntry(tmp, 0, i))) {
					return false;
				}
			}
			
			return true;
		}
		
		void fixDiagonal(std::vector<Element> &v) {
			for (size_t i = 0; i < v.size() - 1; i++) {
				if (_F.isZero(v[i+1])) {
					return;
				}
				
				Element h;
				_F.assign(h, v[i]);
				
				_F.gcd(v[i], v[i+1], h);
				
				_F.mulin(v[i+1], h);
				_F.divin(v[i+1], v[i]);
			}
		}
		
		template<class Matrix>
		void reduceOffDiagonal(Matrix &A) {
			size_t dim = A.rowdim() < A.coldim() ? A.rowdim() : A.coldim();
			
			for (size_t i = 1; i < dim; i++) {
				Element pivot, other;
				A.getEntry(pivot, i, i);
				A.getEntry(other, 0, i);
				
				if (_F.degree(other) < _F.degree(pivot)) {
					continue;
				}
				
				Element q;
				_F.quo(q, other, pivot);
				_F.negin(q);
				
				SubMatrix pivotRow(A, 0, i, 1, A.coldim() - 1);
				SubMatrix otherRow(A, i, i, 1, A.coldim() - 1);
				
				_MD.saxpyin(pivotRow, q, otherRow);
			}
		}
		
		template<class Matrix>
		void hermite(Matrix &A) {
			if (A.rowdim() == 0 || A.coldim() == 0) {
				return;
			}
			
			if (!findPivot(A)) {
				return;
			}
			
			eliminateCol(A);
			SubMatrix B(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			hermite(B);
			reduceOffDiagonal(A);
		}
		
		template<class Matrix>
		void solveHelper(std::vector<Element> &L, Matrix &A) {
			if (A.rowdim() == 0 || A.coldim() == 0) {
				return;
			}
			
			if (!findPivot(A)) {
				size_t dim = A.rowdim() < A.coldim() ? A.rowdim() : A.coldim();
				for (size_t i = 0; i < dim; i++) {
					L.push_back(_F.zero);
				}
				return;
			}
			
			while (!isDiagonalized(A)) {
				eliminateRow(A);
				hermite(A);
			}
			
			Element tmp;
			L.push_back(A.getEntry(tmp, 0, 0));
			SubMatrix B(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			solveHelper(L, B);
		}
		
		template<class Matrix>
		void solveTextBookHelper(std::vector<Element> &L, Matrix &A) {
			if (A.rowdim() == 0 || A.coldim() == 0) {
				return;
			}
			
			if (!findPivot(A)) {
				size_t dim = A.rowdim() < A.coldim() ? A.rowdim() : A.coldim();
				for (size_t i = 0; i < dim; i++) {
					L.push_back(_F.zero);
				}
				return;
			}
			
			while (!isDiagonalized(A)) {
				eliminateCol(A);
				
				Element pivot;
				if (_F.isUnit(A.getEntry(pivot, 0, 0))) {
					break;
				} else {
					eliminateRow(A);
				}
			}
			
			Element tmp;
			L.push_back(A.getEntry(tmp, 0, 0));
			SubMatrix B(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			solveTextBookHelper(L, B);
		}

	public:
		template<class Matrix>
		void solve(std::vector<Element> &L, Matrix &A) {
			solveHelper(L, A);
			fixDiagonal(L);
		}
		
		template<class Matrix>
		void solveTextBook(std::vector<Element> &L, Matrix &A) {
			solveTextBookHelper(L, A);
			fixDiagonal(L);
		}
		
		template<class Matrix>
		void halfSolve(std::vector<Element> &L, Matrix &A) {
			if (A.rowdim() == 0 || A.coldim() == 0) {
				return;
			}
			
			if (!findPivot(A)) {
				size_t dim = A.rowdim() < A.coldim() ? A.rowdim() : A.coldim();
				for (size_t i = 0; i < dim; i++) {
					L.push_back(_F.zero);
				}
				return;
			}
			
			eliminateCol(A);
			
			Element tmp;
			L.push_back(A.getEntry(tmp, 0, 0));
			SubMatrix B(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			halfSolve(L, B);
		}
	};
}

#endif // __LINBOX_smith_form_kannan_bachem_domain_H
