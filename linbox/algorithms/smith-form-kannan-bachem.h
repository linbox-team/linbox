
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
		void printMatrix(const Matrix &A) const {
			std::cout << "[" << std::endl;
			for (size_t i = 0; i < A.rowdim(); i++) {
				std::cout << "\t[";
				for (size_t j = 0; j < A.coldim(); j++) {
					_F.write(std::cout, A.getEntry(i, j)) << ", ";
				}
				std::cout << "]" << std::endl;
			}
			std::cout << "]" << std::endl;
		}
		
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
					if (_F.isZero(A.getEntry(i, j))) {
						continue;
					}
					
					if (i > 0) {
						swapRows(A, 0, i);
					}
					
					if (j > 0) {
						swapCols(A, 0, j);
					}
					
					return true;
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
			for (size_t i = 1; i < A.rowdim(); i++) {
				if (!_F.isZero(A.getEntry(i, 0))) {
					return false;
				}
			}
			
			for (size_t i = 1; i < A.coldim(); i++) {
				if (!_F.isZero(A.getEntry(0, i))) {
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
				
				if (_F.isZero(other)) {
					continue;
				}
				
				Element q;
				_F.quo(q, other, pivot);
				
				if (_F.isZero(q)) {
					continue;
				}
				
				_F.negin(q);
				
				SubMatrix pivotRow(A, i, i, 1, A.coldim() - i);
				SubMatrix otherRow(A, 0, i, 1, A.coldim() - i);
				
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
			
			L.push_back(A.getEntry(0, 0));
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
				
				if (_F.isUnit(A.getEntry(0, 0))) {
					break;
				}
				
				eliminateRow(A);
			}
			
			L.push_back(A.getEntry(0, 0));
			SubMatrix B(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			solveTextBookHelper(L, B);
		}
		
		// Iliopoulos Specific Methods
		
		template<class Matrix>
		bool findZeroCol(Matrix &A) {
			for (size_t i = 1; i < A.coldim(); i++) {
				if (_F.isZero(A.getEntry(0, i))) {
					swapCols(A, 0, i);
					return true;
				}
			}
			
			return false;
		}
		
		template<class Matrix>
		bool findZeroRow(Matrix &A) {
			for (size_t i = 1; i < A.rowdim(); i++) {
				if (_F.isZero(A.getEntry(i, 0))) {
					swapRows(A, 0, i);
					return true;
				}
			}
			
			return false;
		}
		
		void xgcd(Element &g, std::vector<Element> &ts, const std::vector<Element> &as, const Element &d) {
			ts.resize(as.size());
			_F.assign(ts[0], _F.one);
			_F.assign(g, as[0]);
			
			for (size_t i = 1; i < as.size(); i++) {
				Element r, t1, t2;
				_F.gcd(r, t1, t2, g, as[i]);
				
				for (size_t j = 0; j < i; j++) {
					_F.mulin(ts[j], t1);
				}
				
				_F.assign(ts[i], t2);
				_F.assign(g, r);
			}
			
			Element r, t1, t2;
			_F.gcd(r, t1, t2, g, d);
			
			for (size_t i = 0; i < ts.size(); i++) {
				_F.mulin(ts[i], t1);
			}
			
			_F.assign(g, r);
		}
		
		template<class Matrix>
		void reduceMatrix(Matrix &A, const Element &d) {
			for (size_t i = 0; i < A.rowdim(); i++) {
				for (size_t j = 0; j < A.coldim(); j++) {
					Element tmp;
					A.getEntry(tmp, i, j);
					_F.modin(tmp, d);
					A.setEntry(i, j, tmp);
				}
			}
		}
		
		template<class Matrix>
		void makePivotCol(Matrix &A, const Element &d) {
			std::vector<Element> elms;
			for (size_t i = 1; i < A.coldim(); i++) {
				elms.push_back(A.getEntry(0, i));
			}
			
			Element g;
			std::vector<Element> ts;
			xgcd(g, ts, elms, d);
			
			SubMatrix pivotCol(A, 0, 0, A.rowdim(), 1);
			for (size_t i = 1; i < A.coldim(); i++) {
				if (_F.isZero(ts[i - 1])) {
					continue;
				}
				
				SubMatrix otherCol(A, 0, i, A.rowdim(), 1);
				
				_MD.saxpyin(pivotCol, ts[i - 1], otherCol);
			}
			
			A.setEntry(0, 0, g);
		}
		
		template<class Matrix>
		void eliminateRow(Matrix &A, const Element &d) {
			if (A.coldim() == 2) {
				eliminateRow1(A, 1);
				reduceMatrix(A, d);
				return;
			}
			
			if (!findZeroCol(A)) {
				swapCols(A, 0, 1);
				eliminateRow1(A, 1);
				swapCols(A, 0, 1);
			}
			
			makePivotCol(A, d);
			
			Element pivot;
			SubMatrix pivotCol(A, 0, 0, A.rowdim(), 1);
			pivotCol.getEntry(pivot, 0, 0);
			for (size_t i = 1; i < A.coldim(); i++) {
				Element other;
				SubMatrix otherCol(A, 0, i, A.rowdim(), 1);
				otherCol.getEntry(other, 0, 0);
				
				Element q;
				_F.quo(q, other, pivot);
				_F.negin(q);
				
				_MD.saxpyin(otherCol, q, pivotCol);
			}
			
			reduceMatrix(A, d);
		}
		
		template<class Matrix>
		void makePivotRow(Matrix &A, const Element &d) {
			std::vector<Element> elms;
			for (size_t i = 1; i < A.rowdim(); i++) {
				elms.push_back(A.getEntry(i, 0));
			}
			
			Element g;
			std::vector<Element> ts;
			xgcd(g, ts, elms, d);
			
			SubMatrix pivotRow(A, 0, 0, 1, A.coldim());
			for (size_t i = 1; i < A.rowdim(); i++) {
				if (_F.isZero(ts[i - 1])) {
					continue;
				}
				
				SubMatrix otherRow(A, i, 0, 1, A.coldim());
				
				_MD.saxpyin(pivotRow, ts[i - 1], otherRow);
			}
			
			A.setEntry(0, 0, g);
		}
		
		template<class Matrix>
		void eliminateCol(Matrix &A, const Element &d) {
			if (A.rowdim() == 2) {
				eliminateCol1(A, 1);
				reduceMatrix(A, d);
				return;
			}
			
			if (!findZeroRow(A)) {
				swapRows(A, 0, 1);
				eliminateRow1(A, 1);
				swapRows(A, 0, 1);
			}
			
			makePivotRow(A, d);
			
			Element pivot;
			SubMatrix pivotRow(A, 0, 0, 1, A.coldim());
			pivotRow.getEntry(pivot, 0, 0);
			for (size_t i = 1; i < A.rowdim(); i++) {
				Element other;
				SubMatrix otherRow(A, i, 0, 1, A.coldim());
				otherRow.getEntry(other, 0, 0);
				
				Element q;
				_F.quo(q, other, pivot);
				_F.negin(q);
				
				_MD.saxpyin(otherRow, q, pivotRow);
			}
			
			reduceMatrix(A, d);
		}
		
		template<class Matrix>
		void solveIliopoulosHelper(std::vector<Element> &L, Matrix &A, const Element &d) {			
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
				eliminateCol(A, d);
				
				if (_F.isUnit(A.getEntry(0, 0))) {
					break;
				}
				
				eliminateRow(A, d);
			}
			
			L.push_back(A.getEntry(0, 0));
			SubMatrix B(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			solveIliopoulosHelper(L, B, d);
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
			
			L.push_back(A.getEntry(0, 0));
			SubMatrix B(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			halfSolve(L, B);
		}
		
		template<class Matrix>
		void solveIliopoulos(std::vector<Element> &L, Matrix &A, const Element &d) {
			solveIliopoulosHelper(L, A, d);
			fixDiagonal(L);
		}
		
		template<class Matrix>
		void solveAdaptive(std::vector<Element> &L, Matrix &A) {
			std::vector<Element> ds;
			halfSolve(ds, A);
			
			Element d;
			_F.assign(d, ds[0]);
			for (size_t i = 1; i < ds.size(); i++) {
				_F.mulin(d, ds[i]);
			}
			
			eliminateRow(A, d);
			
			solveIliopoulosHelper(L, A, d);
			fixDiagonal(L);
		}
	};
}

#endif // __LINBOX_smith_form_kannan_bachem_domain_H
