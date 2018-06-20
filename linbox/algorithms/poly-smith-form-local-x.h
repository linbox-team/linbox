
/* tests/__LINBOX_smith_form_kannan_bachem.h
 * Copyright (C) 2017 Gavin Harrison,
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
#include <stdlib.h>
#include <algorithm>

#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"
#include "linbox/ring/polynomial-local-x.h"

#ifndef __LINBOX_poly_smith_form_local_x_domain_H
#define __LINBOX_poly_smith_form_local_x_domain_H

namespace LinBox
{
	// Specialized Polynomial Smith Form over Polynomials mod x^n
	// BaseField should be a Polynomial Ring
	template<class BaseField>
	class PolySmithFormLocalXDomain
	{
	public:
		typedef typename BaseField::Element Polynomial;
		typedef typename BaseField::Coeff Coeff;
		
		typedef PolynomialLocalX<BaseField> Field;
		typedef typename Field::Element Element;
		
		typedef MatrixDomain<Field> MatrixDom;
		
		typedef typename MatrixDom::OwnMatrix OwnMatrix;
		typedef typename MatrixDom::Matrix SubMatrix;

	private:
		BaseField _BF;
		Field _F;
		MatrixDom _MD;

	public:
		PolySmithFormLocalXDomain(const BaseField &F, size_t e) : _BF(F), _F(F, e), _MD(_F) {}
		PolySmithFormLocalXDomain(const PolySmithFormLocalXDomain &D) : _BF(D._BF), _F(D._F), _MD(_F) {}

	private:
		
		template<class Matrix1>
		void printMatrix(const Matrix1 &A) const {
			std::cout << "[" << std::endl;
			for (size_t i = 0; i < A.rowdim(); i++) {
				std::cout << "\t[";
				for (size_t j = 0; j < A.coldim(); j++) {
					_F.write(std::cout, A.getEntry(i, j));
					if (j < A.coldim() - 1) {
						std::cout << ", ";
					}
				}
				std::cout << "]";
				if (i < A.rowdim() - 1) {
					std::cout << ",";
				}
				std::cout << std::endl;
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
		
		template<class Matrix>
		bool findPivot(Matrix &A) const {
			int d_row = -1, d_col = -1; // location of divisor
			size_t d_firstNonZero = -1;
			Element divisor;
			bool d_isUnit = false;
			
			for (size_t i = 0; i < A.rowdim() && !d_isUnit; i++) {
				for (size_t j = 0; j < A.coldim() && !d_isUnit; j++) {
					Element tmp;
					_F.assign(tmp, A.getEntry(i, j));
					
					if (_F.isZero(tmp)) {
						continue;
					} 
					
					size_t firstNonZero = _F.firstNonZeroCoeff(tmp);
					if (firstNonZero == 0) {
						d_row = i;
						d_col = j;
						d_isUnit = true;
						continue;
					} 
					
					if (d_row > -1 && d_col > -1 && d_firstNonZero <= firstNonZero) {
						continue;
					}
					
					_F.assign(divisor, tmp);
					d_row = i;
					d_col = j;
					d_firstNonZero = firstNonZero;
				}
			}
			
			if (d_row == -1 && d_col == -1) {
				return false;
			}
			
			swapRows(A, 0, d_row);
			swapCols(A, 0, d_col);
			return true;
		}
		
		template<class Matrix>
		size_t reduceMatrix(Matrix &A) {
			Element pivot;
			A.getEntry(pivot, 0, 0);
			
			if (_F.isUnit(pivot)) {
				return 0;
			}
			
			size_t e = _F.firstNonZeroCoeff(pivot);
			
			for (size_t i = 0; i < A.rowdim(); i++) {
				for (size_t j = 0; j < A.coldim(); j++) {
					Element tmp;
					A.getEntry(tmp, i, j);
					_F.rightShiftIn(tmp, e);
					A.setEntry(i, j, tmp);
				}
			}
			
			_F.setExponent(_F.getExponent() - e);
			
			return e;
		}
		
		template<class Matrix>
		void zeroRow(Matrix &A) const {
			for (size_t i = 1; i < A.coldim(); i++) {
				A.setEntry(0, i, _F.zero);
			}
		}
		
		template<class Matrix>
		void eliminateCol(Element &pivot, Matrix &A) const {
			Element pivot_inv;
			A.getEntry(pivot, 0, 0);
			
			_F.inv(pivot_inv, pivot);
			_F.negin(pivot_inv);
			
			SubMatrix pivotRow(A, 0, 0, 1, A.coldim());
			_MD.mulin(pivotRow, pivot_inv);
			
			for (size_t r = 1; r < A.rowdim(); r++) {
				SubMatrix otherRow(A, r, 0, 1, A.coldim());
				
				Element other;
				otherRow.getEntry(other, 0, 0);
				
				_MD.saxpyin(otherRow, other, pivotRow);
			}
		}
		
		template<class Matrix>
		void solveHelper(std::vector<Element> &L, std::vector<size_t> &es, size_t e, Matrix &A) {
			if (A.rowdim() == 0 || A.coldim() == 0) {
				return;
			}
			
			if (!findPivot(A)) {
				size_t dim = std::min(A.rowdim(), A.coldim());
				for (size_t i = 0; i < dim; i++) {
					L.push_back(_F.zero);
				}
				return;
			}
			
			e += reduceMatrix(A);
			
			Element pivot;
			eliminateCol(pivot, A);
			
			zeroRow(A);
			
			L.push_back(pivot);
			es.push_back(e);
			SubMatrix B(A, 1, 1, A.rowdim() - 1, A.coldim() - 1);
			solveHelper(L, es, e, B);
		}

	public:
		template<class Matrix>
		void solve(std::vector<Polynomial> &result, Matrix &A) {
			size_t initial_exp = _F.getExponent();
			
			std::vector<Element> L;
			std::vector<size_t> es;
			
			solveHelper(L, es, 0, A);
			
			_F.setExponent(initial_exp);
			
			for (size_t i = 0; i < L.size(); i++) {
				Element tmp;
				_F.leftShift(tmp, L[i], es[i]);
								
				Polynomial r;
				_F.denormalize(r, tmp);
				result.push_back(r);
			}
		}
		
		template<class Matrix>
		void solveDet(Polynomial &det, Matrix &A) {
			size_t initial_exp = _F.getExponent();
			
			std::vector<Element> L;
			std::vector<size_t> es;
			solveHelper(L, es, 0, A);
			
			_F.setExponent(initial_exp);
			
			size_t e = 0;
			for (size_t i = 0; i < es.size(); i++) {
				e += es[i];
			}
			
			Element tmp;
			_F.leftShift(tmp, L[0], e);
			for (size_t i = 1; i < L.size(); i++) {
				_F.mulin(tmp, L[i]);
			}
			_F.denormalize(det, tmp);
		}
	}; // end of class PolySmithFormLocalXDomain
}

#endif // __LINBOX_poly_smith_form_local_x_domain_H
