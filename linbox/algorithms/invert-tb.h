
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
#include "linbox/matrix/densematrix/blas-matrix.h"
#include "linbox/matrix/matrixdomain/matrix-domain.h"

#ifndef __LINBOX_invert_textbook_domain_H
#define __LINBOX_invert_textbook_domain_H

namespace LinBox
{
	/**
	 * Assumes that Field is a field, not a ring.
	 * Uses Gauss-Jordan elimination.
	 */
	template<class Field>
	class InvertTextbookDomain
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
		InvertTextbookDomain(const Field &F) : _F(F), _MD(_F) {}
		InvertTextbookDomain(const InvertTextbookDomain &D) : _F(D._F), _MD(D._MD) {}

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
		
		template<class Matrix>
		bool findPivot(Matrix &A, size_t pivot) {
			for (size_t i = pivot; i < A.rowdim(); i++) {
				if (_F.isZero(A.getEntry(i, pivot))) {
					continue;
				}
				
				if (i != pivot) {
					swapRows(A, pivot, i);
				}
				
				return true;
			}
			
			return false;
		}
		
		template<class Matrix>
		bool eliminateCol(Matrix &A, size_t pivotIdx) {
			SubMatrix pivot(A, pivotIdx, pivotIdx, 1, A.coldim() - pivotIdx);
			
			if (!_F.isUnit(pivot.getEntry(0, 0))) {
				return false;
			}
			
			Element inv;
			_F.inv(inv, pivot.getEntry(0, 0));
			_MD.mulin(pivot, inv);
			
			for (size_t i = 0; i < A.rowdim(); i++) {
				if (i == pivotIdx) {
					continue;
				}
				
				SubMatrix other(A, i, pivotIdx, 1, A.coldim() - pivotIdx);
				
				Element q;
				_F.neg(q, other.getEntry(0, 0));
				
				_MD.saxpyin(other, q, pivot);
			}
			
			return true;
		}

	public:
		
		// Returns false if invert fails, i.e., not full rank or element fails to invert
		template<class Matrix1, class Matrix2>
		bool invert(Matrix1 &Ainv, Matrix2 &A) {	
			size_t dim = A.rowdim();
			if (dim != A.coldim() || dim != Ainv.rowdim() || dim != Ainv.coldim()) {
				return false;
			}
			
			OwnMatrix B(_F, dim, 2 * dim);
			SubMatrix Acopy(B, 0, 0, dim, dim);
			_MD.copy(Acopy, A);
			
			SubMatrix Inv(B, 0, dim, dim, dim);
			Inv.zero();
			for (size_t i = 0; i < dim; i++) {
				Inv.setEntry(i, i, _F.one);
			}
			
			for (size_t i = 0; i < dim; i++) {				
				if (!findPivot(B, i)) {
					return false;
				}
				
				if (!eliminateCol(B, i)) {
					return false;
				}
			}
			
			_MD.copy(Ainv, Inv);
			return true;
		}
	};
}

#endif // __LINBOX_invert_textbook_domain_H
