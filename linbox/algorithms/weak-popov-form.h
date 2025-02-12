/* tests/__LINBOX_smith_form_kannan_bachem.h
 * Copyright (C) 2018 Gavin Harrison,
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

#ifndef __LINBOX_weak_popov_form_domain_H
#define __LINBOX_weak_popov_form_domain_H

#include <iostream>
#include <stdlib.h>
#include <algorithm>
#include <utility>
#include <vector>

#include "linbox/matrix/dense-matrix.h"

namespace LinBox
{
	// PolynomialRing = K[x] where K is a field
	template<class PolynomialRing>
	class WeakPopovFormDomain
	{
	public:
		typedef typename PolynomialRing::Element Polynomial;
		typedef typename PolynomialRing::Coeff Coeff;
		typedef typename PolynomialRing::CoeffField CoeffField;

	private:
		PolynomialRing _R;

	public:
	  WeakPopovFormDomain(const PolynomialRing &R) : _R(R) {}
	  WeakPopovFormDomain(const WeakPopovFormDomain &D) : _R(D._R)  {}

	// private:

		template<class Matrix1>
		void printMatrix(const Matrix1 &A) const {
			std::cout << "[" << std::endl;
			for (size_t i = 0; i < A.rowdim(); i++) {
				std::cout << "\t[";
				for (size_t j = 0; j < A.coldim(); j++) {
					_R.write(std::cout, A.getEntry(i, j));
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
		  // SubMatrix row1(M, r1, 0, 1, M.coldim());
		  // SubMatrix row2(M, r2, 0, 1, M.coldim());
		  // row1.swap(row2);
		  // -> use row iterator no submatrix
		  std::swap(*(M.rowBegin()+r1), *(M.rowBegin()+r2));
		}

		template<typename Matrix>
		void swapCols(Matrix &M, size_t c1, size_t c2) const {
		  // SubMatrix col1(M, 0, c1, M.rowdim(), 1);
		  // SubMatrix col2(M, 0, c2, M.rowdim(), 1);
		  // col1.swap(col2);
		  // -> use row iterator no submatrix
		  std::swap(*(M.colBegin()+c1), *(M.colBegin()+c2));
		}

		/**
		 * Gets a list of weak Popov pivots: indexes of the right-most max degree entry of the row.
		 * returns -1 if the row is zero, or the index of the pivot.
		 */
		template<typename Matrix>
		void findPivots(std::vector<long> &pivots, const Matrix &M) const {
			pivots.clear();

			for (size_t i = 0; i < M.rowdim(); i++) {
				long index = -1;
				size_t max_deg = 0;

				for (size_t j = 0; j < M.coldim(); j++) {
					Polynomial tmp;
					M.getEntry(tmp, i, j);

					size_t deg = _R.deg(tmp);
					if (_R.isZero(tmp) || deg < max_deg) {
						continue;
					}

					index = j;
					max_deg = deg;
				}

				pivots.push_back(index);
			}
		}

		bool findMatchingPivots(std::pair<size_t, size_t> &pair, const std::vector<long> &pivots) const {
			for (size_t i = 0; i < pivots.size() - 1; i++) {
				if (pivots[i] < 0) {
					continue;
				}

				for (size_t j = i + 1; j < pivots.size(); j++) {
					if (pivots[i] == pivots[j]) {
						std::pair<size_t, size_t> rv(i, j);
						pair = rv;
						return true;
					}
				}
			}

			return false;
		}

		template<typename Matrix>
		void eliminate(const Matrix &M1, Matrix &M2, const Matrix &V1, Matrix &V2, const Coeff &c, size_t e) const {
			for (size_t i = 0; i < M1.coldim(); i++) {
				Polynomial f1, f2;

				M1.getEntry(f1, 0, i);
				M2.getEntry(f2, 0, i);

				_R.mulCoeffIn(f1, c);
				_R.leftShiftIn(f1, e);

				_R.addin(f2, f1);

				M2.setEntry(0, i, f2);
			}

			for (size_t i = 0; i < V1.coldim(); i++) {
				Polynomial f1, f2;

				V1.getEntry(f1, 0, i);
				V2.getEntry(f2, 0, i);

				_R.mulCoeffIn(f1, c);
				_R.leftShiftIn(f1, e);

				_R.addin(f2, f1);

				V2.setEntry(0, i, f2);
			}
		}

		template<typename Matrix1, typename Matrix2>
		void eliminate(Matrix1 &M, Matrix2 &V, size_t row1, size_t row2, size_t pivot_col) const {
			Polynomial p1, p2;
			M.getEntry(p1, row1, pivot_col);
			M.getEntry(p2, row2, pivot_col);

			using SubMatrix1=DenseSubmatrix<typename Matrix1::Field>;
			using SubMatrix2=DenseSubmatrix<typename Matrix2::Field>;

			SubMatrix1 M1(M, row1, 0, 1, M.coldim());
			SubMatrix1 M2(M, row2, 0, 1, M.coldim());
			SubMatrix2 V1(V, row1, 0, 1, V.coldim());
			SubMatrix2 V2(V, row2, 0, 1, V.coldim());

			size_t d1 = _R.deg(p1);
			size_t d2 = _R.deg(p2);

			Coeff c1, c2;
			_R.leadCoeff(c1, p1);
			_R.leadCoeff(c2, p2);

			if (d1 < d2) {
				size_t e = d2 - d1;
				Coeff tmp;
				_R.getCoeffField().div(tmp, c2, c1);
				_R.getCoeffField().negin(tmp);

				eliminate(M1, M2, V1, V2, tmp, e);
			} else {
				size_t e = d1 - d2;
				Coeff tmp;
				_R.getCoeffField().div(tmp, c1, c2);
				_R.getCoeffField().negin(tmp);

				eliminate(M2, M1, V2, V1, tmp, e);
			}
		}

		template<typename Matrix1, typename Matrix2>
		void extendedSolve(Matrix1 &M, Matrix2 &V) const {
			std::vector<long> pivots;
			findPivots(pivots, M);

			std::pair<size_t, size_t> match;
			bool matchFound = findMatchingPivots(match, pivots);

			while (matchFound) {
				size_t row1 = match.first;
				size_t row2 = match.second;
				size_t pivot_col = (size_t) pivots[row1];

				eliminate(M, V, row1, row2, pivot_col);

				findPivots(pivots, M);
				matchFound = findMatchingPivots(match, pivots);
			}
		}

		template<typename Matrix>
		size_t findZeroRow(const Matrix &M) const {
			for (size_t i = 0; i < M.rowdim(); i++) {
				bool isZero = true;

				for (size_t j = 0; j < M.coldim() && isZero; j++) {
					if (!_R.isZero(M.getEntry(i, j))) {
						isZero = false;
					}
				}

				if (isZero) {
					return i;
				}
			}

			return -1;
		}

		template<typename Matrix>
		void solveDetHelper(Polynomial &det, Matrix &T) const {
			if (T.rowdim() == 1) {
				_R.mulin(det, T.getEntry(0, 0));
				return;
			}

			using SubMatrix=DenseSubmatrix<typename Matrix::Field>;

			SubMatrix M(T, 0, 0, T.rowdim(), T.coldim() - 1);
			SubMatrix V(T, 0, T.coldim() - 1, T.rowdim(), 1);

			extendedSolve(M, V);

			size_t k = findZeroRow(M);

			Polynomial tmp;
			V.getEntry(tmp, k, 0);

			_R.mulin(det, tmp);

			swapRows(T, k, T.rowdim() - 1);

			SubMatrix SubT(T, 0, 0, T.rowdim() - 1, T.coldim() - 1);
			solveDetHelper(det, SubT);
		}

		template<typename Matrix>
		void solveDet(Polynomial &det, const Matrix &T_in) const {
			DenseMatrix<typename Matrix::Field> T(T_in);
			_R.assign(det, _R.one);

			solveDetHelper(det, T);
		}
	}; // end of class WeakPopovFormDomain
}

#endif // __LINBOX_weak_popov_form_domain_H
