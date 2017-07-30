
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
		bool isZero(Vector &P) const {
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
		
		template<typename Vector>
		void vectorToPoly(Polynomial &p, const Vector &v) const {
			_PD.init(p, v);
		}
		
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
		void dxgcd(Polynomial &s, Polynomial &t, Polynomial &u, Polynomial &v, const Polynomial &a, const Polynomial &b) const {
			if (_PD.areEqual(a,b)) {
				_PD.assign(s, _PD.one);
				_PD.assign(t, _PD.zero);
				_PD.assign(u, _PD.one);
				_PD.assign(v, _PD.one);
				return;
			}

			Polynomial g;
			_PD.dxgcd(g,s,t,u,v,a,b);
		}
		
		template<typename PMatrix>
		void makeLeftElim(PMatrix &L, size_t pivotRow, size_t otherRow, const Polynomial &s, const Polynomial &t,
			Polynomial &u, Polynomial &v) {
			
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
		
		// use row pivotRow to eliminate row otherRow in the pivotRow-th column
		template<typename PMatrix>
		void eliminateCol(PMatrix &M, size_t pivotRow, size_t otherRow) {
			Polynomial other;
			_PD.init(other, M(otherRow, pivotRow));
			
			if (_PD.isZero(other)) {
				return;
			}
			
			Polynomial pivot;
			_PD.init(pivot, M(pivotRow, pivotRow));
			
			Polynomial s, t, u, v;
			dxgcd(s, t, u, v, pivot, other);
						
			PMatrix L(_F, M.rowdim(), M.rowdim(), size);
			makeLeftElim(L, pivotRow, otherRow, s, t, u, v);
			
			PMatrix Z(_F, M.coldim(), M.rowdim(), 1);
			_PMD.mul(Z, L, M);
			
			M.setsize(Z.real_degree() + 1);
			M.copy(Z, 0, Z.real_degree());
		}
		
		template<typename PMatrix>
		void eliminateCol(PMatrix &M, size_t pivotRow) {
			for (size_t otherRow = pivotRow + 1; otherRow < M.rowdim(); otherRow++) {
				eliminateCol(M, pivotRow, otherRow);
			}
		}
		
		template<typename PMatrix>
		void makeRightElim(PMatrix &R, size_t pivotCol, size_t otherCol, const Polynomial &s, const Polynomial &t,
			Polynomial &u, Polynomial &v) {
			
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
		void eliminateRow(PMatrix &M, size_t pivotCol, size_t otherCol) {
			Polynomial other;
			_PD.init(other, M(pivotCol, otherCol));
			
			if (_PD.isZero(other)) {
				return;
			}
			
			Polynomial pivot;
			_PD.init(pivot, M(pivotCol, pivotCol));
			
			Polynomial s, t, u, v;
			dxgcd(s, t, u, v, pivot, other);
			
			PMatrix R(_F, M.coldim(), M.coldim(), size);
			makeRightElim(R, pivotCol, otherCol, s, t, u, v);
			
			PMatrix Z(_F, M.coldim(), M.rowdim(), 1);
			_PMD.mul(Z, M, R);
			
			M.setsize(Z.real_degree() + 1);
			M.copy(Z, 0, Z.real_degree());
		}
		
		template<typename PMatrix>
		void eliminateRow(PMatrix &M, size_t pivotCol) {
			for (size_t otherCol = pivotCol + 1; otherCol < M.coldim(); otherCol++) {
				eliminateRow(M, pivotCol, otherCol);
			}
		}
		
		// Makes the bottom right n-by-n into a lower triangular hermite matrix.
		template<typename PMatrix>
		void hermite(PMatrix &M, size_t n) {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			for (size_t pivot = n; pivot < dim; pivot++) {
				if (!findPivot(M, pivot)) {
					break;
				}
				
				eliminateRow(M, pivot);
			}
		}

		// True if pivot row/col is zero for all indexes greater other than pivot
		template<typename PMatrix>
		bool isDiagonalized(const PMatrix &M, int pivot) const {
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
		
	public:
		
		template<typename PMatrix>
		void solve(PMatrix &M) {
			eliminateCol(M, 0);
			
			hermite(M, 0);
		}
		
		template<typename PMatrix>
		void solveTextbook(PMatrix &M) {
			size_t dim = M.rowdim() < M.coldim() ? M.rowdim() : M.coldim();
			
			for (size_t pivot = 0; pivot < dim; pivot++) {
				if(!findPivot(M, pivot)) {
					return;
				}
				
				while(true) {
					if (isDiagonalized(M, pivot)) {
						break;
					}
					
					eliminateCol(M, pivot);
					eliminateRow(M, pivot);
				}
			}
		}
		
		
/*
		void reduceOffDiagonal(Rep &A, int s, int e) const
		{
			for (int i = s; i <= e; i++)
			{
				Element nii,ii;
				A.getEntry(ii, i, i);
				field().normalize(nii,ii);

				Element tmp;
				field().div(tmp, nii, ii);

				if (field().isOne(tmp))
					continue;

				// A[i] = A[i,i] * n
				// where n = normalized(A[i,i]) / A[i,i]
				A.setEntry(i, i, nii);
				for (size_t j = i+1; j < A.coldim(); j++)
				{
					Element ij;
					A.getEntry(ij, i, j);
					field().mulin(ij, tmp);
					A.setEntry(i, j, ij);
				}
			}

			// Ording of reduction here is an improvement to Kannan/Bachem
			// Introduced by Chou/Collins '82
			// Reduce from bottom to top and left to right
			// * 4 5 6
			// 0 * 2 3
			// 0 0 * 1
			// 0 0 0 *
			for (int i = e-1; i >= s; i--)
			{
				for (int j = i+1; j <= e; j++)
				{
					Element jj, ij, tmp;

					A.getEntry(ij, i, j);
					if (field().isZero(ij))
						continue;

					A.getEntry(jj, j, j);
					field().quo(tmp, ij, jj);

					// A[i] = A[i] - quo(A[i,j], A[j,j]) * A[j]
					for (size_t k = j; k < A.coldim(); k++)
					{
						Element ik, jk;

						A.getEntry(ik, i, k);
						A.getEntry(jk, j, k);

						field().mulin(jk, tmp);
						field().subin(ik, jk);

						A.setEntry(i, k, ik);
					}
				}
			}
		}

		// Puts the lower-right n-by-n minor of A into Hermite Normal Form
		void hermite(Rep &A, int n) const
		{
			int dim = (int)A.rowdim();

			for (int i = n; i < dim; i++)
			{
				for (int j = n; j < i; j++)
					eliminateCol(A, j, i);

				if (!findPivot(A, i))
					return;

				reduceOffDiagonal(A, n, i);
			}
		}

		bool isRowDiagonalized(const Rep &A, int n) const
		{
			for (size_t i = n+1; i < A.coldim(); i++)
			{
				Element ni;
				A.getEntry(ni, n, i);
				if (!field().isZero(ni))
					return false;
			}
			return true;
		}

		bool pivotDividesRemaining(Rep &A, int n) const
		{
			Element nn;
			A.getEntry(nn, n, n);

			for (size_t i = n+1; i < A.rowdim(); i++)
			{
				for (size_t j = i; j < A.coldim(); j++)
				{
					Element ij, g;
					A.getEntry(ij, i, j);

					if (field().isZero(ij))
						continue;

					field().gcd(g, nn, ij);

					if (!field().areAssociates(g, nn))
					{
						// Add row i to row n
						for (size_t k = i; k < A.coldim(); k++)
						{
							Element ik;
							A.getEntry(ik, i, k);
							A.setEntry(n, k, ik);
						}
						return false;
					}
				}
			}

			return true;
		}

	public:
		template<class Vector>
		Vector &solve(Vector &S, const Rep &A) const
		{
			size_t dim = A.rowdim();
			linbox_check(A.coldim() == dim && S.size() >= dim);

			Rep B(A);

			for (size_t i = 0; i < dim;)
			{
				if (!findPivot(B, i))
					break;

				for (size_t j = i+1; j < dim; j++)
					eliminateRow(B, i, j);

				hermite(B, i);

				if (!isRowDiagonalized(B, i))
					continue;

				if (!pivotDividesRemaining(B, i))
					continue;

				i++;
				
				//std::cout << i << "/" << dim << std::endl;
			}

			for (size_t i = 0; i < dim; i++)
			{
				Element ii;
				B.getEntry(ii, i, i);
				S.setEntry(i, ii);
			}

			return S;
		}
	*/
	};
}

#endif // __LINBOX_smith_form_kannan_bachem_domain_H
