
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

#include "linbox/ring/polynomial-local-x.h"

#include "linbox/matrix/matrixdomain/matrix-domain.h"

#include "linbox/algorithms/smith-form-local.h"
#include "linbox/algorithms/poly-smith-form-local-x.h"
#include "linbox/algorithms/weak-popov-form.h"
#include "linbox/algorithms/poly-dixon.h"

#ifndef __LINBOX_poly_smith_form_domain_H
#define __LINBOX_poly_smith_form_domain_H

namespace LinBox
{
	template<class Ring>
	class PolySmithFormDomain
	{
	public:
		typedef typename Ring::Element Polynomial;
		typedef typename Ring::Coeff Coeff;
		typedef typename Ring::RandIter RandIter;
		
		typedef MatrixDomain<Ring> MatrixDom;
		
		typedef typename MatrixDom::OwnMatrix Matrix;
		typedef typename MatrixDom::Matrix SubMatrix;

	private:
		Ring _R;
		RandIter _RI;
		MatrixDom _MD;

	public:
		PolySmithFormDomain(const Ring &R) : _R(R), _RI(R), _MD(_R) {}
		PolySmithFormDomain(const PolySmithFormDomain &D) : _R(D._R), _RI(D._RI), _MD(D._MD) {}
		
		template<class Matrix1>
		void local(
			std::vector<Polynomial> &result,
			const Matrix1 &M,
			const Polynomial &f,
			long multiplicity) const {
		
			typedef typename Ring::QuotientRing QuotientRing;
			typedef typename QuotientRing::Element QPolynomial;
		
			SmithFormLocal<QuotientRing> SFD;
			
			Polynomial modulus;
			_R.pow(modulus, f, multiplicity);
			
			QuotientRing QR(_R.quotient(modulus));
			
			typename MatrixDomain<QuotientRing>::OwnMatrix QM(M, QR);
			
			std::list<QPolynomial> L;
			SFD(L, QM, QR);
			
			Hom<Ring, QuotientRing> hom(_R, QR);
			
			size_t j = 0;
			typename std::list<QPolynomial>::const_iterator it;
			for (it = L.begin(); it != L.end(); it++) {
				Polynomial tmp;
				hom.preimage(tmp, *it);
									
				if (_R.isOne(tmp)) {
					// noop
				} else if (_R.isZero(tmp)) {
					_R.mulin(result[j], modulus);
				} else {
					_R.mulin(result[j], tmp);
				}
				j++;
			}
		}
		
		template<class Matrix1>
		void localRank(
			std::vector<Polynomial> &result,
			const Matrix1 &M,
			const Polynomial &f) const {
		
			typedef typename Ring::QuotientRing QuotientRing;
			typedef typename QuotientRing::Element QPolynomial;
			
			SmithFormLocal<QuotientRing> SFD;
			
			QuotientRing QR(_R.quotient(f));
			
			typename MatrixDomain<QuotientRing>::OwnMatrix QM(M, QR);
			
			std::list<QPolynomial> L;
			SFD(L, QM, QR);
			
			typename std::list<QPolynomial>::reverse_iterator it = L.rbegin();
			it++;
			if (QR.isZero(*it)) {
				_R.mulin(result[result.size() - 2], f);
				_R.mulin(result[result.size() - 1], f);
			} else {
				Polynomial f2;
				_R.pow(f2, f, 2);
				_R.mulin(result[result.size() - 1], f2);
			}
		}
		
		template<class Matrix1>
		void localFactored(
			std::vector<Polynomial> &result,
			const Matrix1 &M,
			const Polynomial &sf_factor,
			long multiplicity) const {
		
			std::vector<std::pair<Polynomial, long>> factors;
			_R.factor(factors, sf_factor);
			
			for (size_t i = 0; i < factors.size(); i++) {
				local(result, M, factors[i].first, factors[i].second * multiplicity);
			}
		}
		
		template<class Matrix1>
		void localFactoredRank(
			std::vector<Polynomial> &result,
			const Matrix1 &M,
			const Polynomial &sf_factor) const {
		
			std::vector<std::pair<Polynomial, long>> factors;
			_R.factor(factors, sf_factor);
			
			for (size_t i = 0; i < factors.size(); i++) {
				localRank(result, M, factors[i].first);
			}
		}
		
		template<class Matrix1>
		void solve(
			std::vector<Polynomial> &result,
			const Matrix1 &M,
			const Polynomial &det,
			bool isDet = false) const {
			
			std::vector<std::pair<Polynomial, long>> factors;
			
			_R.squareFree(factors, det);
			
			result.clear();
			for (size_t i = 0; i < M.rowdim(); i++) {
				result.push_back(_R.one);
			}
			
			for (size_t i = 0; i < factors.size(); i++) {
				if (isDet && factors[i].second == 1) {
					_R.mulin(result[result.size() - 1], factors[i].first);
				} else if (isDet && factors[i].second == 2) {
					localFactoredRank(result, M, factors[i].first);
				} else {
					localFactored(result, M, factors[i].first, factors[i].second);
				}
			}
		}
		
		// Methods for Computing a Modulus
		template<class Iterator>
		size_t detLimit(const Iterator &begin, const Iterator &end) const {
			size_t lim = 0;
			for (auto it = begin; it != end; it++) {
				size_t max_deg = 0;
				for (auto elm = it->begin(); elm != it->end(); elm++) {
					max_deg = std::max(max_deg, _R.deg(*elm));
				}
				lim += max_deg;
			}
			
			return lim;
		}
		
		template<class Matrix1>
		size_t detLimit(const Matrix1 &M) const {
			size_t limit1 = detLimit(M.rowBegin(), M.rowEnd());
			size_t limit2 = detLimit(M.colBegin(), M.colEnd());
			
			return std::min(limit1, limit2) + 1;
		}
		
		template<class Matrix1>
		void detLocalX(Polynomial &det, const Matrix1 &M) const {
			size_t exponent = detLimit(M);
		
			PolySmithFormLocalXDomain<Ring> SFD(_R, exponent);
			PolynomialLocalX<Ring> L(_R, exponent);
				
			typename MatrixDomain<PolynomialLocalX<Ring>>::OwnMatrix G(M, L);
			
			SFD.solveDet(det, G);
			_R.monicIn(det);
		}
		
		template<class Matrix1>
		void detPopov(Polynomial &det, const Matrix1 &M) const {
			WeakPopovFormDomain<Ring> PFD(_R);
			PFD.solveDet(det, M);
			_R.monicIn(det);
		}
		
		template<class Matrix1>
		bool dixon(
			Polynomial &minpoly,
			const Matrix1 &M,
			const Polynomial &f,
			size_t max_deg) const {
		
			typedef typename Ring::QuotientRing QuotientRing;
			QuotientRing QR(_R.quotient(f));
			PolyDixonDomain<Ring, QuotientRing> DD(_R, QR);
			
			Matrix Mx(M);
					
			Matrix b(_R, Mx.rowdim(), 1);
			for (size_t i = 0; i < b.rowdim(); i++) {
				Polynomial e;
				_RI.random(e, max_deg);
				
				b.setEntry(i, 0, e);
			}
			
			Matrix x(_R, Mx.rowdim(), 1);
			size_t m = 4 * ((max_deg / _R.deg(f)) + 1);
						
			if (!DD.solve(x, Mx, b, f, m)) {
				return false;
			}
			
			Polynomial fm;
			_R.pow(fm, f, m);
			
			_R.assign(minpoly, _R.one);
			for (size_t i = 0; i < x.rowdim(); i++) {
				Polynomial numer, denom, tmp;
				DD.rat_recon(numer, denom, x.getEntry(i, 0), fm);
				
				_R.lcm(tmp, minpoly, denom);
				_R.assign(minpoly, tmp);
			}
			
			return true;
		}
		
		template<class Matrix1>
		bool dixon(Polynomial &mp, const Matrix1 &M) const {
			size_t m = detLimit(M);
			
			bool success = false;
			for (size_t d = 1; d < 11 && !success; d++) {
				Polynomial f;
				_RI.randomIrreducible(f, d++);
				success = dixon(mp, M, f, m);
			}
			_R.monicIn(mp);
			return success;
		}
	}; // end of class PolySmithFormDomain
}

#endif // __LINBOX_poly_smith_form_domain_H
