/* linbox/algorithms/frobenius-large.h
 * Copyright (C) 2018 Gavin Harrison
 *
 * Written by Gavin Harrison <gavin.har@gmail.com>
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
 */

#ifndef __LINBOX_frobenius_large_H
#define __LINBOX_frobenius_large_H

#include <list>
#include <vector>
#include <math.h> 

#include <algorithm>
#include <iostream>

#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/wiedemann.h"

#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sum.h"
#include "linbox/blackbox/toeplitz.h"

namespace LinBox
{

// source
// Computing the Frobenius Normal Form of a Sparse Matrix (2000)
// URL: https://doi.org/10.1007/978-3-642-57201-2_30

// PolynomialRing = NTL_zz_pX or NTL_zz_pEX
template<class _PolynomialRing>
class FrobeniusLarge {
public:
	typedef _PolynomialRing PolynomialRing;
	typedef typename PolynomialRing::Element Polynomial;
	typedef typename PolynomialRing::Coeff Coeff;
	
	typedef typename PolynomialRing::CoeffField Field;
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	
	typedef MatrixDomain<Field> MatrixDom;
	typedef typename MatrixDom::OwnMatrix Matrix;
	
	typedef Toeplitz<Field, PolynomialRing> Toep;
		
protected:
	Field _F;
	RandIter _RI;
	PolynomialRing _R;
	MatrixDom _MD;
public:
	FrobeniusLarge(const PolynomialRing &R) 
		: _F(R.getCoeffField()), _RI(_F), _R(R), _MD(_F) {}
	
	void randomPolynomial(Polynomial &f, size_t d) const {
		_R.assign(f, _R.zero);
		
		for (size_t i = 0; i <= d; i++) {
			Coeff c;
			_RI.random(c);
			_R.setCoeff(f, i, c);
		}
	}
	
	template<class Blackbox>
	void minpoly(Polynomial &g, const Blackbox &A) {
		typedef BlackboxContainer<Field, Blackbox> Sequence;
		
		Sequence seq(&A, _F, _RI);
		MasseyDomain<Field, Sequence> MasseyDom(&seq, 20);
		
		BlasVector<Field> phi(_F);
		size_t deg;
		MasseyDom.minpoly(phi, deg);
		
		_R.init(g, phi);
	}
	
	template<class Blackbox>
	void kthInvariantFactor(
		Polynomial &fk, 
		const Blackbox &A, 
		const Polynomial &m, 
		size_t k) 
	{
		size_t n = A.rowdim();
		
		Polynomial u, v;
		randomPolynomial(u, n+k-3);
		randomPolynomial(v, n+k-3);
		
		Toep U(_R, u, n, k-1);
		Toep V(_R, v, k-1, n);
		
		Compose<Toep, Toep> B(U, V);
		Sum<Blackbox, Compose<Toep, Toep>> Ak(A, B);
		
		minpoly(fk, Ak);
		_R.gcdin(fk, m);
	}
	
	template<class Blackbox>
	void thresholdSearch(
		std::vector<Polynomial> &fs, 
		std::vector<size_t> &ms,
		const Blackbox &A,
		size_t l,
		const Polynomial &fl,
		size_t m,
		const Polynomial &fm) 
	{
		if (_R.areEqual(fl, fm)) {
			fs.push_back(fl);
			ms.push_back(m-l+1);
			return;
		}
		
		if (l == m-1) {
			if (_R.areEqual(fl, fm)) {
				fs.push_back(fl);
				ms.push_back(2);
			} else {
				fs.push_back(fl);
				fs.push_back(fm);
				ms.push_back(1);
				ms.push_back(1);
			}
			return;
		}
		
		size_t k = (size_t) std::ceil((l + m)/2.0);
		
		Polynomial fk;
		kthInvariantFactor(fk, A, fl, k);
		
		std::vector<Polynomial> gs;
		std::vector<size_t> as;
		thresholdSearch(gs, as, A, l, fl, k, fk);
				
		std::vector<Polynomial> hs;
		std::vector<size_t> bs;
		thresholdSearch(hs, bs, A, k, fk, m, fm);
				
		for (size_t i = 0; i < as.size() - 1; i++) {
			fs.push_back(gs[i]);
			ms.push_back(as[i]);
		}
		fs.push_back(gs[gs.size() - 1]);
		ms.push_back(as[as.size() - 1] + bs[0] - 1);
		for (size_t i = 1; i < bs.size(); i++) {
			fs.push_back(hs[i]);
			ms.push_back(bs[i]);
		}
	}
	
	template<class Blackbox>
	void solve(
		std::vector<Polynomial> &fs,
		std::vector<size_t> &ms,
		const Blackbox &A,
		size_t limit = 0)
	{
		assert(A.rowdim() == A.coldim());
		fs.clear();
		ms.clear();
		
		Polynomial f1, fn;
		minpoly(f1, A);
		
		if (_R.deg(f1) == A.rowdim()) {
			fs.push_back(f1);
			ms.push_back(1);
			return;
		}
		
		size_t n = A.rowdim() - _R.deg(f1) + 2;
		if (0 < limit && limit < n) {
			kthInvariantFactor(fn, A, f1, n);
			thresholdSearch(fs, ms, A, 1, f1, n, fn);
			n = std::min(limit, n);
			return;
		}
		
		thresholdSearch(fs, ms, A, 1, f1, n, _R.one);
	}
};

}

#endif //__LINBOX_frobenius_large_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
