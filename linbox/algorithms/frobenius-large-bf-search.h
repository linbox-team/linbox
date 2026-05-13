/* linbox/algorithms/frobenius-large-bf-search.h
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
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

#ifndef __LINBOX_frobenius_large_butterfly_search_H
#define __LINBOX_frobenius_large_butterfly_search_H

#include <list>
#include <vector>
#include <math.h>

#include <algorithm>
#include <iostream>

#include "linbox/matrix/matrix-domain.h"

#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/wiedemann.h"

#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sum.h"
#include "linbox/blackbox/butterfly.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/matrix/sparse-matrix.h"

namespace LinBox
{

// Butterfly preconditioner + exponential threshold search
// PolynomialRing = NTL_zz_pX or NTL_zz_pEX
template<class _PolynomialRing>
class FrobeniusLargeButterflySearch {
public:
	typedef _PolynomialRing PolynomialRing;
	typedef typename PolynomialRing::Element Polynomial;
	typedef typename PolynomialRing::Coeff Coeff;

	typedef typename PolynomialRing::CoeffField Field;
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	typedef MatrixDomain<Field> MatrixDom;
	typedef typename MatrixDom::OwnMatrix Matrix;

	typedef CekstvSwitch<Field> Switch;
	typedef Butterfly<Field, Switch> BB;

protected:
	Field _F;
	RandIter _RI;
	PolynomialRing _R;
	MatrixDom _MD;

public:
	FrobeniusLargeButterflySearch(const PolynomialRing &R)
		: _F(R.getCoeffField()), _RI(_F), _R(R), _MD(_F) {}

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

		typename Switch::Factory facU(_RI);
		typename Switch::Factory facV(_RI);

		BB U(_F, n, facU);
		BB V(_F, n, facV);

		typedef TransposeOwner<BB> BBt;
		BBt Vt(V);

		typedef SparseMatrix<Field> Diag;
		Diag D(_F, n, n);
		for (size_t i = 0; i < k-1; ++i) {
			Element e;
			do { _RI.random(e); } while (_F.isZero(e));
			D.setEntry(n-i-1, n-i-1, e);
		}

		typedef Compose<Diag, BB> DU_t;
		DU_t DU(D, U);

		typedef Compose<BBt, DU_t> VtDU_t;
		VtDU_t B(Vt, DU);
		typedef Sum<Blackbox, VtDU_t> Ak_t;
		Ak_t Ak(A, B);

		minpoly(fk, Ak);
		_R.gcdin(fk, m);
	}

	template<class Blackbox>
	void exponentialThresholdSearch(
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
			ms.push_back(m - l + 1);
			return;
		}

		if (l == m - 1) {
			fs.push_back(fl); ms.push_back(1);
			fs.push_back(fm); ms.push_back(1);
			return;
		}

		size_t step = 1;
		Polynomial f_probe;
		bool found = false;

		while (l + step < m) {
			kthInvariantFactor(f_probe, A, fl, l + step);
			if (!_R.areEqual(f_probe, fl)) { found = true; break; }
			step *= 2;
		}

		size_t lo, hi;
		Polynomial f_hi;

		if (found) {
			lo = (step == 1) ? l : l + step/2;
			hi = l + step;
			_R.assign(f_hi, f_probe);
		} else {
			lo = l + step/2;
			hi = m;
			_R.assign(f_hi, fm);
		}

		while (hi > lo + 1) {
			size_t mid = (lo + hi + 1) / 2;
			Polynomial f_mid;
			kthInvariantFactor(f_mid, A, fl, mid);
			if (_R.areEqual(f_mid, fl)) lo = mid;
			else { hi = mid; _R.assign(f_hi, f_mid); }
		}

		fs.push_back(fl);
		ms.push_back(lo - l + 1);

		exponentialThresholdSearch(fs, ms, A, hi, f_hi, m, fm);
	}

	/** fs is the distinct invariant factors of A in nonincreasing order by degree.
	 *  ms[i] is the run-length count of fs[i] in the full invariant factor list.
	 *  If limit is positive, only the first limit invariant factors are found.
	 */
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

		Polynomial f1;
		minpoly(f1, A);

		if (_R.deg(f1) == A.rowdim()) {
			fs.push_back(f1);
			ms.push_back(1);
			return;
		}

		size_t n = A.rowdim() - _R.deg(f1) + 2;
		if (0 < limit && limit < n) {
			Polynomial flimit;
			kthInvariantFactor(flimit, A, f1, limit);
			exponentialThresholdSearch(fs, ms, A, 1, f1, limit, flimit);
			return;
		}

		exponentialThresholdSearch(fs, ms, A, 1, f1, n, _R.one);
	}

	/** fs is the invariant factor list of A in nonincreasing order by degree.
	 *  If limit is positive, only the first limit invariants are found.
	 */
	template<class Blackbox>
	void frobeniusInvariants(
		std::vector<Polynomial> &fs,
		const Blackbox &A,
		size_t limit = 0)
	{ solve(fs, A, limit); }

	template<class Blackbox>
	void solve(
		std::vector<Polynomial> &fs,
		const Blackbox &A,
		size_t limit = 0)
	{
		std::vector<Polynomial> fsu;
		std::vector<size_t> ms;
		solve(fsu, ms, A, limit);
		for (size_t i = 0; i < fsu.size(); ++i)
			for (size_t j = 0; j < ms[i]; ++j)
				fs.push_back(fsu[i]);
	}
};

}

#endif //__LINBOX_frobenius_large_butterfly_search_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s