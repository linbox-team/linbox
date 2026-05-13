/* linbox/algorithms/frobenius-large-generic.h
 * Copyright (C) 2018 Gavin Harrison
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Single-file implementation of all Frobenius large-field variants.
 * Combines Toeplitz, Butterfly, and Dense preconditioners with
 * binary and exponential threshold search strategies.
 *
 * The 6 original class names are preserved as template aliases:
 *   FrobeniusLarge              (Toeplitz   + binary search)
 *   FrobeniusLargeSearch        (Toeplitz   + exponential search)
 *   FrobeniusLargeButterfly     (Butterfly  + binary search)
 *   FrobeniusLargeButterflySearch (Butterfly + exponential search)
 *   FrobeniusLargeDense         (Dense      + binary search)
 *   FrobeniusLargeDenseSearch   (Dense      + exponential search)
 *
 * Adding a new preconditioner:
 *   1. Write a new derived class inheriting FrobeniusLargeBase,
 *      implementing only kthInvariantFactorImpl.
 *   2. Add two using aliases (binary + exponential).
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 * ========LICENCE========
 */

#ifndef __LINBOX_frobenius_large_generic_H
#define __LINBOX_frobenius_large_generic_H

#include <list>
#include <vector>
#include <math.h>
#include <algorithm>
#include <iostream>
#include <cassert>

#include "linbox/matrix/matrix-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/wiedemann.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sum.h"
#include "linbox/blackbox/toeplitz.h"
#include "linbox/blackbox/butterfly.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/matrix/sparse-matrix.h"

namespace LinBox
{

// ============================================================
// Search strategy selector
// ============================================================
enum SearchStrategy { BinarySearch, ExponentialSearch };

// ============================================================
// CRTP base class
// Holds all shared infrastructure and search logic.
// Derived must implement kthInvariantFactorImpl.
// ============================================================
template<class _PolynomialRing, class Derived, SearchStrategy _Search>
class FrobeniusLargeBase {
public:
	typedef _PolynomialRing PolynomialRing;
	typedef typename PolynomialRing::Element Polynomial;
	typedef typename PolynomialRing::Coeff Coeff;

	typedef typename PolynomialRing::CoeffField Field;
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;

	typedef MatrixDomain<Field> MatrixDom;
	typedef typename MatrixDom::OwnMatrix Matrix;

protected:
	Field _F;
	RandIter _RI;
	PolynomialRing _R;
	MatrixDom _MD;

public:
	FrobeniusLargeBase(const PolynomialRing &R)
		: _F(R.getCoeffField()), _RI(_F), _R(R), _MD(_F) {}

	// ----------------------------------------------------------
	// Helpers
	// ----------------------------------------------------------
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

	// Dispatches to the derived class preconditioner implementation.
	template<class Blackbox>
	void kthInvariantFactor(
		Polynomial &fk,
		const Blackbox &A,
		const Polynomial &m,
		size_t k)
	{
		static_cast<Derived*>(this)->kthInvariantFactorImpl(fk, A, m, k);
	}

	// ----------------------------------------------------------
	// Binary threshold search
	// ----------------------------------------------------------
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
			ms.push_back(m - l + 1);
			return;
		}
		if (l == m - 1) {
			if (_R.areEqual(fl, fm)) {
				fs.push_back(fl);
				ms.push_back(2);
			} else {
				fs.push_back(fl); ms.push_back(1);
				fs.push_back(fm); ms.push_back(1);
			}
			return;
		}
		size_t k = (size_t)std::ceil((l + m) / 2.0);
		Polynomial fk;
		kthInvariantFactor(fk, A, fl, k);

		std::vector<Polynomial> gs; std::vector<size_t> as;
		thresholdSearch(gs, as, A, l, fl, k, fk);

		std::vector<Polynomial> hs; std::vector<size_t> bs;
		thresholdSearch(hs, bs, A, k, fk, m, fm);

		for (size_t i = 0; i < as.size() - 1; i++) { fs.push_back(gs[i]); ms.push_back(as[i]); }
		fs.push_back(gs[gs.size() - 1]);
		ms.push_back(as[as.size() - 1] + bs[0] - 1);
		for (size_t i = 1; i < bs.size(); i++) { fs.push_back(hs[i]); ms.push_back(bs[i]); }
	}

	// ----------------------------------------------------------
	// Exponential threshold search
	// ----------------------------------------------------------
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
			lo = (step == 1) ? l : l + step / 2;
			hi = l + step;
			_R.assign(f_hi, f_probe);
		} else {
			lo = l + step / 2;
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

	// ----------------------------------------------------------
	// solve (run-length encoded output)
	// ----------------------------------------------------------

	/** fs = distinct invariant factors in nonincreasing order by degree.
	 *  ms[i] = run-length count of fs[i] in the full list.
	 *  If limit > 0, only the first limit invariant factors are found.
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
			if (_Search == ExponentialSearch)
				exponentialThresholdSearch(fs, ms, A, 1, f1, limit, flimit);
			else
				thresholdSearch(fs, ms, A, 1, f1, limit, flimit);
			return;
		}

		if (_Search == ExponentialSearch)
			exponentialThresholdSearch(fs, ms, A, 1, f1, n, _R.one);
		else
			thresholdSearch(fs, ms, A, 1, f1, n, _R.one);
	}

	/** fs = full invariant factor list (with repeats) in nonincreasing order.
	 *  If limit > 0, only the first limit invariants are found.
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

// ============================================================
// Toeplitz preconditioner
// ============================================================
template<class _PolynomialRing, SearchStrategy _Search>
class FrobeniusLargeToeplitzImpl
	: public FrobeniusLargeBase<_PolynomialRing,
	                             FrobeniusLargeToeplitzImpl<_PolynomialRing, _Search>,
	                             _Search>
{
	typedef FrobeniusLargeBase<_PolynomialRing,
	                            FrobeniusLargeToeplitzImpl<_PolynomialRing, _Search>,
	                            _Search> Base;
public:
	typedef typename Base::Polynomial Polynomial;
	typedef typename Base::Field Field;
	typedef Toeplitz<Field, _PolynomialRing> Toep;

	FrobeniusLargeToeplitzImpl(const _PolynomialRing &R) : Base(R) {}

	template<class Blackbox>
	void kthInvariantFactorImpl(
		Polynomial &fk,
		const Blackbox &A,
		const Polynomial &m,
		size_t k)
	{
		size_t n = A.rowdim();
		Polynomial u, v;
		this->randomPolynomial(u, n + k - 3);
		this->randomPolynomial(v, n + k - 3);
		Toep U(this->_R, u, n, k - 1);
		Toep V(this->_R, v, k - 1, n);
		Compose<Toep, Toep> B(U, V);
		Sum<Blackbox, Compose<Toep, Toep>> Ak(A, B);
		this->minpoly(fk, Ak);
		this->_R.gcdin(fk, m);
	}
};

// ============================================================
// Butterfly preconditioner
// ============================================================
template<class _PolynomialRing, SearchStrategy _Search>
class FrobeniusLargeButterflyImpl
	: public FrobeniusLargeBase<_PolynomialRing,
	                             FrobeniusLargeButterflyImpl<_PolynomialRing, _Search>,
	                             _Search>
{
	typedef FrobeniusLargeBase<_PolynomialRing,
	                            FrobeniusLargeButterflyImpl<_PolynomialRing, _Search>,
	                            _Search> Base;
public:
	typedef typename Base::Polynomial Polynomial;
	typedef typename Base::Element Element;
	typedef typename Base::Field Field;
	typedef CekstvSwitch<Field> Switch;
	typedef Butterfly<Field, Switch> BB;

	FrobeniusLargeButterflyImpl(const _PolynomialRing &R) : Base(R) {}

	template<class Blackbox>
	void kthInvariantFactorImpl(
		Polynomial &fk,
		const Blackbox &A,
		const Polynomial &m,
		size_t k)
	{
		size_t n = A.rowdim();

		typename Switch::Factory facU(this->_RI);
		typename Switch::Factory facV(this->_RI);
		BB U(this->_F, n, facU);
		BB V(this->_F, n, facV);

		typedef TransposeOwner<BB> BBt;
		BBt Vt(V);

		typedef SparseMatrix<Field> Diag;
		Diag D(this->_F, n, n);
		for (size_t i = 0; i < k - 1; ++i) {
			Element e;
			do { this->_RI.random(e); } while (this->_F.isZero(e));
			D.setEntry(n - i - 1, n - i - 1, e);
		}

		typedef Compose<Diag, BB> DU_t;
		DU_t DU(D, U);
		typedef Compose<BBt, DU_t> VtDU_t;
		VtDU_t B(Vt, DU);
		typedef Sum<Blackbox, VtDU_t> Ak_t;
		Ak_t Ak(A, B);

		this->minpoly(fk, Ak);
		this->_R.gcdin(fk, m);
	}
};

// ============================================================
// Dense preconditioner
// ============================================================
template<class _PolynomialRing, SearchStrategy _Search>
class FrobeniusLargeDenseImpl
	: public FrobeniusLargeBase<_PolynomialRing,
	                             FrobeniusLargeDenseImpl<_PolynomialRing, _Search>,
	                             _Search>
{
	typedef FrobeniusLargeBase<_PolynomialRing,
	                            FrobeniusLargeDenseImpl<_PolynomialRing, _Search>,
	                            _Search> Base;
public:
	typedef typename Base::Polynomial Polynomial;
	typedef typename Base::Element Element;
	typedef typename Base::Matrix Matrix;

	FrobeniusLargeDenseImpl(const _PolynomialRing &R) : Base(R) {}

	template<class Blackbox>
	void kthInvariantFactorImpl(
		Polynomial &fk,
		const Blackbox &A,
		const Polynomial &m,
		size_t k)
	{
		size_t n = A.rowdim();
		Matrix U(this->_F, n, k - 1);
		Matrix V(this->_F, k - 1, n);

		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < k - 1; ++j) {
				Element e; this->_RI.random(e);
				U.setEntry(i, j, e);
			}
		for (size_t i = 0; i < k - 1; ++i)
			for (size_t j = 0; j < n; ++j) {
				Element e; this->_RI.random(e);
				V.setEntry(i, j, e);
			}

		Compose<Matrix, Matrix> B(U, V);
		Sum<Compose<Matrix, Matrix>, Blackbox> Ak(B, A);

		this->minpoly(fk, Ak);
		this->_R.gcdin(fk, m);
	}
};

// ============================================================
// Public aliases — same names as the original 6 files
// ============================================================
template<class R> using FrobeniusLarge              = FrobeniusLargeToeplitzImpl<R,  BinarySearch>;
template<class R> using FrobeniusLargeSearch        = FrobeniusLargeToeplitzImpl<R,  ExponentialSearch>;
template<class R> using FrobeniusLargeButterfly     = FrobeniusLargeButterflyImpl<R, BinarySearch>;
template<class R> using FrobeniusLargeButterflySearch = FrobeniusLargeButterflyImpl<R, ExponentialSearch>;
template<class R> using FrobeniusLargeDense         = FrobeniusLargeDenseImpl<R,     BinarySearch>;
template<class R> using FrobeniusLargeDenseSearch   = FrobeniusLargeDenseImpl<R,     ExponentialSearch>;

} // namespace LinBox

#endif // __LINBOX_frobenius_large_generic_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s