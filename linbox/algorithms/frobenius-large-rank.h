/* linbox/algorithms/frobenius-large-rank.h
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * EGNS rank-based local Smith form, adapted for FNF computation.
 *
 * Reference:
 *   Elsheikh, Giesbrecht, Novocin, Saunders.
 *   Fast Computation of Smith Forms of Sparse Matrices Over Local Rings.
 *   ISSAC 2012. arXiv:1201.5365v2.
 *
 * REGIME: Suited to "flat" matrices where the max f-exponent e in any
 * invariant factor is small.  Cost per irreducible f of degree d:
 *   O(eta * d * e^2 * n)  field operations.
 * For tall matrices (e = O(n)) the PhiBlackbox grows to O(n^2) x O(n^2)
 * and cost rockets to O(n^5) — use Villard or Gavin instead.
 *
 * FIELD REQUIREMENT: The scalar-Wiedemann rank estimator requires |F| >> n^2
 * to ensure eigenvalue distinctness after diagonal preconditioning.  This is
 * satisfied by word-size primes (e.g. p = 65521 or p = 10000019).  Over
 * GF(2), GF(3), GF(5) the preconditioner fails and blackboxRank will
 * silently underestimate rank.  Defensive guards in localMultiplicities
 * detect the resulting non-monotone rho sequence and bail cleanly rather
 * than crashing with size_t underflow.
 *
 * OUTPUT CONVENTIONS:
 *   solve(fs, ms, A)          — RLE format matching FrobeniusLarge::solve.
 *                               Includes a trailing (1, n_search - total)
 *                               entry to match the FrobeniusLarge convention.
 *   frobeniusInvariants(fs, A, k) — expanded list, top k, NO trailing 1s.
 *                               Matches InvariantFactors::frobeniusInvariants.
 *
 * STILL NEEDED (not yet implemented):
 *   - factorMinpoly for NTL_zz_pEX (extension-field polynomial ring).
 *     Currently only NTL_zz_pX is wired.  The wiring is a one-function swap.
 *   - Block Wiedemann rank to handle small fields robustly.
 *   - Integration stub in test-frobenius-suite.h (add bit7 = EGNS to mask).
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation; either version 2.1 of the License, or
 * (at your option) any later version.
 * ========LICENCE========
 */

#ifndef __LINBOX_frobenius_large_rank_H
#define __LINBOX_frobenius_large_rank_H

#include <vector>
#include <utility>
#include <cassert>
#include <algorithm>
#include <iostream>

#include "linbox/matrix/matrix-domain.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/wiedemann.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/diagonal.h"
#include "linbox/solutions/rank.h"
//#include "linbox/solutions/method.h"

// NTL factoring for factorMinpoly.
// lzz = "long zz" = word-size prime.  Matches NTL_zz_pX.
// For NTL_zz_pEX, swap to <NTL/lzz_pEXFactoring.h> (same interface).
#include <NTL/lzz_pXFactoring.h>

namespace LinBox {

// ============================================================================
// PhiBlackbox
//
// Matrix-free black-box for  phi_l( (xI - A) mod f^l )  over F.
// (EGNS §2.1–2.2, Lemma 7.)
//
// Dimensions: (d*l*n) x (d*l*n) over F.
//
// The underlying matrix is  M = I_n ⊗ C_{f^l}  -  A ⊗ I_{dl}
// where C_{f^l} is the companion matrix of f^l.
//
// Vector layout: v ∈ F^{dl*n} represents n polynomials of degree < dl,
//   v[i * dl + k] = coefficient of x^k in v_i.
//
// apply cost:   dl A-matvecs  +  O(n * dl) field ops.
// applyTranspose: same cost, uses A.applyTranspose.
// ============================================================================

template<class _Field, class _PolynomialRing, class _Blackbox>
class PhiBlackbox {
public:
	typedef _Field          Field;
	typedef _PolynomialRing PolynomialRing;
	typedef _Blackbox       Blackbox;
	typedef typename Field::Element          Element;
	typedef typename PolynomialRing::Element Polynomial;

protected:
	const Blackbox       &_A;
	const PolynomialRing &_R;
	const Field          &_F;
	Polynomial            _fl;          // f^l, monic, degree dl
	std::vector<Element>  _flCoeffs;    // coefficients 0..dl of f^l
	size_t _l;
	size_t _d;                          // deg(f)
	size_t _dl;                         // d * l
	size_t _n;                          // rowdim of A

public:
	PhiBlackbox(const Blackbox &A, const Field &F, const PolynomialRing &R,
	            const Polynomial &f, size_t l)
		: _A(A), _R(R), _F(F), _l(l)
	{
		_d  = (size_t)_R.deg(f);
		_dl = _d * _l;
		_n  = A.rowdim();

		_R.assign(_fl, _R.one);
		for (size_t i = 0; i < _l; i++) _R.mulin(_fl, f);

		_flCoeffs.resize(_dl + 1);
		for (size_t k = 0; k <= _dl; k++) _R.getCoeff(_flCoeffs[k], _fl, k);
		// _flCoeffs[_dl] == 1  (monic)
	}

	size_t rowdim()      const { return _dl * _n; }
	size_t coldim()      const { return _dl * _n; }
	const Field &field() const { return _F; }

	// y = [I_n ⊗ C  -  A ⊗ I] x
	//
	// C = C_{f^l} acts as "multiply by x mod f^l":
	//   (C v)_0   = -a_0 * v_{dl-1}
	//   (C v)_k   = v_{k-1} - a_k * v_{dl-1}   for k = 1..dl-1
	// where a_k = _flCoeffs[k].
	//
	// y_{i,k} = (C v_i)_k  -  (A * slice_k)_i
	template<class OutVec, class InVec>
	OutVec &apply(OutVec &y, const InVec &x) const {
		// Pass 1: deinterleave x into dl coefficient-slices of length n.
		std::vector<BlasVector<Field>> slice;
		slice.reserve(_dl);
		for (size_t k = 0; k < _dl; k++) {
			BlasVector<Field> sk(_F, _n);
			for (size_t i = 0; i < _n; i++) _F.assign(sk[i], x[i * _dl + k]);
			slice.push_back(std::move(sk));
		}

		// Pass 2: A * slice_k for each k  (dl blackbox matvecs — dominates).
		std::vector<BlasVector<Field>> Aslice;
		Aslice.reserve(_dl);
		for (size_t k = 0; k < _dl; k++) {
			BlasVector<Field> ak(_F, _n);
			_A.apply(ak, slice[k]);
			Aslice.push_back(std::move(ak));
		}

		// Pass 3: combine.
		Element tmp, w;
		for (size_t i = 0; i < _n; i++) {
			const Element &lead = slice[_dl - 1][i];

			// k = 0
			_F.mul(tmp, lead, _flCoeffs[0]);
			_F.neg(w, tmp);
			_F.subin(w, Aslice[0][i]);
			_F.assign(y[i * _dl], w);

			// k = 1 .. dl-1
			for (size_t k = 1; k < _dl; k++) {
				_F.assign(w, slice[k - 1][i]);
				_F.mul(tmp, lead, _flCoeffs[k]);
				_F.subin(w, tmp);
				_F.subin(w, Aslice[k][i]);
				_F.assign(y[i * _dl + k], w);
			}
		}
		return y;
	}

	// y = M^T x = [I_n ⊗ C^T  -  A^T ⊗ I] x
	//
	// C^T acts as "reverse shift":
	//   (C^T v)_k     = v_{k+1}              for k = 0..dl-2
	//   (C^T v)_{dl-1} = -Σ_j a_j v_j
	//
	// Requires A to implement applyTranspose.
	template<class OutVec, class InVec>
	OutVec &applyTranspose(OutVec &y, const InVec &x) const {
		// Pass 1: deinterleave.
		std::vector<BlasVector<Field>> slice;
		slice.reserve(_dl);
		for (size_t k = 0; k < _dl; k++) {
			BlasVector<Field> sk(_F, _n);
			for (size_t i = 0; i < _n; i++) _F.assign(sk[i], x[i * _dl + k]);
			slice.push_back(std::move(sk));
		}

		// Pass 2: A^T * slice_k.
		std::vector<BlasVector<Field>> ATslice;
		ATslice.reserve(_dl);
		for (size_t k = 0; k < _dl; k++) {
			BlasVector<Field> ak(_F, _n);
			_A.applyTranspose(ak, slice[k]);
			ATslice.push_back(std::move(ak));
		}

		// Pass 3: combine.
		Element tmp, w;
		for (size_t i = 0; i < _n; i++) {
			// k = 0..dl-2:  (C^T v)_k = v_{k+1}
			for (size_t k = 0; k + 1 < _dl; k++) {
				_F.assign(w, slice[k + 1][i]);
				_F.subin(w, ATslice[k][i]);
				_F.assign(y[i * _dl + k], w);
			}
			// k = dl-1:  (C^T v)_{dl-1} = -Σ_j a_j v_j
			_F.assign(w, _F.zero);
			for (size_t j = 0; j < _dl; j++) {
				_F.mul(tmp, _flCoeffs[j], slice[j][i]);
				_F.subin(w, tmp);
			}
			_F.subin(w, ATslice[_dl - 1][i]);
			_F.assign(y[i * _dl + (_dl - 1)], w);
		}
		return y;
	}
};

// ============================================================================
// FrobeniusLargeRank
//
// EGNS rank-based FNF.
// Provides the same public API as FrobeniusLarge (frobenius-large.h) and
// InvariantFactors (invariant-factors.h).
// ============================================================================

template<class _PolynomialRing>
class FrobeniusLargeRank {
public:
	typedef _PolynomialRing                    PolynomialRing;
	typedef typename PolynomialRing::Element    Polynomial;
	typedef typename PolynomialRing::Coeff      Coeff;
	typedef typename PolynomialRing::CoeffField Field;
	typedef typename Field::Element             Element;
	typedef typename Field::RandIter            RandIter;

protected:
	Field          _F;
	RandIter       _RI;
	PolynomialRing _R;

public:
	FrobeniusLargeRank(const PolynomialRing &R)
		: _F(R.getCoeffField()), _RI(_F), _R(R) {}

	// ================================================================
	// minpoly — scalar Wiedemann / Berlekamp-Massey.
	// Matches FrobeniusLarge::minpoly exactly.
	// ================================================================
	template<class Blackbox>
	void minpoly(Polynomial &g, const Blackbox &A) {
		typedef BlackboxContainer<Field, Blackbox> Sequence;
		Sequence seq(&A, _F, _RI);
		MasseyDomain<Field, Sequence> MD(&seq, 20);
		BlasVector<Field> phi(_F);
		size_t deg;
		MD.minpoly(phi, deg);
		_R.init(g, phi);
	}

	// ================================================================
	// blackboxRank — delegates to LinBox::rank (Wiedemann method).
	//
	// WHY NOT DIAGONAL PRECONDITIONING:
	// PhiBlackbox matrices in the EGNS pipeline are nilpotent:
	//   - Eigenvalue 0: phi_l = I_n⊗J_l(0) - A⊗I_l is nilpotent when A is.
	//   - Eigenvalue λ: phi_1 = -(A-λI) is nilpotent; higher levels also.
	// A diagonal D does NOT break nilpotency — (D·M)^k = 0 when M^k = 0.
	// So BM returns x^k (nilpotency index, not rank), and deg-1 can be
	// far from the true rank, producing wrong rho values.
	//
	// LinBox::rank uses Toeplitz/butterfly preconditioning which DOES
	// break nilpotency (T·M is not nilpotent for generic full-rank T),
	// giving correct rank for nilpotent and near-nilpotent matrices.
	// ================================================================
	template<class Blackbox>
	size_t blackboxRank(const Blackbox &M) {
		size_t r = 0;
		LinBox::rank(r, M, Method::Wiedemann());
		return r;
	}

	// ================================================================
	// localMultiplicities — EGNS Algorithm 1 for one irreducible f.
	//
	// After this call, rs has length e and rs[i] is the number of
	// local invariant factors of (xI - A) at f equal to f^i.
	//
	// Setting e = (exponent of f in minpoly) + 1 captures all multiplicities.
	//
	// Defensive guards: if the rho sequence is non-monotone (rank failure,
	// usually from a small field), rs is zeroed and the function returns
	// false.  The caller should treat a false return as "result unreliable."
	// ================================================================
	template<class Blackbox>
	bool localMultiplicities(
		std::vector<size_t> &rs,
		const Blackbox      &A,
		const Polynomial    &f,
		size_t               e)
	{
		rs.assign(e, 0);
		if (e == 0) return true;

		size_t d = (size_t)_R.deg(f);

		// e rank computations: rho[l-1] = rank( phi_l((xI-A) mod f^l) )
		std::vector<size_t> rho(e);
		for (size_t l = 1; l <= e; l++) {
			PhiBlackbox<Field, PolynomialRing, Blackbox> phi(A, _F, _R, f, l);
			rho[l - 1] = blackboxRank(phi);
		}

		// Defensive check: rho must be non-decreasing (Theorem 4 guarantees
		// strict increase by multiples of d when the rank algorithm is exact).
		for (size_t l = 1; l < e; l++) {
			if (rho[l] < rho[l - 1] || (rho[l] - rho[l - 1]) % d != 0) {
				// Rank underestimation detected — field likely too small.
				// Zero rs and signal failure rather than letting size_t
				// arithmetic underflow to astronomically large values.
				rs.assign(e, 0);
				return false;
			}
		}

		// Triangular solve — EGNS Corollary 5.
		//
		// rho[l] - rho[l-1] = d * S_l   where S_l = sum_{i=0}^{l} r_i
		// (prefix sum of the r_i sequence)
		//
		// So sigma[l] := (rho[l] - rho[l-1]) / d = S_l
		// and r_l = S_l - S_{l-1} = sigma[l] - sigma[l-1]   (forward diff)
		std::vector<size_t> sigma(e);
		sigma[0] = rho[0] / d;
		for (size_t l = 1; l < e; l++)
			sigma[l] = (rho[l] - rho[l - 1]) / d;

		rs[0] = sigma[0];
		for (size_t l = 1; l < e; l++)
			rs[l] = sigma[l] - sigma[l - 1];

		return true;
	}

	// ================================================================
	// solve — full FNF, RLE output matching FrobeniusLarge::solve.
	//
	// fs[i]  = i-th distinct invariant factor (largest first by degree)
	// ms[i]  = run-length count of fs[i]
	//
	// Includes a trailing (1, n_search − total) entry to match the
	// FrobeniusLarge convention where n_search = n − deg(minpoly) + 2.
	// ================================================================
	template<class Blackbox>
	void solve(
		std::vector<Polynomial> &fs,
		std::vector<size_t>     &ms,
		const Blackbox          &A)
	{
		assert(A.rowdim() == A.coldim());
		fs.clear();
		ms.clear();

		// (a) Minimal polynomial.
		Polynomial mp;
		minpoly(mp, A);

		// (b) Factor minpoly = prod f_j^{e_j}.
		std::vector<std::pair<Polynomial, size_t>> factors;
		factorMinpoly(factors, mp);

		if (factors.empty()) return;  // degenerate input (degree-0 minpoly)

		// (c) Per-factor local multiplicities → partitions.
		std::vector<Polynomial>          distinctFs;
		std::vector<std::vector<size_t>> partitions;

		for (size_t fi = 0; fi < factors.size(); fi++) {
			const Polynomial &f = factors[fi].first;
			size_t eAlg         = factors[fi].second + 1;

			std::vector<size_t> rs;
			bool ok = localMultiplicities(rs, A, f, eAlg);
			if (!ok) continue;  // rank failure — skip this factor

			// Build partition (nonincreasing exponents).
			// rs[0] counts unit invariant factors — omit from FNF output.
			std::vector<size_t> p;
			for (size_t i = eAlg; i-- > 1; )
				for (size_t k = 0; k < rs[i]; k++) p.push_back(i);

			if (!p.empty()) {
				distinctFs.push_back(f);
				partitions.push_back(std::move(p));
			}
		}

		// (d) Assemble: invariants[j] = prod_f f^{partitions[fi][j]}.
		size_t kTotal = 0;
		for (auto &p : partitions) kTotal = std::max(kTotal, p.size());

		if (kTotal == 0) return;

		std::vector<Polynomial> invariants(kTotal);
		for (size_t j = 0; j < kTotal; j++) _R.assign(invariants[j], _R.one);

		for (size_t fi = 0; fi < distinctFs.size(); fi++) {
			const Polynomial &f = distinctFs[fi];
			const auto       &p = partitions[fi];

			// Cache f^0 .. f^maxExp.
			size_t maxExp = p.empty() ? 0 : p.front();
			std::vector<Polynomial> fpow(maxExp + 1);
			_R.assign(fpow[0], _R.one);
			for (size_t k = 1; k <= maxExp; k++) {
				_R.assign(fpow[k], fpow[k - 1]);
				_R.mulin(fpow[k], f);
			}

			for (size_t j = 0; j < p.size(); j++)
				_R.mulin(invariants[j], fpow[p[j]]);
		}

		// (e) RLE in nonincreasing-degree order.
		size_t i = 0;
		while (i < kTotal) {
			size_t j = i;
			while (j < kTotal && _R.areEqual(invariants[j], invariants[i])) j++;
			fs.push_back(invariants[i]);
			ms.push_back(j - i);
			i = j;
		}

		// (f) Trailing-1 entry matching FrobeniusLarge::solve convention.
		//
		// Villard's threshold search anchors at n_search = n - deg(f1) + 2,
		// emitting a synthetic trailing (1, n_search - total) entry.
		// We replicate this so the output format is byte-compatible for
		// direct comparison in the test suite.
		// Note: this entry does NOT appear in frobeniusInvariants().
		size_t deg_mp = (size_t)_R.deg(mp);
		size_t n      = A.rowdim();
		if (deg_mp < n) {
			size_t n_search = n - deg_mp + 2;
			size_t total    = 0;
			for (size_t c : ms) total += c;
			if (total < n_search) {
				fs.push_back(_R.one);
				ms.push_back(n_search - total);
			}
		}
	}

	// ================================================================
	// solve — expanded list with repeats (FrobeniusLarge overload).
	// ================================================================
	template<class Blackbox>
	void solve(std::vector<Polynomial> &fs, const Blackbox &A) {
		std::vector<Polynomial> fsu;
		std::vector<size_t>     msu;
		solve(fsu, msu, A);
		fs.clear();
		for (size_t i = 0; i < fsu.size(); i++)
			for (size_t j = 0; j < msu[i]; j++)
				fs.push_back(fsu[i]);
	}

	// ================================================================
	// frobeniusInvariants — matches InvariantFactors::frobeniusInvariants.
	//
	// Returns the full FNF in nonincreasing degree order, truncated to
	// the top k factors when k > 0.  NO trailing 1s (unlike solve above).
	//
	// IMPORTANT: EGNS always computes the full FNF regardless of k.
	// It pays full cost whether k=1 or k=n — it has no "targeted kth"
	// query.  Use Villard (kthInvariantFactor) for targeted queries.
	// ================================================================
	template<class Blackbox>
	std::vector<Polynomial> &frobeniusInvariants(
		std::vector<Polynomial> &fs,
		const Blackbox          &A,
		size_t                   k = 0)
	{
		// Compute full FNF via the RLE solve, then expand — but strip
		// the trailing-1 entry before expanding.
		std::vector<Polynomial> fsu;
		std::vector<size_t>     msu;
		solve(fsu, msu, A);

		fs.clear();

		// Drop trailing 1 entry (the Villard-compatibility artefact).
		size_t limit = fsu.size();
		if (limit > 0 && _R.isOne(fsu.back())) limit--;

		for (size_t i = 0; i < limit; i++)
			for (size_t j = 0; j < msu[i]; j++) {
				if (k > 0 && fs.size() >= k) goto done;
				fs.push_back(fsu[i]);
			}
		done:
		return fs;
	}

protected:
	// ================================================================
	// factorMinpoly — factor monic polynomial into irreducibles.
	//
	// Delegates to NTL::CanZass.
	// Assumes Polynomial == NTL::zz_pX  (i.e. _PolynomialRing = NTL_zz_pX).
	//
	// TODO: For NTL_zz_pEX, replace with the zz_pEX variant of CanZass
	//   from <NTL/lzz_pEXFactoring.h> — same interface, swap type names.
	// ================================================================
	void factorMinpoly(
		std::vector<std::pair<Polynomial, size_t>> &factors,
		const Polynomial                           &p)
	{
		factors.clear();
		if (_R.deg(p) <= 0) return;   // constant polynomial — no factors

		NTL::vec_pair_zz_pX_long ntl_factors;
		NTL::CanZass(ntl_factors, p);

		for (long i = 0; i < ntl_factors.length(); i++) {
			NTL::MakeMonic(ntl_factors[i].a);
			factors.emplace_back(
				ntl_factors[i].a,
				static_cast<size_t>(ntl_factors[i].b));
		}
	}
};

} // namespace LinBox

#endif // __LINBOX_frobenius_large_rank_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s