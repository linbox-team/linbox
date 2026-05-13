#ifndef __LINBOX_frobenius_large_search_H
#define __LINBOX_frobenius_large_search_H

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
#include "linbox/blackbox/toeplitz.h"

namespace LinBox
{

template<class _PolynomialRing>
class FrobeniusLargeSearch {
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
	FrobeniusLargeSearch(const PolynomialRing &R)
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

		// exponential scan from l to find bracket containing first transition
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

		// binary search [lo, hi] for exact transition point
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

#endif //__LINBOX_frobenius_large_search_H