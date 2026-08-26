/* linbox/algorithms/consecutive-invariant-factors.h
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 *
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * ========LICENCE========
 * This file is part of the library LinBox.
 *
 * LinBox is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * LinBox is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 * ========LICENCE========
 */

#ifndef __LINBOX_consecutive_invariant_factors_H
#define __LINBOX_consecutive_invariant_factors_H

#include <cassert>
#include <cstddef>
#include <vector>

#include "linbox/algorithms/invariant-factors.h"
#include "linbox/blackbox/butterfly.h"
#include "linbox/blackbox/compose.h"
#include "linbox/blackbox/sum.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/field/hom.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/sparse-matrix.h"

namespace LinBox
{

/** Selects the rank-shift representation used before the LIF computation. */
enum ConsecutiveLIFShift {
	ConsecutiveLIFAuto,
	ConsecutiveLIFDense,
	ConsecutiveLIFButterfly
};

namespace Protected
{

template<class Field> class ConsecutiveLIFSwitchFactory;

/**
 * Algebraic pass/swap switch used by the representative-minor proof.
 *
 * Its 2 by 2 matrix is
 *
 *       [ 1-a    a  ].
 *       [  a    1-a ]
 *
 * Thus a=0 is the identity and a=1 is the transposition.  These two
 * specializations are the routing property needed by the nested-minor
 * argument; CekstvSwitch does not have this property.
 */
template<class _Field>
class ConsecutiveLIFSwitch {
public:
	typedef _Field Field;
	typedef typename Field::Element Element;
	typedef ConsecutiveLIFSwitch<Field> Self_t;
	typedef ConsecutiveLIFSwitchFactory<Field> Factory;

	ConsecutiveLIFSwitch() {}
	explicit ConsecutiveLIFSwitch(const Element &a) : _a(a) {}

	bool apply(const Field &F, Element &x, Element &y) const
	{
		Element oneMinusA, newX, newY;
		F.sub(oneMinusA, F.one, _a);
		F.mul(newX, oneMinusA, x);
		F.axpyin(newX, _a, y);
		F.mul(newY, _a, x);
		F.axpyin(newY, oneMinusA, y);
		F.assign(x, newX);
		F.assign(y, newY);
		return true;
	}

	bool applyTranspose(const Field &F, Element &x, Element &y) const
	{
		return apply(F, x, y);
	}

	template<typename _Tp1>
	struct rebind {
		typedef ConsecutiveLIFSwitch<_Tp1> other;

		void operator()(other &Ap, const Self_t &A,
		                const _Tp1 &T, const Field &F)
		{
			Hom<Field, _Tp1>(F, T).image(Ap.getData(), A.getData());
		}
	};

	Element &getData() { return _a; }
	const Element &getData() const { return _a; }

private:
	Element _a;
};

/**
 * A reference-owning factory is intentional: the two butterflies and the
 * scalar diagonal consume one sequential random stream, rather than copies
 * of a random iterator starting in the same state.
 */
template<class Field>
class ConsecutiveLIFSwitchFactory {
public:
	typedef typename Field::RandIter RandIter;
	typedef typename Field::Element Element;

	explicit ConsecutiveLIFSwitchFactory(RandIter &r) : _r(r) {}

	ConsecutiveLIFSwitch<Field> makeSwitch()
	{
		Element a;
		_r.random(a);
		return ConsecutiveLIFSwitch<Field>(a);
	}

private:
	RandIter &_r;
};

} // namespace Protected

/**
 * Recover an arbitrary consecutive block of Frobenius invariant factors.
 *
 * The public indexing is decreasing Frobenius order
 *
 *       F_1, F_2, ..., F_n,       F_i = s_{n-i+1},
 *
 * where s_1 divides ... divides s_n.  A request (first,count) returns
 *
 *       F_first, ..., F_{first+count-1}.
 *
 * The caller supplies the exact minimal polynomial F_1.  The algorithm forms
 * one rank-at-most-(first-1) shifted black box, invokes Gavin Harrison's
 * existing LIF routine once for its leading count factors, and removes the
 * new factors by componentwise gcd with F_1.
 *
 * The probability argument for the shift is separate from the projection
 * probability requested from LIF.  Consequently lifProbability controls only
 * the latter; over a small field the caller should use an extension field or
 * repeat the whole shifted computation.
 */
template<class _Field, class _PolynomialRing,
	 class _MatrixDomain = MatrixDomain<_Field> >
class ConsecutiveInvariantFactors {
public:
	typedef _Field Field;
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef _PolynomialRing PolynomialRing;
	typedef typename PolynomialRing::Element Polynomial;
	typedef _MatrixDomain MatrixDom;
	typedef typename MatrixDom::OwnMatrix Matrix;
	typedef InvariantFactors<Field, PolynomialRing, MatrixDom> LIFDomain;

	ConsecutiveInvariantFactors(const Field &F, const PolynomialRing &R)
		: _F(F), _RI(_F), _R(R), _LIF(_F, _R)
	{}

	/** Integer ceiling of log base 2; used as an automatic cost heuristic. */
	static size_t denseButterflyCrossover(size_t n)
	{
		assert(n > 0);
		size_t levels = 0;
		for (size_t m = n - 1; m != 0; m >>= 1) ++levels;
		return levels;
	}

	/** Dense costs O(nh), butterfly costs O(n log n), h=first-1. */
	static ConsecutiveLIFShift automaticShift(size_t n, size_t first)
	{
		assert(n > 0);
		assert(first >= 1 && first <= n);
		return first - 1 <= denseButterflyCrossover(n)
			? ConsecutiveLIFDense : ConsecutiveLIFButterfly;
	}

	/**
	 * Return F_first,...,F_{first+count-1} in that (nonincreasing) order.
	 * Exactly one call to largestInvariantFactors is made.
	 */
	template<class Blackbox>
	std::vector<Polynomial> &consecutiveInvariantFactors(
		std::vector<Polynomial> &result,
		const Blackbox &A,
		const Polynomial &minimalPolynomial,
		size_t first,
		size_t count,
		double lifProbability = 0.99,
		ConsecutiveLIFShift shift = ConsecutiveLIFAuto)
	{
		const size_t n = A.rowdim();
		assert(n == A.coldim());
		assert(first >= 1 && first <= n);
		assert(count >= 1 && count <= n - first + 1);
		assert(0.0 < lifProbability && lifProbability < 1.0);

		// result is allowed to contain minimalPolynomial on entry.
		Polynomial minpolyA;
		_R.assign(minpolyA, minimalPolynomial);

		const size_t h = first - 1;
		if (h == 0) {
			return runLIF(result, A, minpolyA, count,
			              lifProbability);
		}

		if (shift == ConsecutiveLIFAuto)
			shift = automaticShift(n, first);

		if (shift == ConsecutiveLIFDense)
			return runDense(result, A, minpolyA, h, count,
			                lifProbability);

		assert(shift == ConsecutiveLIFButterfly);
		return runButterfly(result, A, minpolyA, h, count,
		                    lifProbability);
	}

	/** Short compatibility spelling. */
	template<class Blackbox>
	std::vector<Polynomial> &frobeniusInvariants(
		std::vector<Polynomial> &result,
		const Blackbox &A,
		const Polynomial &minimalPolynomial,
		size_t first,
		size_t count,
		double lifProbability = 0.99,
		ConsecutiveLIFShift shift = ConsecutiveLIFAuto)
	{
		return consecutiveInvariantFactors(result, A, minimalPolynomial,
		                                   first, count, lifProbability, shift);
	}

private:
	Field _F;
	RandIter _RI;
	PolynomialRing _R;
	LIFDomain _LIF;

	template<class ShiftedBlackbox>
	std::vector<Polynomial> &runLIF(
		std::vector<Polynomial> &result,
		const ShiftedBlackbox &Ashift,
		const Polynomial &minimalPolynomial,
		size_t count,
		double lifProbability)
	{
		std::vector<Polynomial> leading;
		_LIF.largestInvariantFactors(leading, Ashift, count, lifProbability);
		assert(leading.size() >= count);

		result.clear();
		result.reserve(count);
		for (size_t i = 0; i < count; ++i) {
			Polynomial f;
			_R.assign(f, leading[leading.size() - 1 - i]);
			_R.gcdin(f, minimalPolynomial);
			result.push_back(f);
		}
		return result;
	}

	template<class Blackbox>
	std::vector<Polynomial> &runDense(
		std::vector<Polynomial> &result,
		const Blackbox &A,
		const Polynomial &minimalPolynomial,
		size_t h,
		size_t count,
		double lifProbability)
	{
		const size_t n = A.rowdim();
		Matrix U(_F, n, h);
		Matrix V(_F, h, n);

		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < h; ++j) {
				Element e;
				_RI.random(e);
				U.setEntry(i, j, e);
			}

		for (size_t i = 0; i < h; ++i)
			for (size_t j = 0; j < n; ++j) {
				Element e;
				_RI.random(e);
				V.setEntry(i, j, e);
			}

		typedef Compose<Matrix, Matrix> DenseShift;
		DenseShift UV(U, V);
		/*
		 * Keep the dense term first.  Sum sends its first term the caller's
		 * BlasSubvector output, which BlasMatrix supports, and sends the
		 * second term its internal std::vector output, which the input
		 * black box supports.  Reversing these operands makes BlasMatrix
		 * receive std::vector and fails to compile.
		 */
		typedef Sum<DenseShift, Blackbox> ShiftedBlackbox;
		ShiftedBlackbox Ashift(UV, A);
		return runLIF(result, Ashift, minimalPolynomial, count,
		              lifProbability);
	}

	template<class Blackbox>
	std::vector<Polynomial> &runButterfly(
		std::vector<Polynomial> &result,
		const Blackbox &A,
		const Polynomial &minimalPolynomial,
		size_t h,
		size_t count,
		double lifProbability)
	{
		const size_t n = A.rowdim();
		typedef Protected::ConsecutiveLIFSwitch<Field> Switch;
		typedef Butterfly<Field, Switch> ButterflyMatrix;
		typedef TransposeOwner<ButterflyMatrix> ButterflyTranspose;
		typedef SparseMatrix<Field> Diagonal;

		// The proof uses one common scalar d on all h active positions.
		Element d;
		do {
			_RI.random(d);
		} while (_F.isZero(d));

		typename Switch::Factory factory(_RI);
		ButterflyMatrix U(_F, n, factory);
		ButterflyMatrix V(_F, n, factory);
		ButterflyTranspose Vt(V);

		Diagonal D(_F, n, n);
		for (size_t i = 0; i < h; ++i)
			D.setEntry(n - 1 - i, n - 1 - i, d);
		D.finalize();

		typedef Compose<Diagonal, ButterflyMatrix> DU_t;
		DU_t DU(D, U);
		typedef Compose<ButterflyTranspose, DU_t> ButterflyShift;
		ButterflyShift VtDU(Vt, DU);
		typedef Sum<Blackbox, ButterflyShift> ShiftedBlackbox;
		ShiftedBlackbox Ashift(A, VtDU);

		return runLIF(result, Ashift, minimalPolynomial, count,
		              lifProbability);
	}
};

} // namespace LinBox

#endif // __LINBOX_consecutive_invariant_factors_H

// Local Variables:
// mode: c++
// tab-width: 4
// indent-tabs-mode: t
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s