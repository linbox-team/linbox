/* tests/test-consecutive-invariant-factors.C
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * ========LICENCE========
 */

#include "linbox/linbox-config.h"

#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <vector>

#include "linbox/algorithms/consecutive-invariant-factors.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/ring/modular.h"
#include "linbox/ring/ntl.h"
#include "linbox/util/commentator.h"

using namespace LinBox;

template<class PolyRing>
static typename PolyRing::Element xPower(const PolyRing &R, size_t degree)
{
	typedef typename PolyRing::Element Polynomial;
	typedef typename PolyRing::Coeff Coeff;

	Polynomial f(degree + 1);
	for (size_t i = 0; i < degree; ++i)
		R.setCoeff(f, i, (Coeff)0);
	R.setCoeff(f, degree, (Coeff)1);
	return f;
}

template<class PolyRing>
static bool equalLists(
	const PolyRing &R,
	const std::vector<typename PolyRing::Element> &a,
	const std::vector<typename PolyRing::Element> &b)
{
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); ++i)
		if (!R.areEqual(a[i], b[i])) return false;
	return true;
}

template<class PolyRing>
static void writeList(
	std::ostream &out,
	const PolyRing &R,
	const std::vector<typename PolyRing::Element> &fs)
{
	out << "[";
	for (size_t i = 0; i < fs.size(); ++i) {
		if (i != 0) out << ", ";
		R.write(out, fs[i]);
	}
	out << "]";
}

template<class Polynomial>
static std::vector<Polynomial> profileSlice(
	const std::vector<Polynomial> &profile,
	size_t first,
	size_t count)
{
	return std::vector<Polynomial>(
		profile.begin() + static_cast<std::ptrdiff_t>(first - 1),
		profile.begin() + static_cast<std::ptrdiff_t>(first - 1 + count));
}

/*
 * Both stages are probabilistic, and Gavin's computeGenerator creates its own
 * random iterator.  Accept an attempt only after exact comparison with the
 * known Smith profile; a persistent implementation error fails every attempt.
 */
template<class Consecutive, class PolyRing, class Blackbox>
static bool runCase(
	Consecutive &CIF,
	const PolyRing &R,
	const Blackbox &A,
	const typename PolyRing::Element &minpolyA,
	size_t first,
	const std::vector<typename PolyRing::Element> &expected,
	const char *description)
{
	typedef typename PolyRing::Element Polynomial;
	const size_t maxAttempts = 8;
	const double lifProbability = 0.999;
	std::vector<Polynomial> result;

	for (size_t attempt = 0; attempt < maxAttempts; ++attempt) {
		CIF.consecutiveInvariantFactors(
			result, A, minpolyA, first, expected.size(), lifProbability);
		if (equalLists(R, result, expected)) return true;
	}

	std::ostream &report = commentator().report(
		Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
	report << description << " failed after " << maxAttempts
	       << " attempts\nexpected: ";
	writeList(report, R, expected);
	report << "\ncomputed: ";
	writeList(report, R, result);
	report << std::endl;
	return false;
}

int main()
{
	typedef Givaro::Modular<double> Field;
	typedef NTL_zz_pX PolyRing;
	typedef PolyRing::Element Polynomial;
	typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;
	typedef ConsecutiveInvariantFactors<Field, PolyRing> Consecutive;

	const size_t q = 65521U;
	Field F(q);
	PolyRing R(static_cast<int>(q));
	std::srand(0x13579U);

	/*
	 * A = J_3^5 direct-sum J_2^5 direct-sum J_1^5, so n=30 and
	 *
	 * F_1,...,F_5   = x^3,
	 * F_6,...,F_10  = x^2,
	 * F_11,...,F_15 = x,
	 * F_16,...,F_30 = 1.
	 */
	const size_t n = 30;
	SparseMat A(F);
	A.resize(n, n);
	size_t offset = 0;

	for (size_t block = 0; block < 5; ++block) {
		for (size_t i = 1; i < 3; ++i)
			A.setEntry(offset + i, offset + i - 1, F.one);
		offset += 3;
	}
	for (size_t block = 0; block < 5; ++block) {
		A.setEntry(offset + 1, offset, F.one);
		offset += 2;
	}
	offset += 5; // five J_1 blocks
	if (offset != n) return -1;
	A.finalize();

	const Polynomial x1 = xPower(R, 1);
	const Polynomial x2 = xPower(R, 2);
	const Polynomial x3 = xPower(R, 3);
	Polynomial one;
	R.assign(one, R.one);

	std::vector<Polynomial> profile;
	profile.insert(profile.end(), 5, x3);
	profile.insert(profile.end(), 5, x2);
	profile.insert(profile.end(), 5, x1);
	profile.insert(profile.end(), 15, one);

	Consecutive CIF(F, R);
	bool pass = true;
	commentator().start("Consecutive invariant factors test suite",
	                    "consecutive-invariant-factors");

	// h=0: no shift, one ordinary LIF call.
	pass = runCase(CIF, R, A, x3, 1, profileSlice(profile, 1, 3),
	               "unshifted leading block F_1,...,F_3") && pass;

	/*
	 * ceil(log_2(30))=5 and h=3, so Auto selects the dense shift.
	 * This is the exact requested block F_4,...,F_13:
	 * x^3,x^3, x^2,x^2,x^2,x^2,x^2, x,x,x.
	 */
	pass = runCase(CIF, R, A, x3, 4, profileSlice(profile, 4, 10),
	               "dense consecutive block F_4,...,F_13") && pass;

	// h=6>5: Auto selects the corrected pass/swap butterfly shift.
	pass = runCase(CIF, R, A, x3, 7, profileSlice(profile, 7, 7),
	               "butterfly consecutive block F_7,...,F_13") && pass;

	// An interior block that crosses from nonunits into unit factors.
	pass = runCase(CIF, R, A, x3, 13, profileSlice(profile, 13, 5),
	               "interior block F_13,...,F_17 containing units") && pass;

	pass = pass
		&& Consecutive::denseButterflyCrossover(n) == 5
		&& Consecutive::automaticShift(n, 4) == ConsecutiveLIFDense
		&& Consecutive::automaticShift(n, 7) == ConsecutiveLIFButterfly;

	commentator().stop(MSG_STATUS(pass), (const char*)0,
	                   "Consecutive invariant factors test suite");
	return pass ? 0 : -1;
}