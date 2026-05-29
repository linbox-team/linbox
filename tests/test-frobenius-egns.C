/* tests/test-frobenius-egns.C
 * Copyright (C) 2026 Omesh Dhar Dwivedi
 * Written by Omesh Dhar Dwivedi <odd23@drexel.edu>
 *
 * Quick accuracy test for FrobeniusLargeRank (EGNS).
 * Builds flat Jordan matrices with known FNF, runs all three algorithms,
 * compares each against the ground-truth invariants.
 *
 * No dependency on test-frobenius-suite.h.
 * Include order: algorithm headers before NTL ring wrapper.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

// Algorithm headers first — they pull in linbox/ring/ntl.h (via toeplitz.h),
// which forward-declares UnparametricRandIter before the lzz_pX wrappers
// try to specialize it.
#include "linbox/algorithms/frobenius-large-generic.h"
#include "linbox/algorithms/frobenius-large-rank.h"

#include "linbox/ring/ntl/ntl-lzz_pX.h"
#include "linbox/matrix/sparse-matrix.h"

using namespace LinBox;

typedef NTL_zz_pX                                    PolynomialRing;
typedef PolynomialRing::CoeffField                   Field;
typedef PolynomialRing::Coeff                        Coeff;
typedef PolynomialRing::Element                      Polynomial;
typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;


// ---- Normalise: strip trailing 1s from any algorithm's output. ----
static std::vector<Polynomial>
stripOnes(const std::vector<Polynomial> &fs, const PolynomialRing &R)
{
	std::vector<Polynomial> out = fs;
	while (!out.empty() && R.isOne(out.back()))
		out.pop_back();
	return out;
}

// ---- Compare two invariant-factor lists element-by-element. ----
static bool sameInvariants(
	const std::vector<Polynomial> &a,
	const std::vector<Polynomial> &b,
	const PolynomialRing           &R)
{
	if (a.size() != b.size()) return false;
	for (size_t i = 0; i < a.size(); i++)
		if (!R.areEqual(a[i], b[i])) return false;
	return true;
}

// ---- Build (x - eigenval)^exp as a Polynomial. ----
static Polynomial makeXminusLambdaPow(
	const PolynomialRing &R,
	const Field          &F,
	int                   eigenval,
	size_t                exp)
{
	// (x - lambda):  coefficient[0] = -lambda,  coefficient[1] = 1
	Coeff c_one, c_neg;
	F.assign(c_one, F.one);
	F.init(c_neg, -eigenval);          // -lambda mod p

	Polynomial base;
	R.assign(base, R.zero);
	R.setCoeff(base, 0, c_neg);
	R.setCoeff(base, 1, c_one);

	Polynomial result;
	R.assign(result, R.one);
	for (size_t i = 0; i < exp; i++) R.mulin(result, base);
	return result;
}

// ---- Build a flat Jordan matrix and its expected invariant list. ----
//
// blockSizes: sizes in ANY order — the function sorts them (descending)
//             and lays them out on the diagonal.
// Expected invariant factors = one per block, sorted largest first:
//   (x - eigenval)^blockSizes[0],  (x - eigenval)^blockSizes[1], ...
static void buildFlatMatrix(
	SparseMat                &A,
	std::vector<Polynomial>  &expected,
	const Field              &F,
	const PolynomialRing     &R,
	std::vector<size_t>       blockSizes,   // passed by value, sorted here
	int                       eigenval)
{
	std::sort(blockSizes.begin(), blockSizes.end(), std::greater<size_t>());

	Field::Element lam, one;
	F.init(lam, eigenval);
	F.assign(one, F.one);

	size_t row = 0;
	for (size_t bsz : blockSizes) {
		for (size_t j = 0; j < bsz; j++) {
			if (!F.isZero(lam))             A.setEntry(row+j, row+j,   lam);
			if (j + 1 < bsz)               A.setEntry(row+j, row+j+1, one);
		}
		row += bsz;
	}
	A.finalize();

	expected.clear();
	for (size_t bsz : blockSizes)
		expected.push_back(makeXminusLambdaPow(R, F, eigenval, bsz));
}

// ---- Run one algorithm, print result, return pass/fail. ----
template<class Alg>
static bool runOne(
	Alg                            &alg,
	const std::string              &algName,
	const SparseMat                &A,
	const std::vector<Polynomial>  &expected,
	const PolynomialRing           &R)
{
	std::vector<Polynomial> fs;
	alg.frobeniusInvariants(fs, A, 0);
	std::vector<Polynomial> got = stripOnes(fs, R);

	bool pass = sameInvariants(got, expected, R);
	std::cout << "  " << (pass ? "[PASS] " : "[FAIL] ") << algName << "\n";

	if (!pass) {
		std::cout << "    expected (" << expected.size() << "): ";
		for (auto &f : expected) { R.write(std::cout, f); std::cout << " "; }
		std::cout << "\n";
		std::cout << "    got      (" << got.size() << "): ";
		for (auto &f : got)      { R.write(std::cout, f); std::cout << " "; }
		std::cout << "\n";
	}
	return pass;
}


int main()
{
	// Construct Field first then PolynomialRing from it.
	// Passing long to PolynomialRing(integer) is ambiguous between
	// Givaro::Integer's int32_t and int64_t constructors.
	// Matches test-frobenius-suite.C pattern: Field F(p,e); PolyRing R(F).
	int prime = 65521;
	Field F(prime);
	PolynomialRing R(F);

	FrobeniusLarge<PolynomialRing>          FT(R);
	FrobeniusLargeButterfly<PolynomialRing> FB(R);
	FrobeniusLargeRank<PolynomialRing>      FR(R);

	int total = 0, passed = 0;
	auto runCase = [&](const std::string &desc,
	                   SparseMat         &A,
	                   std::vector<Polynomial> &expected)
	{
		std::cout << "--- " << desc << " ---\n";
		bool p1 = runOne(FT, "Toeplitz ", A, expected, R);
		bool p2 = runOne(FB, "Butterfly", A, expected, R);
		bool p3 = runOne(FR, "EGNS     ", A, expected, R);
		total  += 3;
		passed += p1 + p2 + p3;
		std::cout << "\n";
	};

	// ---- Case 1: n=6, eigenvalue 0, blocks (2,2,1,1) ----
	// Expected: [x^2, x^2, x, x]
	{
		SparseMat A(F, 6, 6);
		std::vector<Polynomial> expected;
		buildFlatMatrix(A, expected, F, R, {2, 2, 1, 1}, 0);
		runCase("n=6  eigen=0  blocks=(2,2,1,1)", A, expected);
	}

	// ---- Case 2: n=6, eigenvalue 1, blocks (2,2,1,1) ----
	// Expected: [(x-1)^2, (x-1)^2, (x-1), (x-1)]
	{
		SparseMat A(F, 6, 6);
		std::vector<Polynomial> expected;
		buildFlatMatrix(A, expected, F, R, {2, 2, 1, 1}, 1);
		runCase("n=6  eigen=1  blocks=(2,2,1,1)", A, expected);
	}

	// ---- Case 3: n=9, eigenvalue 2, blocks (3,3,2,1) ----
	// Expected: [(x-2)^3, (x-2)^3, (x-2)^2, (x-2)]
	{
		SparseMat A(F, 9, 9);
		std::vector<Polynomial> expected;
		buildFlatMatrix(A, expected, F, R, {3, 3, 2, 1}, 2);
		runCase("n=9  eigen=2  blocks=(3,3,2,1)", A, expected);
	}

	// ---- Case 4: n=12, eigenvalue 0, blocks (3,3,2,2,1,1) ----
	// Expected: [x^3, x^3, x^2, x^2, x, x]
	{
		SparseMat A(F, 12, 12);
		std::vector<Polynomial> expected;
		buildFlatMatrix(A, expected, F, R, {3, 3, 2, 2, 1, 1}, 0);
		runCase("n=12 eigen=0  blocks=(3,3,2,2,1,1)", A, expected);
	}

	std::cout << "==============================\n";
	std::cout << "Result: " << passed << "/" << total << " passed\n";
	return (passed == total) ? 0 : 1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
//
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s