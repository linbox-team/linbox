/* linbox/algorithms/invariant-factors.h
 * Copyright (C) 2015 Gavin Harrison
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

#ifndef __LINBOX_invariant_factors_H
#define __LINBOX_invariant_factors_H

#include <list>
#include <vector>
#include <math.h>

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

#include "linbox/blackbox/block-compose.h"
#include "linbox/blackbox/fflas-csr.h"

#include "linbox/algorithms/block-coppersmith-domain.h"
#include "linbox/algorithms/block-massey-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/poly-smith-form.h"

namespace LinBox
{

template<class _Field, class _PolynomialRing, class _MatrixDomain = MatrixDomain<_Field>>
class InvariantFactors {
public:
	typedef _Field Field;
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef _MatrixDomain MatrixDom;
	typedef typename MatrixDom::OwnMatrix Matrix;

	typedef _PolynomialRing PolynomialRing;
	typedef typename PolynomialRing::Element Polynomial;
	typedef typename MatrixDomain<PolynomialRing>::OwnMatrix PolyMatrix;

	typedef PolySmithFormDomain<PolynomialRing> SmithFormDom;

protected:
	Field _F;
	PolynomialRing _R;
	SmithFormDom _SFD;

public:
	InvariantFactors(const Field &F, const PolynomialRing &R) : _F(F), _R(R), _SFD(R) {}

public:
	size_t min_block_size(size_t t, double p) const {
		size_t q = _F.cardinality();
		assert ((0.0 < p) && (p < 1.0) && (t >= 1) && (q >= 2));

		double k = (q == 2) ? 3 : 2;
		return ceil(log(k / (1 - sqrt(p)))/log(q)) + t;
	}

//protected:
	template<class Blackbox>
	void computeGenerator(
		std::vector<Matrix> &gen,
		const Blackbox &M,
		size_t b) const
	{
		RandIter RI(_F);
		RandomDenseMatrix<RandIter, Field> RDM(_F, RI);
		MatrixDom MD(_F);

		size_t n = M.rowdim();
		Matrix U(_F, b, n);
		Matrix V(_F, n, b);

		RDM.random(U);
		RDM.random(V);

		typedef BlackboxBlockContainer<Field, Blackbox, MatrixDom> Sequence;
		Sequence blockSeq(&M, _F, U, V);
		BlockCoppersmithDomain<MatrixDom, Sequence> coppersmith(MD, &blockSeq, 10);

		coppersmith.right_minpoly(gen);
	}

	template<class Blackbox>
	void computeGenerator2(
		std::vector<Matrix> &gen,
		const Blackbox &M,
		size_t b) const
	{
		RandIter RI(_F);
		RandomDenseMatrix<RandIter, Field> RDM(_F, RI);
		MatrixDom MD(_F);

		size_t n = M.rowdim();
		Matrix U(_F, b, n);
		Matrix V(_F, n, b);

		RDM.random(U);
		RDM.random(V);

		typedef BlackboxBlockContainer<Field, Blackbox, MatrixDom> Sequence;
		Sequence blockSeq(&M, _F, U, V);
		BlockMasseyDomain<Field, Sequence> coppersmith(&blockSeq, 10);

		coppersmith.left_minpoly(gen);
	}

	void convert(PolyMatrix &G, const std::vector<Matrix> &minpoly) const {
		size_t b = G.rowdim();
		for (size_t i = 0; i < b; i++) {
			for (size_t j = 0; j < b; j++) {
				std::vector<long> coeffs;
				for (size_t k = 0; k < minpoly.size(); k++) {
					long coeff;
					_F.convert(coeff, minpoly[k].getEntry(i, j));
					coeffs.push_back(coeff);
				}

				Polynomial tmp;
				_R.init(tmp, coeffs);
				G.setEntry(i, j, tmp);
			}
		}
	}

public:
	// computes the t largest invariant factors of A with probability of at least p.
	template<class Blackbox>
	std::vector<Polynomial> &largestInvariantFactors(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
		size_t b,
		int earlyTerm = 10) const
	{
		std::vector<Matrix> minpoly;
		computeGenerator(minpoly, A, b);

		PolyMatrix G(_R, b, b);
		convert(G, minpoly);

		Polynomial det;
		_SFD.detLocalX(det, G);

		_SFD.solve(lifs, G, det);

		return lifs;
	}

	template<class Blackbox>
	std::vector<Polynomial> &largestInvariantFactors2(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
		size_t b,
		int earlyTerm = 10) const
	{
		std::vector<Matrix> minpoly;
		computeGenerator2(minpoly, A, b);

		PolyMatrix G(_R, b, b);
		convert(G, minpoly);

		Polynomial det;
		_SFD.detLocalX(det, G);
		_SFD.solve(lifs, G, det);

		return lifs;
	}

	template<class Blackbox>
	std::vector<Polynomial> &largestInvariantFactors(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
		const Polynomial &mod,
		size_t b,
		int earlyTerm = 10) const
	{
		std::vector<Matrix> minpoly;
		computeGenerator(minpoly, A, b);

		PolyMatrix G(_R, b, b);
		convert(G, minpoly);

		Polynomial det;
		_SFD.solve(lifs, G, mod, false);

		return lifs;
	}

	template<class Blackbox>
	std::vector<Polynomial> &largestInvariantFactors(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
		size_t t,
		double p,
		int earlyTerm = 10) const
	{
		size_t b = min_block_size(t, p);

		std::vector<Matrix> minpoly;
		computeGenerator(minpoly, A, b);

		PolyMatrix G(_R, b, b);
		convert(G, minpoly);

		Polynomial det;
		_SFD.detLocalX(det, G);
		_SFD.solve(lifs, G, det);

		return lifs;
	}

   /// lifs is the first k invariant factors of A in nonincreasing order by degree.
	template<class Blackbox>
	std::vector<Polynomial> &frobeniusInvariants(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
      size_t k) const
   {
      // b is a fudge for now.  Todo: make a full FNF based on blockwied
      size_t b = (k == 0 ? A.rowdim() : k+2);
      largestInvariantFactors(lifs, A, b);
      std::reverse(lifs.begin(),lifs.end());
      return lifs;
   }

	template<class Blackbox>
	std::vector<Polynomial> &lifsit(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
		size_t b0,
		size_t s,
		size_t t, // number of iterations
		size_t k) const
	{
		largestInvariantFactors(lifs, A, b0);

		size_t b = s;
		Polynomial mod;
		_R.assign(mod, lifs[k]);
		for (size_t i = 0; i < t; i ++) {
			std::cout << "computing using block size: " << b << std::endl;
			if (_R.isIrreducible(mod)) {
				return lifs;
			}

			std::vector<Polynomial> part;
			largestInvariantFactors(part, A, mod, b);

			for (size_t j = 0; j < part.size() - lifs.size() + k; j++) {
				if (_R.isZero(part[j])) {
					_R.assign(part[j], mod);
				}
			}
			for (size_t j = k; j < lifs.size(); j++) {
				_R.assign(part[part.size() - lifs.size() + j], lifs[j]);
			}

			lifs = part;
			_R.assign(mod, lifs[k]);
			b *= 2;
		}

		return lifs;
	}

	void randomNonzero(const Field &F, Element &elm) const {
		typename Field::RandIter RI(F);
		do {
			RI.random(elm);
		} while (F.isZero(elm));
	}

	template<class SparseMat>
	void randomTriangular(SparseMat &T, size_t s, bool upper, bool randomDiag = false) const {
		Field F(T.field());
		typename Field::RandIter RI(F);

		for (size_t i = 0; i < T.rowdim(); i++) {
			if (randomDiag) {
				Element elm;
				randomNonzero(F, elm);
				T.setEntry(i, i, elm);
			} else {
				T.setEntry(i, i, F.one);
			}
		}

		for (size_t r = 0; r < T.rowdim() - 1; r++) {
			for (size_t k = 0; k < s; k++) {
				size_t c = (rand() % (T.coldim() - r - 1)) + r + 1;

				Element elm;
				randomNonzero(F, elm);

				if (upper) {
					T.setEntry(r, c, elm);
				} else {
					T.setEntry(c, r, elm);
				}
			}
		}

		T.finalize();
	}

	template<class Blackbox>
	std::vector<Polynomial> &precondLifs(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
		size_t b,
		int earlyTerm = 10) const
	{
		size_t n = A.rowdim();
		size_t s = std::ceil(std::log(n)) + 1;

		typedef SparseMatrix<Field, SparseMatrixFormat::CSR> SparseMat;
		SparseMat L(_F, n, n);
		SparseMat U(_F, n, n);
		randomTriangular(L, s, false, false);
		randomTriangular(U, s, true, false);

		typedef FflasCsr<Field> Csr;
		Csr FL(&L);
		Csr FU(&U);

		BlockCompose<Csr, Csr> LU(FL, FU);
		BlockCompose<Csr, BlockCompose<Csr, Csr>> ALU(A, LU);

		return largestInvariantFactors(lifs, ALU, b);
	}

	template<class Blackbox>
	bool det(
		Element &det,
		const Blackbox &A,
		size_t b = 12,
		int earlyTerm = 10) const
	{
		std::vector<Polynomial> lifs;
		precondLifs(lifs, A, b, earlyTerm);

		typename PolynomialRing::CoeffField CF(_R.getCoeffField());
		typename PolynomialRing::Coeff cdet;

		CF.assign(cdet, CF.one);
		size_t total_deg = 0;
		for (size_t i = 0; i < lifs.size(); i++) {
			total_deg += _R.deg(lifs[i]);

			typename PolynomialRing::Coeff c;
			_R.getCoeff(c, lifs[i], 0);
			CF.mulin(cdet, c);
		}

		integer idet;
		CF.convert(idet, cdet);
		_F.init(det, idet);

		return _F.isZero(det) || total_deg == A.rowdim();
	}

	template<class Blackbox>
	bool rank(
		size_t &rank,
		const Blackbox &A,
		size_t b = 12,
		size_t t = 3,
		int earlyTerm = 10) const
	{
		std::vector<Polynomial> lifs;
		precondLifs(lifs, A, b, earlyTerm);

		typename PolynomialRing::CoeffField CF(_R.getCoeffField());
		typename PolynomialRing::Coeff c;

		rank = 0;
		for (size_t i = 0; i < lifs.size(); i++) {
			rank += _R.deg(lifs[i]);

			_R.getCoeff(c, lifs[i], 0);
			if (CF.isZero(c)) {
				rank -= 1;
			}
		}

		_R.getCoeff(c, lifs[t], 0);
		return _R.isOne(lifs[t]) || (_R.deg(lifs[t]) == 1 && CF.isZero(c));
	}
};

}

#endif //__LINBOX_invariant_factors_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
