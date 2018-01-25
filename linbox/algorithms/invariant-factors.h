/* linbox/algorithms/coppersmith-invariant-factors.h
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

#include "linbox/ring/polynomial-local-x.h"
#include "linbox/algorithms/poly-smith-form-local-x.h"
#include "linbox/algorithms/smith-form-local.h"

#include "linbox/algorithms/block-coppersmith-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/matrix/random-matrix.h"

namespace LinBox
{

template<class _Field, class _PolynomialRing, class _QuotientRing>
class InvariantFactors {
public:
	typedef _Field Field;
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef MatrixDomain<Field> MatrixDom;
	typedef typename MatrixDom::OwnMatrix Matrix;
	
	typedef _PolynomialRing PolynomialRing;
	typedef typename PolynomialRing::Element Polynomial;
	typedef MatrixDomain<PolynomialRing> PolyMatrixDom;
	typedef typename PolyMatrixDom::OwnMatrix PolyMatrix;
	
	typedef _QuotientRing QuotientRing;
	typedef typename QuotientRing::Element QPolynomial;
	typedef MatrixDomain<QuotientRing> QuotMatrixDom;
	typedef typename QuotMatrixDom::OwnMatrix QuotMatrix;
	typedef SmithFormLocal<QuotientRing> SmithFormLocalDom;
	
	typedef PolynomialLocalX<PolynomialRing> LocalRing;
	typedef MatrixDomain<LocalRing> LocalMatrixDom;
	typedef typename LocalMatrixDom::OwnMatrix LocalMatrix;
	typedef PolySmithFormLocalXDomain<PolynomialRing> LocalSmithFormDom;
		
protected:
	Field _F;
	PolynomialRing _R;
	
public:
	InvariantFactors(Field &F, PolynomialRing &R) : _F(F), _R(R) {}

//protected:

	size_t min_block_size(size_t t, double p) const {
		size_t q = _F.cardinality();
		assert (0.0 < p < 1.0 && t >= 1 && q >= 2);
		
		double k = (q == 2) ? 3 : 2;
		return ceil(log(k / (1 - sqrt(p)))/log(q)) + t;
	}

	template<class Blackbox>
	void computeGenerator(
		std::vector<Matrix> &gen,
		const Blackbox &M,
		size_t b,
		int earlyTerm) const
	{
		RandIter RI(_F);
		RandomDenseMatrix<RandIter, Field> RDM(_F, RI);
		MatrixDom MD(_F);
		
		size_t n = M.rowdim();
		Matrix U(_F, b, n);
		Matrix V(_F, n, b);
		
		RDM.random(U);
		RDM.random(V);
		
		typedef BlackboxBlockContainer<Field, Blackbox> Sequence;
		Sequence blockSeq(&M, _F, U, V);
		BlockCoppersmithDomain<MatrixDom, Sequence> coppersmith(MD, &blockSeq, earlyTerm);
		
		coppersmith.right_minpoly(gen);
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
	
	size_t detLimit(const PolyMatrix &M, size_t dim) const {
		size_t limit1 = 0;
		for (size_t i = 0; i < M.rowdim(); i++) {
			size_t max_degree = 0;
			for (size_t j = 0; j < M.coldim(); j++) {
				size_t deg = _R.deg(M.getEntry(i, j));
				if (deg > max_degree) {
					max_degree = deg;
				}
			}
			
			limit1 += max_degree;
		}
		
		size_t limit2 = 0;
		for (size_t i = 0; i < M.coldim(); i++) {
			size_t max_degree = 0;
			for (size_t j = 0; j < M.rowdim(); j++) {
				size_t deg = _R.deg(M.getEntry(j, i));
				if (deg > max_degree) {
					max_degree = deg;
				}
			}
			
			limit2 += max_degree;
		}
				
		return std::min(std::min(limit1, limit2), dim) + 1;
	}
	
	void localX(Polynomial &det, const PolyMatrix &M, size_t exponent) const {
		LocalSmithFormDom SFD(_R, exponent);
		LocalRing L(_R, exponent);
				
		LocalMatrix G(M, L);
		
		SFD.solveDet(det, G);
		NTL::MakeMonic(det);
	}
	
	void local(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &f,
		long multiplicity) const {
	
		SmithFormLocalDom SFD;
		
		Polynomial modulus;
		_R.pow(modulus, f, multiplicity);
		
		QuotientRing QR(_R, modulus);
		
		QuotMatrix QM(M, QR);
		
		std::list<QPolynomial> L;
		SFD(L, QM, QR);
		
		Hom<PolynomialRing, QuotientRing> hom(_R, QR);
		
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
	
	void factoredLocal(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &sf_factor,
		long multiplicity) const {
	
		std::vector<std::pair<Polynomial, long>> factors;
		_R.factor(factors, sf_factor);
		
		for (size_t i = 0; i < factors.size(); i++) {
			local(result, M, factors[i].first, factors[i].second * multiplicity);
		}
	}
	
	void factoredLocal(
		std::vector<Polynomial> &result,
		const PolyMatrix &M,
		const Polynomial &det) const {
	
		std::vector<std::pair<Polynomial, long>> factors;
		
		_R.squareFree(factors, det);
		
		result.clear();
		for (size_t i = 0; i < M.rowdim(); i++) {
			result.push_back(_R.one);
		}
		
		for (size_t i = 0; i < factors.size(); i++) {
			if (factors[i].second == 1) {
				_R.mulin(result[result.size() - 1], factors[i].first);
			} else {
				factoredLocal(result, M, factors[i].first, factors[i].second);
			}
		}
	}
	
public:
	// computes the t largest invariant factors of A with probability of at least p.
	template<class Blackbox>
	std::vector<Polynomial> &largestInvariantFactors(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
		size_t t,
		double p,
		int earlyTerm = 10) const {
	
		size_t b = min_block_size(t, p);
	
		std::vector<Matrix> minpoly;
		computeGenerator(minpoly, A, b, earlyTerm);
		
		PolyMatrix G(_R, b, b);
		convert(G, minpoly);
		
		Polynomial det;
		size_t limit = detLimit(G, A.rowdim());
		localX(det, G, limit);
		
		factoredLocal(lifs, G, det);		
		return lifs;
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
