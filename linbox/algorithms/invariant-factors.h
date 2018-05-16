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

#include "linbox/algorithms/poly-smith-form.h"
#include "linbox/algorithms/block-coppersmith-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/blackbox-block-container-smmx.h"
#include "linbox/algorithms/blackbox-block-container-spmv.h"
#include "linbox/matrix/random-matrix.h"

namespace LinBox
{

template<class _Field, class _PolynomialRing>
class InvariantFactors {
public:
	typedef _Field Field;
	typedef typename Field::Element Element;
	typedef typename Field::RandIter RandIter;
	typedef MatrixDomain<Field> MatrixDom;
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

//protected:

	size_t min_block_size(size_t t, double p) const {
		size_t q = _F.cardinality();
		assert (0.0 < p < 1.0 && t >= 1 && q >= 2);
		
		double k = (q == 2) ? 3 : 2;
		return ceil(log(k / (1 - sqrt(p)))/log(q)) + t;
	}
	
	template<class Sequence>
	void computeGenerator(
		std::vector<Matrix> &gen,
		Sequence &blockSeq,
		int earlyTerm = 10) const
	{
		MatrixDom MD(_F);
		BlockCoppersmithDomain<MatrixDom, Sequence> coppersmith(MD, &blockSeq, earlyTerm);
		coppersmith.right_minpoly(gen);
		std::cout << "spmm:" << blockSeq.spmmTime() << "\t";
	}

	template<class Blackbox>
	void computeGenerator(
		std::vector<Matrix> &gen,
		const Blackbox &M,
		size_t b,
		int earlyTerm = 10) const
	{
		RandIter RI(_F);
		RandomDenseMatrix<RandIter, Field> RDM(_F, RI);
		MatrixDom MD(_F);
		
		size_t n = M.rowdim();
		Matrix U(_F, b, n);
		Matrix V(_F, n, b);
		
		RDM.random(U);
		RDM.random(V);
		
		//typedef BlackboxBlockContainer<Field, Blackbox> Sequence;
		typedef BlackboxBlockContainerSmmx<Field, Blackbox> Sequence;
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
		_SFD.detLocalX(det, G);
		_SFD.solve(lifs, G, det);	
		
		return lifs;
	}
	
	// computes the t largest invariant factors of A with probability of at least p.
	template<class Blackbox>
	Element &det(
		Element &d, // det(A)
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
		_SFD.detLocalX(det, G);
		
		// get the constant coefficient of det and convert it to type Element
		typename PolynomialRing::Coeff det0;
		_R.getCoeff(det0, det, 0);
		
		integer tmp;
		_R.getCoeffField().convert(tmp, det0);
		_F.init(d, tmp);
		
		return d;
	}
	
	// computes the t largest invariant factors of A with probability of at least p.
	template<class Blackbox>
	std::vector<Polynomial> &largestInvariantFactors(
		std::vector<Polynomial> &lifs,
		const Blackbox &A,
		const Blackbox &PreR,
		size_t t,
		double p,
		int earlyTerm = 10) const {
	
		size_t b = min_block_size(t, p);
	
		std::vector<Matrix> minpoly;
		computeGenerator(minpoly, A, PreR, b, earlyTerm);
		
		PolyMatrix G(_R, b, b);
		convert(G, minpoly);
		
		Polynomial det;
		_SFD.detLocalX(det, G);
		_SFD.solve(lifs, G, det);	
		
		return lifs;
	}
	
	// computes the t largest invariant factors of A with probability of at least p.
	template<class Blackbox>
	Element &det(
		Element &d, // det(A)
		const Blackbox &A,
		const Blackbox &PreR,
		size_t t,
		double p,
		int earlyTerm = 10) const {
	
		size_t b = min_block_size(t, p);
	
		std::vector<Matrix> minpoly;
		computeGenerator(minpoly, A, PreR, b, earlyTerm);
		
		PolyMatrix G(_R, b, b);
		convert(G, minpoly);
		
		Polynomial det;
		_SFD.detLocalX(det, G);
		
		// get the constant coefficient of det and convert it to type Element
		typename PolynomialRing::Coeff det0;
		_R.getCoeff(det0, det, 0);
		
		integer tmp;
		_R.getCoeffField().convert(tmp, det0);
		_F.init(d, tmp);
		
		return d;
	}
	
	// computes the t largest invariant factors of A with probability of at least p.
	template<class Blackbox>
	std::vector<Polynomial> &largestInvariantFactors(
		std::vector<Polynomial> &lifs,
		const Blackbox &PreL,
		const Blackbox &A,
		const Blackbox &PreR,
		size_t t,
		double p,
		int earlyTerm = 10) const {
	
		size_t b;
		if (t > 1000) t = b = t-1000; else b = min_block_size(t, p);
	
		std::vector<Matrix> minpoly;
		computeGenerator(minpoly, PreL, A, PreR, b, earlyTerm);
		
		PolyMatrix G(_R, b, b);
		convert(G, minpoly);
		
		Polynomial det;
		_SFD.detLocalX(det, G);
		_SFD.solve(lifs, G, det);	
		
		return lifs;
	}
	
	// computes the t largest invariant factors of A with probability of at least p.
	template<class Blackbox>
	Element &det(
		Element &d, // det(A)
		const Blackbox &PreL,
		const Blackbox &A,
		const Blackbox &PreR,
		size_t t,
		double p,
		int earlyTerm = 10) const {
	
		size_t b = min_block_size(t, p);
	
		std::vector<Matrix> minpoly;
		computeGenerator(minpoly, PreL, A, PreR, b, earlyTerm);
		
		PolyMatrix G(_R, b, b);
		convert(G, minpoly);
		
		Polynomial det;
		_SFD.detLocalX(det, G);
		
		// get the constant coefficient of det and convert it to type Element
		typename PolynomialRing::Coeff det0;
		_R.getCoeff(det0, det, 0);
		
		integer tmp;
		_R.getCoeffField().convert(tmp, det0);
		_F.init(d, tmp);
		
		return d;
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
