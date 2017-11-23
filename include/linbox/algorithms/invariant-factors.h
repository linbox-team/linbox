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

#include "linbox/algorithms/block-coppersmith-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/matrix/random-matrix.h"

#include <givaro/givpoly1.h>
#include <linbox/ring/givaro-poly.h>
#include <linbox/algorithms/smith-form-kannan-bachem.h>
#include <linbox/algorithms/smith-form-iliopoulos2.h>

#include <givaro/givtimer.h>

namespace LinBox
{

template<class _Field, class _Blackbox>
class InvariantFactors {
public:
	typedef _Field Field;
	typedef _Blackbox Blackbox;
	typedef MatrixDomain<Field> Domain;
	typedef typename Domain::OwnMatrix Block;
	typedef typename Field::RandIter RandIter;
	typedef GivaroPoly<Field,Givaro::Dense> PolyRing;
	typedef typename PolyRing::Element PolyElement;
	typedef MatrixDomain<PolyRing> PolyMatDom;
	typedef typename PolyMatDom::OwnMatrix PolyBlock;
	typedef RandomDenseMatrix<RandIter, Field> RandomMatrix;
	typedef SmithFormKannanBachemDomain<PolyMatDom> SmithKbDomain;
	typedef BlackboxBlockContainer<Field,Blackbox> Sequence;
	typedef BlockCoppersmithDomain<Domain, Sequence> CoppersmithDomain;
	typedef IliopoulosDomain<PolyRing> IliopoulosDom;
	
protected:
	Domain _MD;
	Field _F;
	RandIter _RI;
	RandomMatrix _RDM;
	PolyRing _R;
	PolyMatDom _PMD;
	SmithKbDomain _SFKB;
	IliopoulosDom _SFI;
	Givaro::Timer timer;
	
public:
	InvariantFactors(Field &F, PolyRing &R) :
		_MD(F),
		_F(F),
		_RI(F),
		_RDM(F, _RI),
		_R(R),
		_PMD(R),
		_SFKB(_PMD),
		_SFI(R)
	{
	}

//protected:
	void computeGenerator(
		std::vector<Block> &gen,
		const Blackbox &M,
		size_t b,
		int earlyTerm)
	{
		size_t n = M.rowdim();
		Block U(_F, b, n);
		Block V(_F, n, b);
		
		_RDM.random(U);
		_RDM.random(V);
		
		Sequence blockSeq(&M, _F, U, V);
		CoppersmithDomain coppersmith(_MD, &blockSeq, earlyTerm);
		
		coppersmith.right_minpoly(gen);
	}
	
	void convertSequenceToPolyMatrix(
		PolyBlock &MM,
		const std::vector<Block> &gen)
	{
		PolyElement temp;
		size_t d = gen.size();
		_R.domain().init(temp, Givaro::Degree(d-1));
		
		size_t b = MM.rowdim();
		for (uint32_t i = 0; i < b; i++) {
			for (uint32_t j = 0; j < b; j++) {
				for (uint32_t k = 0; k < d; k++) {
					_R.domain().setEntry(temp, gen[k].getEntry(i,j), Givaro::Degree(k));
				}
				MM.setEntry(i,j,temp);
			}
		}
	}
	
	void modMatrix(PolyBlock &M, const PolyElement &d)
	{
		size_t n = M.coldim();
		
		for (size_t i = 0; i < n; i++) {
			for (size_t j = 0; j < n; j++) {
				PolyElement tmp;
				M.getEntry(tmp, i, j);
				_R.modin(tmp, d);
				M.setEntry(i, j, tmp);
			}
		}
	}
	
	template <class PolyRingVector>
	void computeSmithForm(PolyRingVector &diag, const PolyBlock &M, size_t b)
	{
		diag.resize(b);
		_SFKB.solve(diag, M);
		
		for (uint32_t i = 0; i < diag.size(); i++) {
			_R.normalizeIn(diag[i]);
		}
	}
	
	template <class PolyRingVector>
	void computeSmithForm(
		PolyRingVector &diag,
		const PolyBlock &M,
		const PolyElement &d,
		size_t b)
	{
		diag.resize(b);
		_SFI.smithForm(diag, M, d);
	}
	
	template <class PolyRingVector>
	void computeFactors(
		PolyRingVector &diag,
		const Blackbox &M,
		size_t b,
		int earlyTerm = 10)
	{
		timer.clear();
		timer.start();
		std::vector<Block> gen;
		computeGenerator(gen, M, b, earlyTerm);
		timer.stop();
		
		std::cout << "Time to compute first generator: "  << timer.usertime();
		std::cout << std::endl;
		
		timer.clear();
		
		timer.start();
		PolyBlock MM(_R, b, b);
		convertSequenceToPolyMatrix(MM, gen);
		timer.stop();
		
		//std::cout << "Time to convert first matrix to poly matrix: "  << timer.usertime();
		//std::cout << std::endl;
		
		timer.clear();
		
		timer.start();
		computeSmithForm(diag, MM, b);
		timer.stop();
		
		std::cout << "Time to compute smith form kb: " << timer.usertime();
		std::cout << std::endl;
	}
	
	template <class PolyRingVector>
	void computeFactors(
		PolyRingVector &diag,
		const Blackbox &M,
		const PolyElement &d,
		size_t b,
		int earlyTerm = 10)
	{
		timer.clear();
		timer.start();
		std::vector<Block> gen;
		computeGenerator(gen, M, b, earlyTerm);
		timer.stop();
		
		std::cout << "Time to compute second generator: " << timer.usertime();
		std::cout << std::endl;
		
		timer.clear();
		timer.start();
		PolyBlock MM(_R, b, b);
		convertSequenceToPolyMatrix(MM, gen);
		timer.stop();
		
		//std::cout << "Time to convert second matrix to poly matrix: " << timer.usertime();
		//std::cout << std::endl;
		
		timer.clear();
		timer.start();
		modMatrix(MM, d);
		computeSmithForm(diag, MM, d, b);
		timer.stop();
		
		std::cout << "Time to compute smith form iliopoulos: " << timer.usertime();
		std::cout << std::endl;
	}

public:
	template <class PolyRingVector>
	void solve(
		PolyRingVector &diag,
		const Blackbox &M,
		size_t b,
		size_t b1,
		size_t r,
		int earlyTerm = 10)
	{
		// Compute first b1 factors
		PolyRingVector partialResult;
		computeFactors(partialResult, M, b1, earlyTerm);
		
		// get r-th factor
		PolyElement d;
		_R.assign(d, partialResult[b1 - r]);
		
		// Compute factors mod r-th factor
		computeFactors(diag, M, d, b, earlyTerm);
		
		// Fill in zeros before r-th factor with r-th factor
		for (size_t i = 0; i < b - r; i++) {
			if (_R.isZero(diag[i])) {
				_R.assign(diag[i], d);
			}
		}
		
		// Fill in remaining factors with original values
		for (size_t i = b - r; i < b; i++) {
			diag[i] = partialResult[i - b + b1];
		}
	}
};

}

#endif //__LINBOX_invariant_factors_H
