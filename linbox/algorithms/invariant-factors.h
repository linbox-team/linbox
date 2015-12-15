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
	typedef Givaro::Poly1Dom<Field,Givaro::Dense> PolyDom;
	typedef GivaroPoly<PolyDom> PolyRing;
	typedef typename PolyRing::Element PolyElement;
	typedef MatrixDomain<PolyRing> PolyMatDom;
	typedef typename PolyMatDom::OwnMatrix PolyBlock;
	typedef RandomDenseMatrix<RandIter, Field> RandomMatrix;
	typedef SmithFormKannanBachemDomain<PolyMatDom> SmithKbDomain;
	typedef BlackboxBlockContainer<Field,Blackbox> Sequence;
	typedef BlockCoppersmithDomain<Domain, Sequence> CoppersmithDomain;
	
protected:
	Domain _MD;
	Field _F;
	RandIter _RI;
	RandomMatrix _RDM;
	PolyDom _PD;
	PolyRing _R;
	PolyMatDom _PMD;
	SmithKbDomain _SFKB;
	
public:
	
	InvariantFactors(Field &F) :
		_MD(F),
		_F(F),
		_RI(F),
		_RDM(F, _RI),
		_PD(F, "x"),
		_R(_PD),
		_PMD(_R),
		_SFKB(_PMD)
	{
	}
	
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
		_PD.init(temp, d-1);
		
		size_t b = MM.rowdim();
		for (uint32_t i = 0; i < b; i++) {
			for (uint32_t j = 0; j < b; j++) {
				for (uint32_t k = 0; k < d; k++) {
					_PD.setEntry(temp, gen[k].getEntry(i,j), k);
				}
				MM.setEntry(i,j,temp);
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
	void computeFactors(
		PolyRingVector &diag,
		const Blackbox &M,
		size_t b,
		int earlyTerm = 10)
	{
		std::vector<Block> gen;
		computeGenerator(gen, M, b, earlyTerm);
		
		PolyBlock MM(_R, b, b);
		convertSequenceToPolyMatrix(MM, gen);
		
		computeSmithForm(diag, MM, b);
	}
};

}

#endif //__LINBOX_invariant_factors_H
