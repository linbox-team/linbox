/* linbox/algorithms/coppersmith-invariant-factors.h
 * Copyright (C) 2014 Alex Stachnik
 *
 * Written by Alex Stachnik <stachnik@udel.edu>
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

template<class Field_,class Blackbox_,class Field2_=Field_>
class InvariantFactors {
public:
	typedef Field_ Field;
	typedef Blackbox_ Blackbox;
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
		
	typedef BlackboxBlockContainer<Field,Blackbox> BBC;
	typedef BlockCoppersmithDomain<MatrixDomain<Field>,BBC> BCD;
	
protected:
	Domain MD_;
	Field F_;
	RandIter RI;
	RandomMatrix RDM;
	PolyDom PD;
	PolyRing R;
	PolyMatDom PMD;
	SmithKbDomain SFKB;
	
public:
	
	InvariantFactors(Field& F) :
		MD_(F),
		F_(F),
		RI(F),
		RDM(F, RI),
		PD(F, "x"),
		R(PD),
		PMD(R),
		SFKB(PMD)
	{
	}
	
	template <class PolyRingVector>
	size_t computeFactors(
		PolyRingVector& diag,
		const Blackbox& M,
		size_t b,
		int earlyTerm=10)
	{
		size_t n = M.rowdim();
		Block U(F_, b, n);
		Block V(F_, n, b);
		
		RDM.random(U);
		RDM.random(V);
		
		BBC blockSeq(&M, F_, U, V);
		BCD coppersmith(MD_, &blockSeq, earlyTerm);
		
		std::vector<size_t> deg;
		std::vector<typename MatrixDomain<Field2_>::OwnMatrix > gen;
		deg=coppersmith.right_minpoly(gen);
		
		size_t d=gen.size();
		PolyBlock MM(R,b,b);
		PolyElement temp;
		PD.init(temp, d-1);
		
		for (uint32_t i = 0; i < b; ++i) {
			for (uint32_t j = 0; j < b; ++j) {
				for (uint32_t k = 0; k < d; ++k) {
					PD.setEntry(temp,gen[k].getEntry(i,j),k);
				}
				MM.setEntry(i,j,temp);
			}
		}
		
		diag.resize(b);
		SFKB.solve(diag,MM);
		
		for (uint32_t i = 0; i < diag.size(); ++i) {
			R.normalizeIn(diag[i]);
		}
		
		return diag.size();
	}
};

}


#endif //__LINBOX_invariant_factors_H
