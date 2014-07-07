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

#ifndef __LINBOX_coppersmith_invariant_factors_H
#define __LINBOX_coppersmith_invariant_factors_H

#include "linbox/algorithms/block-coppersmith-domain.h"
#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/matrix/random-matrix.h"

#include <givaro/givpoly1.h>
#include <linbox/ring/givaro-poly.h>
#include <linbox/algorithms/smith-form-kannan-bachem.h>

namespace LinBox
{

using namespace LinBox;

template<class Field_,class Blackbox_>
class CoppersmithInvariantFactors {
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

	CoppersmithInvariantFactors():
		M_(NULL), n_(0), b_(0) {}

	void init(Field& F, const Blackbox& M, size_t b) {
		F_=F;
		M_=M;
		n_=M.rowdim();
		MD_=Domain(F);
		U_=Block(F,b,M.rowdim());
		V_=Block(F,M.rowdim(),b);
		b_=b;

		RandIter RI(F_);
		RandomDenseMatrix<RandIter,Field> RDM(F_,RI);
		RDM.randomFullRank(U_);
		RDM.randomFullRank(V_);
	}

	CoppersmithInvariantFactors(Field& F, const Blackbox& M, size_t b):
		MD_(F), F_(F), M_(&M), n_(M.rowdim()), b_(b),
		U_(F,b,M.rowdim()), V_(F,M.rowdim(),b)
	{
		RandIter RI(F_);
		RandomDenseMatrix<RandIter,Field> RDM(F_,RI);
		RDM.randomFullRank(U_);
		RDM.randomFullRank(V_);
	}

	template <class Mat1, class Mat2>
	CoppersmithInvariantFactors(Field& F, const Blackbox& M, size_t b,
	                            const Mat1& U, const Mat2& V):
		MD_(F), F_(F), M_(&M), n_(M.rowdim()), b_(b),
		U_(F,b,M.rowdim()), V_(F,M.rowdim(),b)
	{
		MD_.copy(U_,U);
		MD_.copy(V_,V);
	}

	template <class PolyRingVector>
	size_t computeFactors(PolyRingVector& diag, int earlyTerm=10)
	{
		typedef BlackboxBlockContainer<Field,Blackbox> BBC;
		typedef BlockCoppersmithDomain<Domain,BBC> BCD;
		BBC blockSeq(M_,F_,U_,V_);
		BCD coppersmith(MD_,&blockSeq,earlyTerm);

		std::vector<size_t> deg;
		std::vector<Block> gen;
		deg=coppersmith.right_minpoly(gen);
		commentator().report(Commentator::LEVEL_IMPORTANT,PROGRESS_REPORT)
			<<"Finished computing minpoly"<<std::endl;
		PolyDom PD(F_,"x");
		PolyRing R(PD);
		PolyMatDom PMD(R);
		size_t d=gen.size();
		PolyBlock MM(R,b_,b_);
		PolyElement temp;
		PD.init(temp,d-1);
		for (int i=0;i<b_;++i) {
			for (int j=0;j<b_;++j) {
				for (int k=0;k<d;++k) {
					PD.setEntry(temp,gen[k].getEntry(i,j),k);
				}
				MM.setEntry(i,j,temp);
			}
		}
		SmithFormKannanBachemDomain<PolyMatDom> SFKB(PMD);
		diag.resize(b_);
		SFKB.solve(diag,MM);
		for (int i=0;i<diag.size();++i) {
			R.normalize(diag[i],diag[i]);
		}

		return diag.size();
	}

protected:

	Domain MD_;

	Field F_;

	const Blackbox *M_;

	size_t n_,b_;

	Block U_,V_;
};

}


#endif //__LINBOX_coppersmith_invariant_factors_H
