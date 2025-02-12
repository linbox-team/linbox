/* linbox/algorithms/cia.h
 * Copyright(C) LinBox
 *
 *  Written by Clement Pernet <clement.pernet@imag.fr>
 *
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
 *.
 */

#ifndef __LINBOX_cia_H
#define __LINBOX_cia_H

#include <givaro/givpoly1factor.h>
#include "linbox/ring/modular.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/solutions/minpoly.h"

namespace LinBox
{

	/*! @ingroup algorithms
	 * Algorithm computing the integer characteristic polynomial
	 * of a dense matrix.
	 *
	 * @bib [Dumas-Pernet-Wan ISSAC05]
	 *
	 *
	 */
	template < class Polynomial, class Blackbox >
	Polynomial& cia (Polynomial & P, const Blackbox & A,
			 const Method::DenseElimination  & M)
	{
		commentator().start ("Integer Givaro::Dense Charpoly ", "CIA");

		using Ring = typename Blackbox::Field;
		Ring ZZ = A.field();
		typedef Givaro::Modular<double> Field;
		typedef typename Blackbox::template rebind<Field>::other FBlackbox;
		typedef PolynomialRing<Ring> IntPolyDom;
		typedef PolynomialRing<Field> FieldPolyDom;
		typedef typename IntPolyDom::Element IntPoly;
		typedef typename FieldPolyDom::Element FieldPoly;

		IntPolyDom IPD(ZZ);

		/* Computation of the integer minimal polynomial */
		IntPoly intMinPoly(ZZ);
		minpoly (intMinPoly, A, RingCategories::IntegerTag(), M);

        if (intMinPoly.size() == A.coldim()+1){
            commentator().stop ("done", NULL, "CIA");
            return P = intMinPoly;
        }
//         IPD.write(std::cerr<<"Minpoly = ", intMinPoly) << std::endl;

		/* Factorization over the integers */
		std::vector<IntPoly> intFactors;
		std::vector<uint64_t> mult;
		IPD.factor (intFactors, mult, intMinPoly);
		size_t nf = intFactors.size();

		/* One modular characteristic polynomial computation */
		PrimeIterator<IteratorCategories::HeuristicTag> primeg (FieldTraits<Field>::bestBitSize(A.coldim()));
		++primeg;
		Field F(*primeg);
		FBlackbox fbb(F, A.rowdim(), A.coldim());
		FieldPoly fieldCharPoly(F);
		MatrixHom::map(fbb, A);
		charpoly (fieldCharPoly, fbb, M);
		/* Determination of the multiplicities */
		FieldPolyDom FPD (F);
		std::vector<FieldPoly> fieldFactors (nf);
		integer tmp_convert; // PG 2005-08-04
		for (size_t i = 0; i < nf; ++i){
			size_t d= intFactors[i].size();
			fieldFactors[i].resize(d);
			for (size_t j = 0; j < d; ++j)
				//F.init ((fieldFactors[i])[j], (*intFactors[i])[j]);
				F.init ((fieldFactors[i])[j], ZZ.convert(tmp_convert,(intFactors[i])[j]));// PG 2005-08-04
		}

		FieldPoly currPol = fieldCharPoly;
		FieldPoly r,tmp,q;
		std::vector<long> multip (nf);
		for (size_t i = 0; i < nf; ++i) {
			FieldPoly currFact = fieldFactors[i];
			r.clear();
			int m=0;
			q=currPol;
			do{
				currPol = q;
				FPD.divmod (q, r, currPol, currFact);
				m++;
			} while (FPD.isZero (r));
			multip[i] = m-1;
		}

		IntPoly intCharPoly (ZZ);
		IPD.assign (intCharPoly, IPD.one);
		for (size_t i = 0; i < nf; ++i){
			IPD.pow( P, intFactors[i], multip[i] );
			IPD.mulin( intCharPoly, P );
		}
		//for (size_t i = 0; i < nf; ++i)
		//delete intFactors[i];
		commentator().stop ("done", NULL, "CIA");

		return P = intCharPoly;
	}
}

#endif // __LINBOX_cia_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
