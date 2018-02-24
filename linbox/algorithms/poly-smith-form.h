
/* tests/__LINBOX_smith_form_kannan_bachem.h
 * Copyright (C) 2017 Gavin Harrison,
 *
 * Written by Gavin Harrison <gmh33@drexel.edu>,
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

#include <iostream>
#include <stdlib.h>
#include <algorithm>

#include "linbox/matrix/matrixdomain/matrix-domain.h"

#include "linbox/algorithms/poly-dixon.h"

#ifndef __LINBOX_poly_smith_form_domain_H
#define __LINBOX_poly_smith_form_domain_H

namespace LinBox
{
	template<class Ring>
	class PolySmithFormDomain
	{
	public:
		typedef typename Ring::Element Polynomial;
		typedef typename Ring::Coeff Coeff;
		typedef typename Ring::RandIter RandIter;
		
		typedef MatrixDomain<Ring> MatrixDom;
		
		typedef typename MatrixDom::OwnMatrix Matrix;
		typedef typename MatrixDom::Matrix SubMatrix;

	private:
		Ring _R;
		RandIter _RI;
		MatrixDom _MD;

	public:
		PolySmithFormDomain(const Ring &R) : _R(R), _RI(R), _MD(_R) {}
		PolySmithFormDomain(const PolySmithFormDomain &D) : _R(D._R), _RI(D._RI), _MD(D._MD) {}
		
		template<class Matrix1>
		bool dixon(
			Polynomial &minpoly,
			const Matrix1 &M,
			const Polynomial &f,
			size_t max_deg) const {
		
			typedef typename Ring::QuotientRing QuotientRing;
			QuotientRing QR(_R.quotient(f));
			PolyDixonDomain<Ring, QuotientRing> DD(_R, QR);
			
			Matrix Mx(M);
					
			Matrix b(_R, Mx.rowdim(), 1);
			for (size_t i = 0; i < b.rowdim(); i++) {
				Polynomial e;
				_RI.random(e, max_deg);
				
				b.setEntry(i, 0, e);
			}
			
			Matrix x(_R, Mx.rowdim(), 1);
			size_t m = 4 * ((max_deg / _R.deg(f)) + 1);
						
			if (!DD.solve(x, Mx, b, f, m)) {
				return false;
			}
			
			Polynomial fm;
			_R.pow(fm, f, m);
			
			_R.assign(minpoly, _R.one);
			for (size_t i = 0; i < x.rowdim(); i++) {
				Polynomial numer, denom, tmp;
				DD.rat_recon(numer, denom, x.getEntry(i, 0), fm);
				
				_R.lcm(tmp, minpoly, denom);
				_R.assign(minpoly, tmp);
			}
			
			return true;
		}
		
		template<class Matrix1>
		void dixon(Polynomial &minpoly, const Matrix1 &M, size_t m) const {
			for (size_t d = 1; d < 11; d++) {
				Polynomial f;
				_RI.randomIrreducible(f, d++);
				if (dixon(minpoly, M, f, m)) break;
			}
		}
	}; // end of class PolySmithFormDomain
}

#endif // __LINBOX_poly_smith_form_domain_H
