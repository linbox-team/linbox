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

#include "linbox/algorithms/invert-tb.h"

#ifndef __LINBOX_poly_dixon_domain_H
#define __LINBOX_poly_dixon_domain_H

namespace LinBox {
	template<class Ring, class Field>
	class PolyDixonDomain {
		typedef typename Ring::Element Element;
		typedef typename Field::Element FElement;
		
		typedef MatrixDomain<Ring> MatrixDom;
		typedef typename MatrixDom::OwnMatrix Matrix;
		typedef typename MatrixDom::Matrix SubMatrix;
		
		typedef MatrixDomain<Field> FMatrixDom;
		typedef typename FMatrixDom::OwnMatrix FMatrix;
		typedef typename FMatrixDom::Matrix FSubMatrix;
		
		private:
			Ring _R;
			Field _F;
			MatrixDom _MD;
			FMatrixDom _FMD;
			InvertTextbookDomain<Field> _ID;
			
		public:
			PolyDixonDomain(const Ring &R, const Field &F) : _R(R), _F(F), _MD(_R), _FMD(_F), _ID(_F) {}
			PolyDixonDomain(const PolyDixonDomain &D) : _R(D._R), _F(D._F), _MD(_R), _FMD(D._MD), _ID(_F) {}
			
		private:
			// Methods
			
		public:
			void rat_recon(Element &numer, Element &denom, const Element &s, const Element &h) {
				size_t N = _R.deg(h) / 2;
				
				Element v0, v1;
				_R.assign(v0, h);
				_R.assign(v1, _R.zero);
				
				Element w0, w1;
				_R.assign(w0, s);
				_R.assign(w1, _R.one);
				
				while (_R.deg(w0) > N) {
					Element q, z0, z1;
					
					_R.quo(q, v0, w0);
					_R.negin(q);
					
					_R.axpy(z0, q, w0, v0);
					_R.axpy(z1, q, w1, v1);
					
					_R.assign(v0, w0);
					_R.assign(w0, z0);
					
					_R.assign(v1, w1);
					_R.assign(w1, z1);
				}
				
				_R.assign(numer, w0);
				_R.assign(denom, w1);
			}
			
			void divin(Matrix &b, const Element &divisor) {
				for (size_t i = 0; i < b.rowdim(); i++) {
					for (size_t j = 0; j < b.coldim(); j++) {
						Element z;
						_R.div(z, b.getEntry(i, j), divisor);
						b.setEntry(i, j, z);
					}
				}
			}
			
			// Solve x such that Ax = y (mod f^m)
			bool solve(Matrix &x, const Matrix &A, const Matrix &b,
				const Element &f, size_t m) 
			{
				size_t r = A.rowdim();
				
				FMatrix Af(A, _F);
				FMatrix Ai(_F, A.rowdim(), A.coldim());
				
				bool success = _ID.invert(Ai, Af);
				if (!success) {
					return false;
				}
				
				x.zero();
				
				Matrix bi(b);
				FMatrix xi(_F, r, 1);
				Element fi;
				_R.assign(fi, _R.one);
				for (size_t i = 0; i < m; i++) {
					FMatrix bif(bi, _F);
					_FMD.mul(xi, Ai, bif);
					Matrix xir(xi, _R);
					
					Matrix Axir(_R, r, 1);
					_MD.mul(Axir, A, xir);
					
					Matrix fixir(_R, r, 1);
					_MD.mul(fixir, xir, fi);
					_R.mulin(fi, f);
					_MD.addin(x, fixir);
					
					_MD.subin(bi, Axir);
					divin(bi, f);
				}
				
				return true;
			}
	};
}

#endif // __LINBOX_poly_dixon_domain_H
