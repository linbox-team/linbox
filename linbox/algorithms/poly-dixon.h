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

#include "linbox/ring/givaro-poly.h"

#ifndef __LINBOX_poly_dixon_domain_H
#define __LINBOX_poly_dixon_domain_H

namespace LinBox {
	template<class Field>
	class PolyDixonDomain {
		typedef typename Field::Element Element;
		
		typedef MatrixDomain<Field> MatrixDom;
		
		typedef typename MatrixDom::OwnMatrix OwnMatrix;
		typedef typename MatrixDom::Matrix SubMatrix;
		
		private:
			Field _F;
			MatrixDom _MD;
		
		public:
			PolyDixonDomain(const Field &F) : _F(F), _MD(_F) {}
			PolyDixonDomain(const PolyDixonDomain &D) : _F(D._F), _MD(D._MD) {}
			
		private:
			// Methods
			
		public:
			void rat_recon(Element &numer, Element &denom, const Element &s, const Element &h) {
				int64_t N = _F.degree(h).value() / 2;
				
				Element v0, v1;
				_F.assign(v0, h);
				_F.assign(v1, _F.zero);
				
				Element w0, w1;
				_F.assign(w0, s);
				_F.assign(w1, _F.one);
				
				while (_F.degree(w0).value() > N) {
					Element q, z0, z1;
					
					_F.quo(q, v0, w0);
					_F.negin(q);
					
					_F.axpy(z0, q, w0, v0);
					_F.axpy(z1, q, w1, v1);
					
					_F.assign(v0, w0);
					_F.assign(w0, z0);
					
					_F.assign(v1, w1);
					_F.assign(w1, z1);
				}
				
				_F.assign(numer, w0);
				_F.assign(denom, w1);
			}
			
			// Solve x such that Ax = y
			template<class Matrix1, class Matrix2, class Matrix3>
			void solve(Matrix1 &x, const Matrix2 &A, const Matrix3 &y) {
				
			}
	};
}

#endif // __LINBOX_poly_dixon_domain_H