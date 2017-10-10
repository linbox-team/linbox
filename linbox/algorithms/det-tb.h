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

#include <vector>
#include <list>

#ifndef __LINBOX_det_textbook_domain_H
#define __LINBOX_det_textbook_domain_H

namespace LinBox {
	template<class Field>
	class DetTextbookDomain {
		typedef typename Field::Element Element;
		
		typedef MatrixDomain<Field> MatrixDom;
		
		typedef typename MatrixDom::OwnMatrix OwnMatrix;
		typedef typename MatrixDom::Matrix SubMatrix;
	
	private:
		Field _F;
		MatrixDom _MD;
	
	public:
		DetTextbookDomain(const Field &F) : _F(F), _MD(_F) {}
		DetTextbookDomain(const DetTextbookDomain &D) : _F(D._F), _MD(D._MD) {}
		
	private:
		int sgn(const std::vector<size_t> &p) const {
			std::vector<bool> visited(p.size(), false);
			int sgn = 1;
			
			for (size_t k = 0; k < p.size(); k++) {
				if (!visited[k]) {
					size_t next = k;
					size_t len = 0;
					
					while (!visited[next]) {
						len++;
						visited[next] = true;
						next = p[next];
					}
					
					if (len % 2 == 0) {
						sgn *= -1;
					}
				}
			}
			
			return sgn;
		}
		
		void perm(std::vector<size_t> &p, size_t len, size_t k) const {
			std::list<size_t> values;
			for (size_t i = 0; i < len; i++) {
				values.push_back(i);
			}
			
			p.resize(0);
			for (size_t i = 0; i < len - 1; i++) {
				size_t j = k % values.size();
				k = k / values.size();
				
				auto it = values.begin();
				std::advance(it, j);
				
				p.push_back(*it);
				values.erase(it);
			}
			
			p.push_back(*values.begin());
		}
		
		size_t factorial(size_t n) const {
			size_t rv = 1;
			for (size_t i = n; i > 1; i--) {
				rv *= i;
			}
			return rv;
		}
		
	public:
		template<class Matrix>
		void solve(Element &det, const Matrix &M) const {
			if (M.rowdim() != M.coldim()) {
				return;
			}
			
			size_t max = factorial(M.rowdim());
			
			_F.assign(det, _F.zero);
			for (size_t i = 0; i < max; i++) {
				std::vector<size_t> p;
				perm(p, M.rowdim(), i);
				
				Element tmp;
				_F.assign(tmp, _F.one);
				for (size_t j = 0; j < M.rowdim(); j++) {
					_F.mulin(tmp, M.getEntry(j, p[j]));
				}
				
				if (sgn(p) == -1) {
					_F.negin(tmp);
				}
				
				_F.addin(det, tmp);
			}
		}
	};
}

#endif // __LINBOX_det_textbook_domain_H