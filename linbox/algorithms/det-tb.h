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
			template<class Matrix>
			void solve_2by2(Element &det, const Matrix &M) const {
				_F.mul(det, M.getEntry(0, 0), M.getEntry(1, 1));
				
				Element tmp;
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.subin(det, tmp);
			}
			
			template<class Matrix>
			void solve_3by3(Element &det, const Matrix &M) const {
				_F.mul(det, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(det, M.getEntry(2, 2));
				
				Element tmp;
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.addin(det, tmp);
				
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.addin(det, tmp);
				
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.subin(det, tmp);
				
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.subin(det, tmp);
				
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.subin(det, tmp);
			}
			
			template<class Matrix>
			void solve_4by4(Element &det, const Matrix &M) const {
				Element tmp;
				
				// +m[1, 1]*m[2, 2]*m[3, 3]*m[4, 4]
				_F.mul(det, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(det, M.getEntry(2, 2));
				_F.mulin(det, M.getEntry(3, 3));
				
				// -m[1, 1]*m[2, 2]*m[3, 4]*m[4, 3]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.subin(det, tmp);
				
				// -m[1, 1]*m[2, 3]*m[3, 2]*m[4, 4]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.subin(det, tmp);
				
				// +m[1, 1]*m[2, 3]*m[3, 4]*m[4, 2]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.addin(det, tmp);
				
				// +m[1, 1]*m[2, 4]*m[3, 2]*m[4, 3]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.addin(det, tmp);
				
				// -m[1, 1]*m[2, 4]*m[3, 3]*m[4, 2]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.subin(det, tmp);
				
				// -m[1, 2]*m[2, 1]*m[3, 3]*m[4, 4]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.subin(det, tmp);
				
				// +m[1, 2]*m[2, 1]*m[3, 4]*m[4, 3]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.addin(det, tmp);
				
				// +m[1, 2]*m[2, 3]*m[3, 1]*m[4, 4]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.addin(det, tmp);
				
				// -m[1, 2]*m[2, 3]*m[3, 4]*m[4, 1]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.subin(det, tmp);
				
				// -m[1, 2]*m[2, 4]*m[3, 1]*m[4, 3]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.subin(det, tmp);
				
				// +m[1, 2]*m[2, 4]*m[3, 3]*m[4, 1]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.addin(det, tmp);
				
				// +m[1, 3]*m[2, 1]*m[3, 2]*m[4, 4]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.addin(det, tmp);
				
				// -m[1, 3]*m[2, 1]*m[3, 4]*m[4, 2]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.subin(det, tmp);
				
				// -m[1, 3]*m[2, 2]*m[3, 1]*m[4, 4]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.subin(det, tmp);
				
				// +m[1, 3]*m[2, 2]*m[3, 4]*m[4, 1]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.addin(det, tmp);
				
				// +m[1, 3]*m[2, 4]*m[3, 1]*m[4, 2]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.addin(det, tmp);
				
				// -m[1, 3]*m[2, 4]*m[3, 2]*m[4, 1]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.subin(det, tmp);
				
				// -m[1, 4]*m[2, 1]*m[3, 2]*m[4, 3]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.subin(det, tmp);
				
				// +m[1, 4]*m[2, 1]*m[3, 3]*m[4, 2]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.addin(det, tmp);
				
				// +m[1, 4]*m[2, 2]*m[3, 1]*m[4, 3]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.addin(det, tmp);
				
				// -m[1, 4]*m[2, 2]*m[3, 3]*m[4, 1]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.subin(det, tmp);
				
				// -m[1, 4]*m[2, 3]*m[3, 1]*m[4, 2]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.subin(det, tmp);
				
				// +m[1, 4]*m[2, 3]*m[3, 2]*m[4, 1]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.addin(det, tmp);
			}
			
			template<class Matrix>
			void solve_5by5(Element &det, const Matrix &M) const {
				Element tmp;
				
				// +m[1, 1]*m[2, 2]*m[3, 3]*m[4, 4]*m[5, 5]
				_F.mul(det, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(det, M.getEntry(2, 2));
				_F.mulin(det, M.getEntry(3, 3));
				_F.mulin(det, M.getEntry(4, 4));
				
				// -m[1, 1]*m[2, 2]*m[3, 3]*m[4, 5]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// -m[1, 1]*m[2, 2]*m[3, 4]*m[4, 3]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 1]*m[2, 2]*m[3, 4]*m[4, 5]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// +m[1, 1]*m[2, 2]*m[3, 5]*m[4, 3]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 1]*m[2, 2]*m[3, 5]*m[4, 4]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// -m[1, 1]*m[2, 3]*m[3, 2]*m[4, 4]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 1]*m[2, 3]*m[3, 2]*m[4, 5]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// +m[1, 1]*m[2, 3]*m[3, 4]*m[4, 2]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 1]*m[2, 3]*m[3, 4]*m[4, 5]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 1]*m[2, 3]*m[3, 5]*m[4, 2]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 1]*m[2, 3]*m[3, 5]*m[4, 4]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// +m[1, 1]*m[2, 4]*m[3, 2]*m[4, 3]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 1]*m[2, 4]*m[3, 2]*m[4, 5]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// -m[1, 1]*m[2, 4]*m[3, 3]*m[4, 2]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 1]*m[2, 4]*m[3, 3]*m[4, 5]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// +m[1, 1]*m[2, 4]*m[3, 5]*m[4, 2]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// -m[1, 1]*m[2, 4]*m[3, 5]*m[4, 3]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 1]*m[2, 5]*m[3, 2]*m[4, 3]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 1]*m[2, 5]*m[3, 2]*m[4, 4]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// +m[1, 1]*m[2, 5]*m[3, 3]*m[4, 2]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 1]*m[2, 5]*m[3, 3]*m[4, 4]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 1]*m[2, 5]*m[3, 4]*m[4, 2]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// +m[1, 1]*m[2, 5]*m[3, 4]*m[4, 3]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 0), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// -m[1, 2]*m[2, 1]*m[3, 3]*m[4, 4]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 2]*m[2, 1]*m[3, 3]*m[4, 5]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// +m[1, 2]*m[2, 1]*m[3, 4]*m[4, 3]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 2]*m[2, 1]*m[3, 4]*m[4, 5]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// -m[1, 2]*m[2, 1]*m[3, 5]*m[4, 3]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 2]*m[2, 1]*m[3, 5]*m[4, 4]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// +m[1, 2]*m[2, 3]*m[3, 1]*m[4, 4]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 2]*m[2, 3]*m[3, 1]*m[4, 5]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// -m[1, 2]*m[2, 3]*m[3, 4]*m[4, 1]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 2]*m[2, 3]*m[3, 4]*m[4, 5]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 2]*m[2, 3]*m[3, 5]*m[4, 1]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 2]*m[2, 3]*m[3, 5]*m[4, 4]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 2]*m[2, 4]*m[3, 1]*m[4, 3]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 2]*m[2, 4]*m[3, 1]*m[4, 5]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// +m[1, 2]*m[2, 4]*m[3, 3]*m[4, 1]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 2]*m[2, 4]*m[3, 3]*m[4, 5]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 2]*m[2, 4]*m[3, 5]*m[4, 1]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// +m[1, 2]*m[2, 4]*m[3, 5]*m[4, 3]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 2]*m[2, 5]*m[3, 1]*m[4, 3]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 2]*m[2, 5]*m[3, 1]*m[4, 4]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// -m[1, 2]*m[2, 5]*m[3, 3]*m[4, 1]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 2]*m[2, 5]*m[3, 3]*m[4, 4]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 2]*m[2, 5]*m[3, 4]*m[4, 1]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// -m[1, 2]*m[2, 5]*m[3, 4]*m[4, 3]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 1), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// +m[1, 3]*m[2, 1]*m[3, 2]*m[4, 4]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 3]*m[2, 1]*m[3, 2]*m[4, 5]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// -m[1, 3]*m[2, 1]*m[3, 4]*m[4, 2]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 3]*m[2, 1]*m[3, 4]*m[4, 5]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// +m[1, 3]*m[2, 1]*m[3, 5]*m[4, 2]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 3]*m[2, 1]*m[3, 5]*m[4, 4]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 3]*m[2, 2]*m[3, 1]*m[4, 4]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 3]*m[2, 2]*m[3, 1]*m[4, 5]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// +m[1, 3]*m[2, 2]*m[3, 4]*m[4, 1]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 3]*m[2, 2]*m[3, 4]*m[4, 5]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 3]*m[2, 2]*m[3, 5]*m[4, 1]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 3]*m[2, 2]*m[3, 5]*m[4, 4]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 3]*m[2, 4]*m[3, 1]*m[4, 2]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 3]*m[2, 4]*m[3, 1]*m[4, 5]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 3]*m[2, 4]*m[3, 2]*m[4, 1]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 3]*m[2, 4]*m[3, 2]*m[4, 5]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 3]*m[2, 4]*m[3, 5]*m[4, 1]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// -m[1, 3]*m[2, 4]*m[3, 5]*m[4, 2]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 3]*m[2, 5]*m[3, 1]*m[4, 2]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 3]*m[2, 5]*m[3, 1]*m[4, 4]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// +m[1, 3]*m[2, 5]*m[3, 2]*m[4, 1]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 3]*m[2, 5]*m[3, 2]*m[4, 4]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 3]*m[2, 5]*m[3, 4]*m[4, 1]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// +m[1, 3]*m[2, 5]*m[3, 4]*m[4, 2]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 2), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// -m[1, 4]*m[2, 1]*m[3, 2]*m[4, 3]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 4]*m[2, 1]*m[3, 2]*m[4, 5]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// +m[1, 4]*m[2, 1]*m[3, 3]*m[4, 2]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 4]*m[2, 1]*m[3, 3]*m[4, 5]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 4]*m[2, 1]*m[3, 5]*m[4, 2]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// +m[1, 4]*m[2, 1]*m[3, 5]*m[4, 3]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// +m[1, 4]*m[2, 2]*m[3, 1]*m[4, 3]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 4]*m[2, 2]*m[3, 1]*m[4, 5]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// -m[1, 4]*m[2, 2]*m[3, 3]*m[4, 1]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 4]*m[2, 2]*m[3, 3]*m[4, 5]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 4]*m[2, 2]*m[3, 5]*m[4, 1]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// -m[1, 4]*m[2, 2]*m[3, 5]*m[4, 3]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 4]*m[2, 3]*m[3, 1]*m[4, 2]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.subin(det, tmp);
				
				// +m[1, 4]*m[2, 3]*m[3, 1]*m[4, 5]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// +m[1, 4]*m[2, 3]*m[3, 2]*m[4, 1]*m[5, 5]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 4));
				_F.addin(det, tmp);
				
				// -m[1, 4]*m[2, 3]*m[3, 2]*m[4, 5]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 4));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 4]*m[2, 3]*m[3, 5]*m[4, 1]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// +m[1, 4]*m[2, 3]*m[3, 5]*m[4, 2]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 4));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 4]*m[2, 5]*m[3, 1]*m[4, 2]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// -m[1, 4]*m[2, 5]*m[3, 1]*m[4, 3]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 4]*m[2, 5]*m[3, 2]*m[4, 1]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// +m[1, 4]*m[2, 5]*m[3, 2]*m[4, 3]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 4]*m[2, 5]*m[3, 3]*m[4, 1]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// -m[1, 4]*m[2, 5]*m[3, 3]*m[4, 2]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 3), M.getEntry(1, 4));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// +m[1, 5]*m[2, 1]*m[3, 2]*m[4, 3]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 5]*m[2, 1]*m[3, 2]*m[4, 4]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// -m[1, 5]*m[2, 1]*m[3, 3]*m[4, 2]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 5]*m[2, 1]*m[3, 3]*m[4, 4]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// +m[1, 5]*m[2, 1]*m[3, 4]*m[4, 2]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// -m[1, 5]*m[2, 1]*m[3, 4]*m[4, 3]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 0));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 5]*m[2, 2]*m[3, 1]*m[4, 3]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 5]*m[2, 2]*m[3, 1]*m[4, 4]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// +m[1, 5]*m[2, 2]*m[3, 3]*m[4, 1]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 5]*m[2, 2]*m[3, 3]*m[4, 4]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 5]*m[2, 2]*m[3, 4]*m[4, 1]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// +m[1, 5]*m[2, 2]*m[3, 4]*m[4, 3]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 1));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 5]*m[2, 3]*m[3, 1]*m[4, 2]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.addin(det, tmp);
				
				// -m[1, 5]*m[2, 3]*m[3, 1]*m[4, 4]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// -m[1, 5]*m[2, 3]*m[3, 2]*m[4, 1]*m[5, 4]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 3));
				_F.subin(det, tmp);
				
				// +m[1, 5]*m[2, 3]*m[3, 2]*m[4, 4]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 3));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
				
				// +m[1, 5]*m[2, 3]*m[3, 4]*m[4, 1]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// -m[1, 5]*m[2, 3]*m[3, 4]*m[4, 2]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 2));
				_F.mulin(tmp, M.getEntry(2, 3));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 5]*m[2, 4]*m[3, 1]*m[4, 2]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.subin(det, tmp);
				
				// +m[1, 5]*m[2, 4]*m[3, 1]*m[4, 3]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 0));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.addin(det, tmp);
				
				// +m[1, 5]*m[2, 4]*m[3, 2]*m[4, 1]*m[5, 3]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 2));
				_F.addin(det, tmp);
				
				// -m[1, 5]*m[2, 4]*m[3, 2]*m[4, 3]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 1));
				_F.mulin(tmp, M.getEntry(3, 2));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.subin(det, tmp);
				
				// -m[1, 5]*m[2, 4]*m[3, 3]*m[4, 1]*m[5, 2]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 0));
				_F.mulin(tmp, M.getEntry(4, 1));
				_F.subin(det, tmp);
				
				// +m[1, 5]*m[2, 4]*m[3, 3]*m[4, 2]*m[5, 1]
				_F.mul(tmp, M.getEntry(0, 4), M.getEntry(1, 3));
				_F.mulin(tmp, M.getEntry(2, 2));
				_F.mulin(tmp, M.getEntry(3, 1));
				_F.mulin(tmp, M.getEntry(4, 0));
				_F.addin(det, tmp);
			}
			
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
			void solve_general(Element &det, const Matrix &M) const {
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
			
			template<class Matrix>
			void solve(Element &det, const Matrix &M) const {
				if (M.rowdim() != M.coldim()) {
					return;
				}
				
				if (M.rowdim() == 2) {
					solve_2by2(det, M);
				} else if (M.rowdim() == 3) {
					solve_3by3(det, M);
				} else if (M.rowdim() == 4) {
					solve_4by4(det, M);
				} else if (M.rowdim() == 5) {
					solve_5by5(det, M);
				} else {
					solve_general(det, M);
				}
			}
	};
}

#endif // __LINBOX_det_textbook_domain_H