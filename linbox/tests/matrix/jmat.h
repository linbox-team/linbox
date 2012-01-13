/* Copyright (c) LinBox
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



namespace LinBox{

size_t lambda = 1000000000;

template<class Ring, class Matrix>
Matrix& jordanform(const Ring& R, Matrix& Mat, size_t n){
	size_t val = 1;
	typename Ring::Element one, lam;
	R.init(one, val);
	R.init(lam, lambda);

	Mat.setEntry(0, 0, one);

	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			if(i == j) Mat.setEntry(i, j, lam);
			if(j+1 == i) Mat.setEntry(i, j, one);
		}
	}

	return Mat;
}

} //LinBox
