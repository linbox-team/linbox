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

void isPower2(size_t n){
	size_t p = 1;
	while(p < n) p <<= 1;
	if(p != n){} //  THROW ERROR };
}

template<class Ring, class Matrix>
Matrix& hadamard(const Ring& R, Matrix& Mat, size_t n){
	size_t val = 1;
	typename Ring::Element tmp;

	Mat.setEntry(0, 0, R.init(tmp, val));

	for(size_t k = 1; k < n; k*=2){
		for(size_t i = 0; i < k; ++i){
			for(size_t j = 0; j < k; ++j){
				Mat.getEntry(tmp, i, j);
				Mat.setEntry(i+k, j, tmp);
				Mat.setEntry(i, j+k, tmp);
				Mat.setEntry(i+k, j+k, R.negin(tmp));
			}
		}
	}

	return Mat;
}

}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

