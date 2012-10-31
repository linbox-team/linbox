
/*
 * Copyright (c) LinBox
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

template<class Ring, class Matrix>
Matrix& minmat(const Ring& R, Matrix& Mat, size_t n){
	size_t val;
	typename Ring::Element tmp;

	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			val = i < j ? i : j;
			R.init(tmp, val+1);
			Mat.setEntry(i, j, tmp);
		}
	}

	return Mat;
}

template<class Ring, class Matrix>
Matrix& maxmat(const Ring& R, Matrix& Mat, size_t n){
	size_t val;
	typename Ring::Element tmp;

	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			val = i > j ? i : j;
			R.init(tmp, val+1);
			Mat.setEntry(i, j, tmp);
		}
	}

	return Mat;
}

template<class Ring, class Matrix>
Matrix& qlehmer(const Ring& R, Matrix& Mat, size_t n){
	size_t val;
	typename Ring::Element tmp, tmp2;

	for(size_t i = 0; i < n; ++i){
		for(size_t j = 0; j < n; ++j){
			val = i > j ? i : j;
			val++;
			R.init(tmp, val);
			R.init(tmp2, val);
			R.mulin(tmp, tmp2);
			Mat.setEntry(i, j, tmp);
		}
	}

	return Mat;
}

} //LinBox

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

