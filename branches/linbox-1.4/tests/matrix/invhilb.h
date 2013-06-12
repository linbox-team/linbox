
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

#include<cstdlib>
#include "linbox/integer.h"

//#include<iostream>
//using namespace std;
namespace LinBox{

template<class Ring, class Matrix>
Matrix& invhilb(const Ring& R, Matrix& Mat, int n){
	typedef typename Ring::Element Element;

	integer p = n;
	integer r;
	Element val;
	int i, j, i_ind, j_ind;
	for(i_ind=0, i = 1; i <= n; i++, i_ind++){
		if(i>1) p = ((n-i+1)*p*(n+i-1))/((i-1)*(i-1));
		r = p*p;
		R.init(val, (r/(2*i-1)));
		Mat.setEntry((size_t)i_ind,(size_t) i_ind, val);
		for(j_ind=i, j = i+1; j <= n; j++, j_ind++){
			r = -1*((n-j+1)*r*(n+j-1))/((j-1)*(j-1));
			R.init(val, (r/(i+j-1)));
			Mat.setEntry((size_t)i_ind, (size_t)j_ind, val);
			Mat.setEntry((size_t)j_ind, (size_t)i_ind, val);
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

