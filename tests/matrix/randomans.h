//  randomans.h  - Random Almost NonSingular matrix

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
#include<vector>
#include "coin.h"

//#include<iostream>
//using namespace std;

namespace LinBox{

template<class Ring, class Matrix>
Matrix& randomAns(const Ring& R, Matrix& Mat, size_t n, size_t epr){
	epr = epr > n ? n : epr;

	int val, pos, neg;
	vector<int> usedV(epr);
	typename Ring::Element tmp;

	srand((unsigned)time(NULL));

	//  build first n-1 rows
	for(size_t i = 0; i < n-1; ++i) {
		//std::cerr << "in row " << i << endl;
		//  reset used vector
		for(size_t q = 0; q < epr; ++q)
			usedV[q] = -1;

		//  epr elements per row
		for(size_t k = 0; k < epr; ++k){
			neg = 1;
			if(rand()%2) neg = -1;
			//  generate value in [0, ceiling)
			val = neg*(rand()%CEILING);
			//  choose random location for value
			do{
				pos = int(rand()%n);  // pos in [0, n)
			}
			while(used(usedV, pos));
			usedV[k] = pos;  //  record location

			//std::cerr << "\t set value " << val << " in pos " << pos << std::endl;
			//  finally, set entry
			Mat.setEntry(i, pos, R.init(tmp, val));
		}
	}

	//  to store each element in our final row
	vector<typename Ring::Element> tmps(n);
	for(size_t i=0; i<n; ++i)
		R.init(tmps[i], 0);

	//  build last row incorporate every row
	for(size_t i = 0; i < n-1; ++i){
		neg = 1;
		if(rand()%2) neg = -1;
		val = neg*(rand()%COMB_CEILING);
		//  val is the multiplier for this row
		for(size_t j = 0; j < n; ++j){
			R.axpyin(tmps[j], val, Mat.getEntry(i, j));
		}
	}
	for(size_t i=0; i<n; ++i)
		Mat.setEntry(n-1, i, tmps[i]);

	//  add one to random value
	size_t wildcard = rand()%n;
	R.init(tmp, Mat.getEntry(n-1, wildcard));
	R.addin(tmp, 1);
	Mat.setEntry(n-1, wildcard, tmp);

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

