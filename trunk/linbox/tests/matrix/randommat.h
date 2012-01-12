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
 * License along with this library; if not, write to the Free Software Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========

 */

#include "coin.h"
#include<cstdlib>
#include<vector>

//#include<iostream>
//using namespace std;

namespace LinBox{

bool used(vector<int> &array, int value){
	vector<int>::iterator i;
	for(i = array.begin(); i != array.end(); ++i)
		if( (*i) == value ) return true;  //  already used location
	return false;  //  used
}

template<class Ring, class Matrix>
Matrix& randomMat(const Ring& R, Matrix& Mat, size_t n, size_t epr){

	epr = epr > n ? n : epr;

	int val, pos, neg;
	vector<int> usedV(epr);
	typename Ring::Element tmp;

	//srand(time(NULL));

	for(size_t i = 0; i < n; ++i) {
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
			//Mat.setEntry(i, pos, R.init(tmp, 1));
			Mat.setEntry(i, pos, R.init(tmp, val));
		}
	}

	return Mat;
}

}
