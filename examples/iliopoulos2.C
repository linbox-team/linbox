/* linbox/algorithms/coppersmith-invariant-factors.h
 * Copyright (C) 2015 Gavin Harrison
 *
 * Written by Gavin Harrison <gavin.har@gmail.com>
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

#include "linbox/linbox-config.h"

#include <iostream>
#include <vector>

#include "linbox/ring/modular.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/matrix/matrix-domain.h"
#include "linbox/matrix/random-matrix.h"

#include "linbox/algorithms/smith-form-iliopoulos2.h"

#include <givaro/zring.h>

using namespace LinBox;
using namespace Givaro;

typedef ZRing<Integer> Field;
typedef typename ZRing<Integer>::Element Element;
typedef IliopoulosDomain<Field> IliopoulosDom;
typedef MatrixDomain<Field> Domain;
typedef typename Domain::OwnMatrix Matrix;
typedef typename Field::RandIter RandIter;
typedef RandomDenseMatrix<RandIter, Field> RandomMatrix;

Field Z;
Domain MD(Z);
RandIter RI(Z);
RandomMatrix RDM(Z, RI);
IliopoulosDom ID(Z);

void printMatrix(Matrix &A) {
	size_t m = A.rowdim();
	size_t n = A.coldim();
	
	std::cout << "matrix(ZZ, " << m << "," << n << ", [" << std::endl;
	for (size_t i = 0; i < m; i++) {
		Element tmp;
		Z.write(std::cout << "[", A.getEntry(tmp, i, 0));
		for (size_t j = 1; j < n; j++) {
			Z.write(std::cout << ", ", A.getEntry(tmp, i, j));
		}
		std::cout << "]," << std::endl;
	}
	std::cout << "])" << std::endl;
}

void randomL(Matrix &L) {
	size_t n = L.rowdim();
	
	RDM.random(L);
	
	for (size_t i = 0; i < n; i++) {
		L.setEntry(i, i, Z.one);
		
		for (size_t j = i+1; j < n; j++) {
			L.setEntry(i, j, Z.zero);
		}
	}
}

void randomU(Matrix &U) {
	size_t n = U.rowdim();
	
	RDM.random(U);
	
	for (size_t i = 0; i < n; i++) {
		U.setEntry(i, i, Z.one);
		
		for (size_t j = 0; j < i; j++) {
			U.setEntry(i, j, Z.zero);
		}
	}
}

int main(int argc, char** argv)
{
	Matrix A(Z, 4, 4);
	
	//RDM.random(A);
	A.setEntry(0,0,Integer(2));
	A.setEntry(1,1,Integer(6));
	A.setEntry(2,2,Integer(12));
	printMatrix(A);
	
	Matrix L1(Z, 4, 4);
	randomL(L1);
	
	Matrix U1(Z, 4, 4);
	randomU(U1);
	
	Matrix L2(Z, 4, 4);
	randomL(L2);
	
	Matrix U2(Z, 4, 4);
	randomU(U2);
	
	Matrix B(Z, 4, 4);
	
	MD.mul(B, U1, A);
	MD.leftMulin(B, L1);
	MD.rightMulin(B, U2);
	MD.rightMulin(B, L2);
	
	printMatrix(B);
	
	ID.smithFormIn(B, Integer(2 * 6 * 12));
	printMatrix(B);
	
	return 0;
}


