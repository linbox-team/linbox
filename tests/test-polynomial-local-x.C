/* tests/test-polynomial-local-x.C
 * Copyright (C) 2017 The LinBox group,
 *
 * Written by Gavin Harrison <gmh33@drexel.edu>
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

#include "linbox/linbox-config.h"

#include <iostream>

#include <linbox/ring/ntl.h>
#include <linbox/ring/polynomial-local-x.h>
#include "linbox/ring/modular.h"
#include <givaro/givquotientdomain.h>

#include "test-field.h"

using namespace LinBox;

typedef NTL_zz_pX Field;
typedef Field::Coeff Coeff;
typedef Field::CoeffField CoeffField;
typedef Field::Element Polynomial;
typedef PolynomialLocalX<Field> Ring;
typedef Ring::Element Element;

int main (int argc, char **argv)
{
	int p = 3;
	// int e = 1; // currently not used
	
	CoeffField CF(p);
	Field F(p);
	Ring R(F, 5);
	
	R.write(std::cout) << std::endl;
	R.write(std::cout, R.one) << std::endl;
	
	std::vector<integer> v;
	
	Element x, a, b;
	
	R.init(x, v = {0, 1});
	R.init(a, v = {0, 1, 3});
	R.init(b, v = {0, 0, 10, 7});
	
	R.write(std::cout << "a = ", a) << std::endl;
	R.writeNormalized(std::cout << "a = ", a) << std::endl;
	R.write(std::cout << "b = ", b) << std::endl;
	R.writeNormalized(std::cout << "b = ", b) << std::endl;
	
	Element q;
	R.div(q, b, a);
	R.write(std::cout << "q = ", q) << std::endl;
	R.writeNormalized(std::cout << "q = ", q) << std::endl;
	
	Element qa;
	R.mul(qa, q, a);
	R.write(std::cout << "q * a == ", qa) << std::endl;
	R.writeNormalized(std::cout << "q * a == ", qa) << std::endl;
	
	R.div(q, a, x);
	R.write(std::cout << "a / x == ", q) << std::endl;
	R.writeNormalized(std::cout << "a / x == ", q) << std::endl;
	
	R.mul(q, a, b);
	R.write(std::cout << "a * b == ", q) << std::endl;
	R.writeNormalized(std::cout << "a * b == ", q) << std::endl;
	
	R.add(q, a, b);
	R.write(std::cout << "a + b == ", q) << std::endl;
	R.writeNormalized(std::cout << "a + b == ", q) << std::endl;
	
	Element x1, y;
	R.init(a, v = {1}); // 1
	R.init(x1, v = {1, 1, 1}); // 1+x+x^2
	R.init(y, v = {2, 0, 1}); // 2+x^2
	
	Element z;
	R.axpy(z, a, x1, y);
	R.write(std::cout << "a * x + y == ", z) << std::endl;
	R.writeNormalized(std::cout << "a * x + y == ", z) << std::endl;
	
	R.add(z, x1, y);
	R.write(std::cout << "x + y == ", z) << std::endl;
	R.writeNormalized(std::cout << "x + y == ", z) << std::endl;
	
	return 0;
}