
/** 
 * examples/modular-int.C
 *
 * Copyright (C) 2005, 2010 D. Saunders, Z. Wang
 *
 * This file is part of LinBox.
 *
 *   LinBox is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU Lesser General Public License as
 *   published by the Free Software Foundation, either version 2 of
 *   the License, or (at your option) any later version.
 *
 *   LinBox is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU Lesser General Public License for more details.
 *
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with LinBox.  If not, see 
 *   <http://www.gnu.org/licenses/>.
 */

/** \file examples/fields/modular-int.C
\brief  Example of arithmetic in the Modular<int> finite field.
 */
//by  Zhendong wan wan@udel.edu


#include <iostream>

/* the header file for modular<int>*/
#include "linbox/field/modular-int.h"

int main (int argc, char **argv)
{
	/* construct the Z/101Z field */
	LinBox::Modular<int> F(101);

	/* declare local variable a, b, c, x, y*/
	LinBox::Modular<int>::Element a, b, c, x, y;

	// initializtion
	F.init(a, 92);

	F.init(b, 50);

	F.init(c, 56);

	F.init(x, 33);

	F.init(y, 45);

	/* +, -, *, / arithmetic operation */
	//compute c <- a + b  mod 101;
	F.add(c,a,b);
        //output value of a
        F.write(std::cout,a);
        std::cout << " + ";
        //output value of b
        F.write(std::cout,b);
        std::cout << " = ";
        //output value of c
        F.write(std::cout,c);
        std::cout << " mod 101" << "\n";

	//compute c <- a - b mod 101;
	F.sub(c,a,b);
	//output value of a
	F.write(std::cout,a);
	std::cout << " - ";
	//output value of b
	F.write(std::cout,b);
	std::cout << " = ";
	//output value of c
	F.write(std::cout,c);
	std::cout << " mod 101" << "\n";

	//compute c <- a * b mod 101;
	F.mul(c,a,b);
        //output value of a
        F.write(std::cout,a);
        std::cout << " * ";
        //output value of b
        F.write(std::cout,b);
        std::cout << " = ";
        //output value of c
        F.write(std::cout,c);
        std::cout << " mod 101" << "\n";

	//compute c <- a / b mod 101;
	F.div(c,a,b);
        //output value of a
        F.write(std::cout,a);
        std::cout << " / ";
        //output value of b
        F.write(std::cout,b);
        std::cout << " = ";
        //output value of c
        F.write(std::cout,c);
        std::cout << " mod 101" << "\n";

	//compute c = 1 /a mod 101;
	F.inv(c,a);
	std::cout << "1 / " << a << " = ";
	//output value of c
        F.write(std::cout,c);
        std::cout << " mod 101" << "\n";

	/* Inplace operation */ 
	//compute a *= b;
	F.mulin(a,b);
	
	//compute a /= b;
	F.divin(a,b);

	//compute a <- 1 / a/;
	F.invin(a);

	/* a * x + y operation */
	//compute c <- a * x + y;
	F.axpy(c,a,x,y);

	//compute c += a * x;
	F.axpyin(c,a,x);

	/* compare operation */
	// test if a == 1 mod 101; 
	//output value of a
	F.write(std::cout,a);
	std::cout << " == 1 mod 101 " << (F.isOne(a) ? "true" : "false") << "\n";

	//test if a == 0 mod 101;
        //output value of a
        F.write(std::cout,a);
	std::cout << " == 0 mod 101 " << (F.isZero(a) ? "true" : "false") << "\n" ;

	// test if a == b mod 101;
	//output value of a
        F.write(std::cout,a);
	std::cout << " == ";
	//output value of b 
        F.write(std::cout,b);
	std::cout << " mod 101 " << (F.areEqual(a,b) ? "true" : "false") << "\n";
	
	//For performance, we recommend to use all possible operation, for example, using axpy, instead of mul, then addin.

	return 0;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
