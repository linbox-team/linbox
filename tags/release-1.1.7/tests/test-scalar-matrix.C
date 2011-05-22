/* Copyright (C) LinBox
 *
 *
 * using generic testBlackbox  -bds
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/scalar-matrix.h"

#include "test-blackbox.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 20;
	static integer q = 2147483647U;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	commentator.start("Scalar black box test suite", "Scalar");

	typedef Modular<uint32> Field;

	Field F (q);
	Field::Element d; 
	F.init (d, -1);

	typedef ScalarMatrix <Field> Blackbox;
	// C Pernet: why is the default constructor even public?
	// Blackbox A; // Test the default constructor
	//pass = pass && testBlackbox(A);

	Blackbox B (F, n, d); // Test a small one.
	pass = pass && testBlackbox(B);

	//Blackbox C (F, 100000, d); // Test a large one.
	//pass = pass && testBlackbox(C);

	commentator.stop(MSG_STATUS(pass));
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
