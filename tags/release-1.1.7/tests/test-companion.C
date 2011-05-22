/* Copyright (C) LinBox
 *
 * using generic testBlackbox  -bds
 *
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
#include "linbox/blackbox/companion.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;
	static int iterations = 1; // was 10

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.",        TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",          TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	commentator.start("Companion matrix black box test suite", "companion");

	typedef Modular<uint32> Field;
	typedef vector <Field::Element> Vector;
	typedef vector <Field::Element> Polynomial;
	typedef Companion<Field> Blackbox;

	Field F (q);
	Field::Element d; 
	F.init (d, -1);
	Polynomial p(n+1, d);  
	F.init (d, 1); F.assign(p[n], d);

	Blackbox A (F, p);

	pass = pass && testBlackbox(A);

	Blackbox B (F, n);

	pass = pass && testBlackbox(B);

	commentator.stop("companion matrix black box test suite");

	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
