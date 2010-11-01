/* Copyright (C) LinBox
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
#include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/transpose.h"

#include "test-common.h"
#include "test-generic.h"


template<class Vec>
bool eqVec(const Vec& a, const Vec& b)
{
	bool good = true;
	for (typename Vec::size_type i = 0; i < a.size(); ++i) good = good && (a[i] == b[i]);
	return good;
}

using namespace LinBox;

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t n = 200;
	static integer q = 2147483647U;
	static int iterations = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	srand (time (NULL));

	commentator.start("triplesbb black box test suite", "triplesbb");

	typedef Modular<uint32> Field;
	typedef Field::Element Element;
	typedef vector <Element> Vector;
	typedef TriplesBB<Field> Blackbox;

	Field F (q);
	Element d; 
	F.init (d, -1);

	// set up the matrix
	std::vector<Element> values; 
	std::vector<size_t> rowP; 
	std::vector<size_t> colP; 
	for(int i = 1; i < 12; ++i)
	{
		values.push_back(i);
		rowP.push_back(i%n + 1);
		colP.push_back(i/n + 1);
	}

	Blackbox A(F, values, rowP, colP, n, n);

	pass = pass && testBlackbox(A);

	Transpose<Blackbox> B(&A);

	pass = pass && testBlackbox(B);

	Vector x(n), y(n), z(n);
	for(int i = 0; i < 5; ++i) x[i] = i;

	A.apply(y, x);
	B.applyTranspose(z, x);
	pass = pass && eqVec(y, z);

	B.apply(y, x);
	A.applyTranspose(z, x);
	pass = pass && eqVec(y, z);

	Blackbox C(F, n, n);
	for(size_t i = 0; i < rowP.size(); ++i) C.addEntry(values[i], rowP[i], colP[i]);
	pass = pass && testBlackbox(C);

	commentator.stop("triplesbb black box test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
