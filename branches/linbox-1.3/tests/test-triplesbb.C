/* Copyright (C) LinBox
 *
 * using generic testBlackbox  -bds
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file  tests/test-triplesbb.C
 * @ingroup tests
 *
 * @brief no doc
 *
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/vector/vector-domain.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	ofstream report;

	bool pass = true;

	static size_t m = 4;
	static size_t n = 20;
	static integer q = 2147483647U;
	q = 101;
	static int iterations = 1;

	static Argument args[] = {
		{ 'm', "-m M", "Set dimension of test matrices to NxN.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	srand ((unsigned)time (NULL));

	commentator().start("TriplesBB black box test suite", "triplesbb");

	//typedef Modular<uint32_t> Field;
	typedef Modular<double> Field;
	typedef Field::Element Element;
	typedef vector <Element> Vector;
	typedef TriplesBB<Field> Blackbox;

	Field F (q);

	// set up the matrix
	Blackbox A(F, m, n);
	Element d;

	for(int i = 0; i < (int)m; ++i)
	{
		F.init(d, i);
		A.setEntry((size_t)i,(size_t) (2*i)%n, d);
		F.init(d, 1);
		A.setEntry((size_t)i, (size_t)(2*i+1)%n, d);
	}

	pass = pass && testReadWrite(A);
	pass = pass && testBlackbox(A);

/*
	Transpose<Blackbox> B(&A);

	pass = pass && testBlackbox(B);

	Vector x(n), y(n), z(n);
	int k = (n > 5 ? 5 : n);
	for(int i = 0; i < k; ++i) F.init(x[i], i);

	VectorDomain<Field> VD(F);
	A.apply(y, x);
	B.applyTranspose(z, x);
	pass = pass && VD.areEqual(y, z);

	B.apply(y, x);
	A.applyTranspose(z, x);
	pass = pass && VD.areEqual(y, z);
*/

	commentator().stop("TriplesBB black box test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

