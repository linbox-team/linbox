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


#include "linbox/field/modular.h"
#include "linbox/blackbox/triplesbb.h"
#include "test-blackbox.h"
#include "test-common.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	// ofstream report;

	bool pass = true;

	static size_t m = 4;
	static size_t n = 20;
	static size_t nnz = 0;
	static integer q = 2147483647U;
	q = 101;

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of test matrix to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set col dimension of test matrix to N.", TYPE_INT,     &n },
		{ 'z', "-n NNZ", "Set number of nonzero entries in test matrix.", TYPE_INT,     &nnz },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	
	if (nnz == 0) nnz = m*n/10;

	srand ((unsigned)time (NULL));

	commentator().start("TriplesBB black box test suite", "triplesbb");

	//typedef Modular<uint32_t> Field;
	typedef Modular<double> Field;
	typedef Field::Element Element;
	typedef TriplesBB<Field> Blackbox;

	Field F (q);

	// set up the matrix
	Blackbox A(F, m, n);
	Element d;

	for(int i = 0; i < (int)nnz; ++i)
	{
		F.init(d, rand());
		size_t ii = rand()%m;
		size_t jj = rand()%n;
		A.setEntry(ii,jj, d);
	}

	pass = pass && testBlackbox(A);

	commentator().stop("TriplesBB black box test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

