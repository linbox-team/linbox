/* tests/test-direct-sum.C
 * Copyright (C) LinBox
 * Written by David Saunders
 *
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

/*! @file  tests/test-direct-sum.C
 * @ingroup tests
 * @brief no doc.
 * @test NO DOC
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/direct-sum.h"

#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;
	static int iterations1 = 100;
	static int iterations2 = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		{ 'i', "-i I", "Perform each test for I iterations.",   TYPE_INT,     &iterations1 },
		{ 'j', "-j J", "Apply test matrix to J vectors.",         TYPE_INT,     &iterations2 },
		END_OF_ARGUMENTS
	};

	typedef Modular<uint32_t> Field;
	typedef vector<Field::Element> Vector;

	parseArguments (argc, argv, args);
	Field F (q);
	Field::Element k;

	commentator().start("DirectSum black box test suite", "direct sum");
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	F.init(k, 5);
	ScalarMatrix<Field> B(F, 10, k);
	F.init(k, 2);
	ScalarMatrix<Field> C(F, 5, k);

	DirectSum<ScalarMatrix<Field>, ScalarMatrix<Field> > A(&B, &C);
	pass = pass && testBlackbox(A);
	DirectSum<ScalarMatrix<Field>, ScalarMatrix<Field> > D(B, C);
	pass = pass && testBlackbox(D);

	commentator().stop("DirectSum black box test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

