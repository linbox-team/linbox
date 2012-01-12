/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* tests/test-blackbox-block-container.C
 *
 * Written by bds
 *
 * -----------------------------------------------------
 *
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


/*! @file  tests/test-blackbox-block-container.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>
#include <cstdio>

#include "linbox/algorithms/blackbox-block-container.h"

#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"

#include "test-common.h"

using namespace LinBox;


int main (int argc, char **argv)
{

	bool pass = true;

	static size_t n = 40;
	static size_t k = 4;
	static integer q = 65519U;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'k', "-k K", "Set dimension of blocks matrices to KxK.", TYPE_INT,     &k },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator.start("block container test", "bbbc");
	ostream& report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "over Modular<double>" << endl;
	typedef Modular<double> Field;
	typedef BlasMatrix<Field> Block;
	typedef BlasMatrix<Field> Blackbox;
	Modular<double> F (q);

	Blackbox A(F, n, n);
	Block U(F, k, n);
	Block V(F, n, k);
	Block W(F, k, k);

	// Solved a segfault in init.
	BlackboxBlockContainer<Field, Blackbox> B(&A, F, U, V);

	// A more thorough test should be constructed.

	if (pass) commentator.stop("block container test pass");
	else commentator.stop("block container test FAIL");
	return pass ? 0 : -1;
}
