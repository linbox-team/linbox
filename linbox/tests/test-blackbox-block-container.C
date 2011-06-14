/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* tests/test-blackbox-block-container-base.C
 *
 * Written by bds
 *
 * -----------------------------------------------------
 *
 * This file is part of LinBox, licensed under the GNU Lesser General
 * Public License. See COPYING for more information.
 */


/*! @file  tests/test-blackbox-block-container-base.C
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

	commentator.start("blackbox block container test suite", "bbc");
	ostream& report = commentator.report (Commentator::LEVEL_NORMAL, INTERNAL_DESCRIPTION);
	report << "over Modular<double>" << endl;
	typedef Modular<double> Field;
	typedef BlasMatrix<double> Block;
	typedef BlasBlackbox<Field> Blackbox;
	Modular<double> F (q);

	Blackbox A(F, n, n);
	Block U(k, n);
	Block V(n, k);
	Block W(k, k);

	// Solved a segfault in init.
	BlackboxBlockContainer<Field, Blackbox> B(&A, F, U, V);

	// A more thorough test should be constructed.

	commentator.stop("blackbox block container test suite");
	return pass ? 0 : -1;
}
