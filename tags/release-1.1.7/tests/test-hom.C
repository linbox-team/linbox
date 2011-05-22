/* tests/test-hom.C
 * Copyright (C) LinBox
 *
 * Written by Dave Saunders <saunders@cis.udel.edu>
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include "linbox/field/modular.h"
#include "linbox/field/hom.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static integer q = 65521U;
	static size_t n = 10000;
	static int iterations = 10;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16 modulus.", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start("Hom test suite", "Hom");
	bool pass = true;

	Modular<uint32> F_uint32 ((uint32) q);
	Modular<uint16> F_uint16 ((uint16) q);
	Hom<Modular<uint16>, Modular<uint32> > iso(F_uint16, F_uint32);

	uint16 x=2, y;
	uint32 z=2, w;
	iso.image(w, x);
	pass = pass && F_uint32.areEqual(z, w);
	iso.preimage(y, z);
	pass = pass && F_uint16.areEqual(x, y);

	/* for image field!
	uint32 x, y, z, w;
	iso.smul(x, 2, 3);
	iso.mul(y, 2, 3);
	pass = pass && iso.areEqual(x, y);
	iso.smulin(x, 5);
	iso.mulin(y, 5);
	pass = pass && iso.areEqual(x, y);

	iso.saxpy(z, 7, x, 11);
	iso.axpy(w, 7, x, 11);
	pass = pass && iso.areEqual(z, w);
	iso.saxpyin(z, 7, x);
	iso.axpyin(w, 7, x);
	pass = pass && iso.areEqual(z, w);
	*/

	commentator.stop("Hom test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
