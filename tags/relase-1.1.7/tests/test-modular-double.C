
/* tests/test-modular-double.C
 * Brazenly stolen by bds from 
 * tests/test-modular-short.C
 * Brazenly stolen by Zhendong Wan (Copyright (C) 2003) from 
 * tests/test-modular.C
 * Copyright (C) 2001, 2002 Bradford Hovinen,
 * Copyright (C) 2002 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Rename from test-large-modular.C to test-modular.C; made other updates in
 * accordance with changes to Modular interace.
 * ------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "linbox/field/modular-double.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static size_t n = 10000;
	static int iterations = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	commentator.start("Modular<double> field test suite", "Modular<double>");
	bool pass = true;

	//Modular<double> F2 (2); 
	Modular<double> F3 (3); 
	Modular<double> F5 (5); 
	Modular<double> F7 (7); 
	Modular<double> F11 (11); 
	Modular<double> F (32749); 
	Modular<double> G (65521); 
	//Modular<double> H (1099511627689); 

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
	std::ostream& report = commentator.report();
	report << "Field F2" << std::endl;
	//if (!runFieldTests (F2,  "Modular<double>",  iterations, n, false)) pass = false;
	report << "Field F3" << std::endl;
	if (!runFieldTests (F3,  "Modular<double>",  iterations, n, false)) pass = false;
	report << "Field F5" << std::endl;
	if (!runFieldTests (F5,  "Modular<double>",  iterations, n, false)) pass = false;
	report << "Field F7" << std::endl;
	if (!runFieldTests (F7,  "Modular<double>",  iterations, n, false)) pass = false;
	report << "Field F11" << std::endl;
	if (!runFieldTests (F11,  "Modular<double>",  iterations, n, false)) pass = false;
	report << "Field F" << std::endl;
	if (!runFieldTests (F,  "Modular<double>",  iterations, n, false)) pass = false;
	report << "Field G" << std::endl;
	if (!runFieldTests (G,  "Modular<double>",  iterations, n, false)) pass = false;

	commentator.stop("Modular<double> field test suite");
	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
