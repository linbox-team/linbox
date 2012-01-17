/* Copyright (C) LinBox
 *
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


/*! @file  tests/test-modular-balanced-float.C
 * @ingroup tests
 * @brief  tests only runFieldTests on modular-balanced-float
 * @test   tests only runFieldTests on modular-balanced-float
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "linbox/field/modular-balanced.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/*! @bug testRandomIterator fails
 */
int main (int argc, char **argv)
{
	static size_t n = 1000;
	static int iterations = 1;
	static int trials = 100000;
	static int categories = 100;
	static int hist_level = 1;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator.start("ModularBalanced<float> field test suite", "ModularBalanced<float>");
	bool pass = true;

	ModularBalanced<float> F_int (4093);//(2147483629);//(2147483647);
	ModularBalanced<float> G_int (2011);
	ModularBalanced<float> H_int (3);
	integer k = FieldTraits<ModularBalanced<float> >::maxModulus() ;
	prevprime(k,k);
	ModularBalanced<float> I_int(k);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F_int,  "ModularBalanced<float>",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (F_int,  "ModularBalanced<float>", trials, categories, hist_level)) pass = false;

	if (!runFieldTests (G_int,  "ModularBalanced<float>",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (G_int,  "ModularBalanced<float>", trials, categories, hist_level)) pass = false;

	if (!runFieldTests (H_int,  "ModularBalanced<float>",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (H_int,  "ModularBalanced<float>", trials, categories, hist_level)) pass = false;

	if (!runFieldTests (I_int,  "ModularBalanced<float>",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (I_int,  "ModularBalanced<float>", trials, categories, hist_level)) pass = false;




	commentator.stop("ModularBalanced<float> field test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

