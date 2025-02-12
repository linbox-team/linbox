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

#include <linbox/linbox-config.h>
#include <givaro/modular-balanced.h>
#include <givaro/givintprime.h>

#include "linbox/linbox-config.h"
#include "linbox/util/commentator.h"
#include "linbox/util/debug.h"


#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

/*! @bug testRandomIterator fails
 */
int main (int argc, char **argv)
{
	static size_t n = 1000;
	static unsigned int iterations = 1;
	static int trials = 10000;
	static int categories = 1000;
	static int hist_level = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("Givaro::ModularBalanced<float> field test suite", "Givaro::ModularBalanced<float>");
	bool pass = true;

	Givaro::IntPrimeDom IPD;
	integer k; 
	IPD.nextprime(k,Givaro::ModularBalanced<float>::minCardinality()-1);
	Givaro::ModularBalanced<float> Fmin (k);
	Givaro::ModularBalanced<float> F5 (5);
	Givaro::ModularBalanced<float> F7 (7);

	IPD.prevprime(k,Givaro::ModularBalanced<float>::maxCardinality()+1);
	Givaro::ModularBalanced<float> Fmax(k);
	IPD.prevprime(k,k/2);
	Givaro::ModularBalanced<float> Fmid(k);

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator().report() << "Fmin modulus is " << Fmin.cardinality() << std::endl;
	if (!runFieldTests (Fmin,  "Givaro::ModularBalanced<float> Fmin",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (Fmin,  "Givaro::ModularBalanced<float>", trials, categories, hist_level)) pass = false;

	commentator().report() << "F5 modulus is " << F5.cardinality() << std::endl;
	if (!runFieldTests (F5,  "Givaro::ModularBalanced<float> F5",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (F5,  "Givaro::ModularBalanced<float>", trials, categories, hist_level)) pass = false;

	commentator().report() << "F7 modulus is " << F7.cardinality() << std::endl;
	if (!runFieldTests (F7,  "Givaro::ModularBalanced<float> F7",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (F7,  "Givaro::ModularBalanced<float>", trials, categories, hist_level)) pass = false;

	commentator().report() << "Fmax modulus is " << Fmax.cardinality() << std::endl;
	commentator().report() << "maxCardinality is " << Givaro::ModularBalanced<float>::maxCardinality() << std::endl;
	if (!runFieldTests (Fmax,  "Givaro::ModularBalanced<float> Fmax",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (Fmax,  "Givaro::ModularBalanced<float>", trials, categories, hist_level)) pass = false;

	commentator().report() << "Fmid modulus is " << Fmid.cardinality() << std::endl;
	if (!runFieldTests (Fmid,  "Givaro::ModularBalanced<float> Fmid",  iterations, n, false)) pass = false;
	// if (!testRandomIterator (Fmid,  "Givaro::ModularBalanced<float>", trials, categories, hist_level)) pass = false;


	commentator().stop("Givaro::ModularBalanced<float> field test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
