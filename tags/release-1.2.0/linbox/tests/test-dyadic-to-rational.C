/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* tests/test-template.C
 * LinBox copyright/license applies. See COPYING.
 */

#include "linbox/linbox-config.h"
#include "linbox/util/commentator.h"
#include "test-common.h"

#define TestingDyadicToRational 
#include "linbox/algorithms/dyadic-to-rational.h" // can assume commentator and config loaded
#undef TestingDyadicToRational 

int main (int argc, char **argv) { 

	// customize optional args
	size_t n = 10;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.",                 TYPE_INT,     &n },
		END_OF_ARGUMENTS
	};
	parseArguments (argc, argv, args);

	commentator.getMessageClass (TIMING_MEASURE).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (10);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	// run test
	LinBox::commentator.start("DyadicToRational unit/regression test");
	bool pass = testDyadicToRational(n);
	LinBox::commentator.stop("DyadicToRational unit/regression test");
	return pass ? 0 : -1;

	/* parseArguments sets flags n and q, and arg structure ac, av, where ac, av refer to remaining args after flag processing.
	 Any -n <n>, -q <q>, -h flag, and first non-option arg will have been stripped from av.
	 The first non-option arg is used as file name for commentator output.
	 from av.  Ac is argc less the arguments thus stripped and handled.
	 -h, by convention sets ac = -1, indicating that testOBJECT should just explain the args on cout.

	 note: options are deliberately limited to these three: 
	  -h is conventional for giving a brief help with args (and is widely used).
	  -n sets a size parameter, a frequent device.
	  -q often determines the field of computation but may be used in other ways.

	 Generally this is as many run time parameters as a test user can reasonably remember and exploit.  
	 If additional control of testing is desired, usually the easiest and least error prone method 
	 is to edit the test function and recompile.  However, the test function can use the ac and av args 
	 to exploit further parameters. This ac/av use discouraged absent strong reason for it.

	 One use pattern of tests is that the unit designer tests with many parameters during development.
	 These may turn levels of debugging on and off, explore special cases, etc.).  However, when the unit 
	 is deployed and the test serves only as a unit/regression test, these features are seldom needed, 
	 if ever.  They may be left there just in case, but it is suggested NOT to document them in the -h help.
	*/

}
