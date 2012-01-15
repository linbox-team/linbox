/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s


/* tests/test-template.C
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
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */


#include "linbox/linbox-config.h"
#include "linbox/util/commentator.h"
#include "test-common.h"

#include "linbox/algorithms/dyadic-to-rational.h" // can assume commentator and config loaded

using namespace LinBox;

#include <vector>
#include <cmath>
#include "linbox/field/PID-integer.h"
#include "linbox/util/timer.h"
#include "linbox/util/commentator.h"

int test1(size_t k, size_t dxa, size_t denBs)
{

/* Check reconstruction of i/k when it is well approximated
by something over dxa and when the denominator bound (for k) is denBs.
Values of numerator i from -k-2 to k+2 are checked.
When dxa^2 <= denBs, 2 should be returned.
*/
	typedef PID_integer Ring; Ring Z;
	typedef Ring::Element Int;

//std::cout << "test1(k " << k << ", dxa " << dxa << ", denBs " << denBs << ")" << std::endl;
	// k is our denominator
	Int denB = denBs; // denominator bound
	// approximate(i*k^2)/k for i = 0..k-1
	size_t kp = k+2;
	size_t kp2 = 2*kp;
	std::vector<Int> nx(kp2);
	for (size_t i = 0; i < kp2 ; ++i)
		Z.init(nx[i],floor(((double(i)-double(kp))*(double)dxa)/(double)k + 0.5));
	// |nx[i]/dxa - (i-kp)/k| <= 1/2dxa
	std::vector<Int> n(kp2);
	Int d;
	bool pass = true;
	bool claim;
	int ret = 2;

	// check of individual reconstructions
	Int dx; Z.init(dx, dxa);
	int c;
	// individual reconstructions
	for (size_t i = 0; i < kp2 ; ++i) {
	    bool loopclaim = true;
//std::cout << nx[i] << " " << dx << " " << denB << std::endl;
		c = dyadicToRational(Z, n[i], d, nx[i], dx, denB);
//std::cout << " c " << c << " n " << n[i] << " d " << d << " nx " << nx[i] << " dx " << dx << " denB " << denB << std::endl;
		loopclaim = ((c > 0)  && (n[i]*k == ((int(i)-int(kp))*d)));
		if ( c == 2 ) loopclaim = loopclaim && d*denB < dx;
		if (c < ret) ret = c;
		if (! loopclaim) ret = 0;
//if (! claim) std::cout << "F " << pass << claim << ret << std::endl;
		//if (! loopclaim)
//std::cout << "F2 " << loopclaim << " i " << i << " nx/dx " << nx[i] << "/" << dx << ", n/d " << n[i] << "/" << d << std::endl;
		pass = pass && loopclaim;
	}
//std::cout << "result, pass " << pass << " ret " << ret << std::endl;

#if 1
	// check vector reconstruction
	c = dyadicToRational(Z, n, d, nx, dx, denB);
	commentator.report() << "In test1 dyadicToRational returns " << c << std::endl;
	claim = 0 < c;
	if (claim) {
		for (size_t i = 0; i < k ; ++i) claim = claim && (n[i] == (int(i)-int(kp)));
	}
	pass = pass && claim;
	if (!claim) {
		commentator.report() << "first example fails" << std::endl;
		commentator.report() << "data for first example" << std::endl;
		for (size_t i = 0; i < 10 ; ++i)
			commentator.report() << nx[i] << std::endl;
		commentator.report() << dx << " den in" << std::endl;
		commentator.report() << "results for first example" << std::endl;
		for (size_t i = 0; i < 10 ; ++i)
			commentator.report() << n[i] << std::endl;
		commentator.report() << d << " den out" << std::endl;
	}

	if (c < ret) ret = c;
	if (! claim) ret = 0;
	pass = pass && claim;
//std::cout << "v " << pass << claim << ret << std::endl;
#endif

	//return pass;
	return ret;
}

bool testDyadicToRational(size_t k = 10, bool benchmarking = false)
{
	typedef PID_integer Ring; Ring Z;
	typedef Ring::Element Int;
	bool pass = true;
	bool claim = false;
	size_t pow2 = 1; // upperbound for k.
	for (size_t i = k; i > 0; i >>= 1) pow2 *= 2;

	UserTimer clock;
	double totaltime = 0;

	clock.clear(); clock.start();
	claim = 1 <= test1(k, pow2*pow2, k); // some 1's and some 2's
	if (!claim) commentator.report() << "failure: 1 !=test1(k, k*k, k)" << std::endl;
	pass = pass && claim;
	claim = 1 == test1(k, k*k, k*k); // all 1's
	if (!claim) commentator.report() << "failure: 1 !=test1(k, k*k, k*k)" << std::endl;
	pass = pass && claim;
	claim = 2 == test1(k, k*k+2*k + 1, k+1); // all 2's
	if (!claim) commentator.report() << "failure: 2 !=test1(k, (k + 1)^2, k+1)" << std::endl;
	pass = pass && claim;
	clock.stop(); totaltime += clock.time();

#if 1
// special case 1
	Int B; Z.init(B, 1000000000);
	Int B2; Z.init(B2); B2 = B*B;
	Int denB = 4*B+294967296; // 2^32

	Int d;
	Int dxs = denB*denB; // 2^64
	size_t k2 = 10;
	std::vector<Int> nx;
	std::vector<Int> n;
	nx.resize(k2);
	n.resize(k2);

	nx[0] =-143*B2-298423624*B-962150784;
	nx[1] = 239*B2+120348615*B+509085366;
	nx[2] =  -4*B2-959983787*B-562075119;
	nx[3] =  27*B2+  8864641*B+551149627;
	nx[4] =  62*B2+971469325*B+838237476;
	nx[5] = 190*B2+559070838*B+297135961;
	nx[6] = 176*B2+172593329*B+811309753;
	nx[7] = -70*B2-861003759*B-845628342;
	nx[8] =-228*B2-416339507*B-338896853;
	nx[9] = -14*B2-398832745*B-762391791;

	claim = 0 < dyadicToRational(Z, n, d, nx, dxs, denB);

	if (!claim) commentator.report() << "in special case 1 failure claimed" << std::endl;

	pass = pass && claim;

	std::vector<Int> ntrue(k2);
	Int dentrue = 691617936;
    ntrue[0] = -5*B-372642434;
	ntrue[1] =  8*B+965263534;
	ntrue[2] =  -185963102;
	ntrue[3] =  1*B+ 12634812;
	ntrue[4] =  2*B+360969365;
	ntrue[5] =  7*B+144570919;
	ntrue[6] =  6*B+605183272;
	ntrue[7] = -2*B-656769182;
	ntrue[8] = -8*B-563941509;
	ntrue[9] =     -539850878;
	claim = (d == dentrue);
	for (size_t i = 0; i < k2 ; ++i) claim = claim && (n[i] == ntrue[i]);
	pass = pass && claim;

	if (!claim)
	{
	commentator.report() << "data for failing special case 1" << std::endl;
	for (size_t i = 0; i < k2 ; ++i)
		commentator.report() << nx[i] << std::endl;
	commentator.report() << dxs << " den in" << std::endl;
	commentator.report() << "results for failing special case 1" << std::endl;
	for (size_t i = 0; i < k2 ; ++i)
		commentator.report() << n[i] << std::endl;
	commentator.report() << d << " den out" << std::endl;
	}

#endif
#if 1
// case where false should be returned.
	denB = 1023*B+948656640;

	dxs = 4*B2+611686018*B+427387904;
	size_t k3 = 10;
	nx.resize(k3);
	n.resize(k3);

	nx[0] =  -4*B2-474720817*B-626139615;
	nx[1] =  -9*B2-632772311*B-132715070;
	nx[2] = -19*B2-805041562*B-739831073;
	nx[3] =  35*B2+521355378*B+297487606;
    nx[4] =  27*B2+922294617*B+624925795;
    nx[5] =   1*B2+494454325*B+592253092;
    nx[6] = -27*B2-985233904*B-197462327;
    nx[7] = -20*B2-336729946*B-917106131;
    nx[8] = -42*B2-807924857*B-450940124;
    nx[9] = -27*B2-863559911*B-142533799;


	// this should fail
	claim = 2 > dyadicToRational(Z, n, d, nx, dxs, denB);
//std::cout << "d " << d << " dxs " << dxs << " denB " << denB << std::endl;

    pass = pass && claim;
	if (claim) commentator.report() << "third ratrecon falsely claims success" << std::endl;
#endif

// done
	if (benchmarking) commentator.report() << "vec size" << k << ", rat recon time: " << clock << " totaltime " << totaltime << std::endl;
	return pass;
}


int main (int argc, char **argv)
{

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
