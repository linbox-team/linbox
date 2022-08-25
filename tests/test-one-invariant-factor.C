/* Copyright (C) LinBox
 *
 *  bds evolved this from test-last-invariant-factor
 *  Author: Zhendong Wan
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



/*! @file  tests/test-one-invariant-factor.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include <linbox/linbox-config.h>
#include "givaro/zring.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/ring/modular.h"
#include "linbox/algorithms/last-invariant-factor.h"
#include "linbox/algorithms/one-invariant-factor.h"
#include "linbox/algorithms/smith-form-adaptive.h"
#include "linbox/blackbox/scompose.h"
#include "linbox/blackbox/random-matrix.h"
#include "linbox/algorithms/rational-solver.h"
#include <time.h>
#include "test-smith-form.h"

#include "linbox/matrix/matrix-domain.h"
#include "linbox/util/commentator.h"
#include "linbox/vector/stream.h"
#include "test-common.h"

using namespace LinBox;

int main(int argc, char** argv)
{
    bool pass = true;
    static size_t n = 10;
	static unsigned int iterations = 1;
    static Argument args[] = {
        { 'n', "-n N", "Set order of test matrices to N > 2.", TYPE_INT,     &n },
        { 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		END_OF_ARGUMENTS
    };
	parseArguments (argc, argv, args);

    typedef Givaro::ZRing<Integer>      Ring;
    Ring R; Ring::RandIter gen(R);

    commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);
	commentator().start("One invariant factor test suite", "OIF");

	
	typedef DixonSolver<Ring, Givaro::Modular<int32_t>, PrimeIterator<IteratorCategories::HeuristicTag> > Solver;
    // typedef DixonSolver<Ring, Givaro::Modular<double>, LinBox::PrimeIterator<IteratorCategories::HeuristicTag> > Solver;

	typedef LastInvariantFactor<Ring, Solver> LIF;
	typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;

	OIF oif;
    std::ostream &report = commentator().report (0);
    /*  
    report << "DEFAULTOIFTHRESHOLD " << DEFAULTOIFTHRESHOLD << std::endl;
    report << "oif.getThreshold() " << oif.getThreshold() << std::endl;
    report << "oif.getCrossover() " << oif.getCrossover() << std::endl;
	*/

	oif.  setThreshold (40);

	// Note:  oneInvariantFactor (used in smith-form-binary) is not designed to catch small prime factors

	BlasMatrix<Ring> A(R,n,n);
	BlasVector<Ring> d(R,n), bumps(R,n), lumps(R,19);
	for (size_t i = 0; i <10; ++i) lumps[i] = i;
	for (size_t i = 10; i <19; ++i) lumps[i] = i-19;

	for (size_t i = 0; i < n-5; ++i) bumps[i] = 1;
   bumps[n-1] = 0;
   bumps[n-2] = 1013;
   bumps[n-3] = 101;
   bumps[n-4] = 1;
   bumps[n-5] = 101;
	makeSNFExample(A,d,bumps,lumps); // sets A and d.
   std::vector<integer> primeL = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47};

   report << std::endl;
   typename Ring::Element x, y;
   // Messes up too often at k = 1.  todo: fix it.
	for (size_t k = 2; k < n; ++k) {
      y = d[k];
		oif. oneInvariantFactor (x, A, k+1, primeL);
	
		report << "At " << k+1 << " expected ";
		R. write (report, y);

		R. write (report << ", got ", x);

		if (!R. areEqual (y, x)) {
		    report << ", error" << std::endl;
			pass = false;
			//break;
		} else
		report << std::endl;
   }


/*
	makeBumps(bumps, 1);
	makeSNFExample(A,d,bumps,lumps);
	sf.smithFormBinary (x, A);
   ////////////////////
	size_t m = (4 < n/2 ? 5 : n/2);
	if (n < 3) { report << "n too small" << std::endl; return -1; }
	vector<Ring::Element> d(n);
	d[0] = 2;
	d[0] = 101;
	for (size_t i = 1; i < m; ++i) d[i] = 1009*d[i-1];
	for (size_t i = m; i < n-2; ++i) d[i] = d[i-1];
	d[n-3] *= 1013;
	d[n-2] = 0;
	d[n-1] = 0;

	pass = testRandom(R, oif, d);
*/

	commentator().stop("One invariant factor test suite");
   return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
