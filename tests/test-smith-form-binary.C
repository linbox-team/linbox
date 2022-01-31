/* Copyright (C) LinBox
 *
 *  Author: Zhendong Wan, mods -bds
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


/*! @file  tests/test-smith-form-binary.C
 * @ingroup tests
 * @brief Test the EGV divide and conquer SNF alg.
 */

#include "linbox/linbox-config.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/last-invariant-factor.h"
#include "linbox/algorithms/one-invariant-factor.h"
#include "linbox/algorithms/smith-form-binary.h"
#include "linbox/algorithms/smith-form-adaptive.h"
#include "linbox/blackbox/scompose.h"
#include "linbox/blackbox/random-matrix.h"
#include "linbox/algorithms/rational-solver.h"
//#include <time.h>
#include <givaro/modular.h>

#include "linbox/util/commentator.h"
//#include "linbox/vector/stream.h"
//#include "test-field.h"
using namespace LinBox;

#include "test-smith-form.h"

typedef Givaro::ZRing<Integer> PIR;
typedef Givaro::Modular<double> Field;
using PrimeGenerator = PrimeIterator<IteratorCategories::HeuristicTag>;
typedef DixonSolver<PIR, Field, PrimeGenerator> Solver;
typedef LastInvariantFactor<PIR, Solver> LIF;
typedef OneInvariantFactor<PIR, LIF, SCompose, RandomMatrix>  OIF;
typedef SmithFormBinary<PIR, OIF > SF;

template<typename Vect>
void spaceBumps(size_t m, size_t n, Vect& bumps, Vect& lumps) {
    size_t k = std::min(m,n);

    lumps.resize(19);

    for (size_t i = 0; i <10; ++i) lumps[i] = i;
    for (size_t i = 10; i <19; ++i) lumps[i] = i-19;
    for (size_t i = 0; i < k-5; ++i) bumps[i] = 1;
    bumps[k-1] = 0;
    bumps[k-2] = 1013;
    bumps[k-3] = 401;
    bumps[k-4] = 1;
    bumps[k-5] = 401;

//    // Fails on some makeSNFExample builds of A. as a hack, check against adaptive. Todo: fix it.
//    if (not pass) {
//	   DenseMatrix<PIR> B(R,m,n);
//	   makeSNFExample(B,d,bumps,lumps);
//	   sf.smithFormBinary (x, B);
//	   BlasVector<PIR> y(R,k);
//	   SmithFormAdaptive::smithForm (y, B);
//        pass = checkSNFExample(y,x);
//    }

}

template<typename Mat, typename Vect>
bool checkBumpsLumps(SF& sf,
                     Mat& A,
                     Vect& d, Vect& x,
                     const Vect& bumps, const Vect& lumps) {
    makeSNFExample(A,d,bumps,lumps);
    sf.smithFormBinary (x, A);
    return checkSNFExample(d,x);
}

int main(int argc, char** argv)
{
	static size_t m(10);
	static size_t n(10);

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of test matrices to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set col dimension  of test matrices to N.", TYPE_INT,     &n },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

//    commentator().setReportStream(std::clog);
//    commentator().setMaxDepth(-1);
//    commentator().setMaxDetailLevel(-1);

	commentator().start("SmithFormBinary test suite", "SmithFormBinary");

    PIR R;
    DenseMatrix<PIR> A(R,m,n);
    size_t k = std::min(m,n);
    DenseVector<PIR> d(R,k), x(R,k), bumps(R,k), lumps(R,k);


    SF sf;
    sf. setOIFThreshold (40);
    sf. setLIFThreshold (30);

	bool pass, gpass(true);

    {
        makeBumps(bumps, 1);
        pass = checkBumpsLumps(sf,A,d,x,bumps,lumps);
        gpass &= pass;
    }
    std::clog << "Bumps 1, \t" << (pass? "PASSED." : "ERROR.") << std::endl;

    {
        makeBumps(bumps, 2);
        pass = checkBumpsLumps(sf,A,d,x,bumps,lumps);
        gpass &= pass;
    }
    std::clog << "Bumps 2, \t" << (pass? "PASSED." : "ERROR.") << std::endl;

    {
        makeBumps(bumps, 3);
        pass = checkBumpsLumps(sf,A,d,x,bumps,lumps);
        gpass &= pass;
    }
    std::clog << "Bumps 3, \t" << (pass? "PASSED." : "ERROR.") << std::endl;

    {
        spaceBumps(m,n,bumps,lumps);
        pass = checkBumpsLumps(sf,A,d,x,bumps,lumps);
        gpass &= pass;
    }
    std::clog << "Space Bumps,\t" << (pass? "PASSED." : "ERROR.") << std::endl;

    commentator().stop("SmithFormBinary test suite");
    return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
