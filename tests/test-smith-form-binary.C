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

int main(int argc, char** argv)
{
	bool pass = true;

	static size_t m =2;
	static size_t n =5;

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of test matrices to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set col dimension  of test matrices to N.", TYPE_INT,     &n },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("SmithFormBinary test suite", "SmithFormBinary");

	commentator().report() << std::endl << "EGV++ algorithm test suite with LinBox/Givaro ZRing:\n";
	{
//		typedef Givaro::IntegerDom Ring;
		typedef Givaro::ZRing<Integer> PIR;
		PIR R;

		typedef Givaro::Modular<int64_t> Field;
		typedef DixonSolver<PIR, Field, PrimeIterator<IteratorCategories::HeuristicTag> > Solver;
		typedef LastInvariantFactor<PIR, Solver> LIF;
		typedef OneInvariantFactor<PIR, LIF, SCompose, RandomMatrix>  OIF;
		typedef SmithFormBinary<PIR, OIF > SF;

		SF sf;
		sf. setOIFThreshold (30);
		sf. setLIFThreshold (30);

	size_t k = std::min(m,n);
	DenseMatrix<PIR> A(R,m,n);
	BlasVector<PIR> d(R,k), x(R,k), bumps(R,k), lumps(R,19);
	for (uint64_t i = 0; i <10; ++i) lumps[i] = i;
	for (uint64_t i = 10; i <19; ++i) lumps[i] = i-19;

	makeBumps(bumps, 0);
	makeSNFExample(A,d,bumps,lumps);
	sf.smithFormBinary (x, A);
	//pass = pass and checkSNFExample(d,x);
	pass = checkSNFExample(d,x) and pass;

	makeBumps(bumps, 1);
	makeSNFExample(A,d,bumps,lumps);
	sf.smithFormBinary (x, A);
	//pass = pass and checkSNFExample(d,x);
	pass = checkSNFExample(d,x) and pass;

	makeBumps(bumps, 2);
	makeSNFExample(A,d,bumps,lumps);
	sf.smithFormBinary (x, A);
	//pass = pass and checkSNFExample(d,x);
	pass = checkSNFExample(d,x) and pass;

	makeBumps(bumps, 3);
	makeSNFExample(A,d,bumps,lumps);
	sf.smithFormBinary (x, A);
	//pass = pass and checkSNFExample(d,x);
	pass = checkSNFExample(d,x) and pass;

	}

/* I don't think this section adds significant value, so deleting.
	{
		typedef Givaro::ZRing<Givaro::Integer> Ring;

		Ring R; Ring::RandIter gen(R);

		commentator().report() << std::endl << "EGV++ algorithm test suite with ZZ :\n";
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (5);

		//RandomDenseStream<Ring> s1 (R, gen, n, (unsigned int)iterations);
		RandomDenseStream<Ring> s1 (R, gen, n, 1);

		typedef Givaro::Modular<int32_t> Field;
		typedef DixonSolver<Ring, Field, PrimeIterator<IteratorCategories::HeuristicTag> > Solver;
		typedef LastInvariantFactor<Ring, Solver> LIF;
		typedef OneInvariantFactor<Ring, LIF, SCompose, RandomMatrix>  OIF;
		typedef SmithFormBinary<Ring, OIF, MatrixRank<Ring, Field > > SF;

		SF sf;
		sf. setOIFThreshold (30);
		sf. setLIFThreshold (30);

		if (!testRandom(R, sf, s1)) pass = false;
	}
*/

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
