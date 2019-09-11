/* Copyright (C) LinBox
 *
 *  Author: Zhendong Wan
 *  mods: bds
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


/*! @file  tests/test-smith-form-adaptive.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */




#include <linbox/linbox-config.h>

#include "givaro/zring.h"
#include <time.h>
#include "linbox/randiter/random-prime.h"
#include "linbox/util/commentator.h"
#include "linbox/vector/stream.h"
#include "linbox/algorithms/smith-form-adaptive.h"
using namespace LinBox; 

#include "test-smith-form.h"

int main(int argc, char** argv)
{

	bool pass = true;
	static size_t m = 3;
	static size_t n = 35;
	static Argument args[] = {
		{ 'm', "-m M", "Set order of test matrices to M.", TYPE_INT,  &m },
		{ 'n', "-n N", "Set order of test matrices to N.", TYPE_INT,  &n },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("Smith form adaptive algorithm test suite", "EGV++");

	//!@bug should be tried on NTZ_LL too
	typedef Givaro::ZRing<Integer> PIR;
	PIR R;

	size_t k = std::min(m,n);
	DenseMatrix<PIR> A(R,m,n);
	BlasVector<PIR> d(R,k), x(R,k), bumps(R,k), lumps(R,19);
	for (uint64_t i = 0; i <10; ++i) lumps[i] = i;
	for (uint64_t i = 10; i <19; ++i) lumps[i] = i-19;

	makeBumps(bumps, 0);
	makeSNFExample(A,d,bumps,lumps);
	SmithFormAdaptive::smithForm (x, A);
	pass = pass and checkSNFExample(d,x);

	makeBumps(bumps, 1);
	makeSNFExample(A,d,bumps,lumps);
	SmithFormAdaptive::smithForm (x, A);
	pass = pass and checkSNFExample(d,x);

	makeBumps(bumps, 2);
	makeSNFExample(A,d,bumps,lumps);
	SmithFormAdaptive::smithForm (x, A);
	pass = pass and checkSNFExample(d,x);

	makeBumps(bumps, 3);
	makeSNFExample(A,d,bumps,lumps);
	SmithFormAdaptive::smithForm (x, A);
	pass = pass and checkSNFExample(d,x);


	commentator().stop(MSG_STATUS(pass));
	return pass ? 0 : -1;

}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
