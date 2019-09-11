/* Copyright (C) LinBox
 *
 *  Author: Zhendong Wan
 *  Mods: bds
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

/*! @file tests/test-smith-form.C
 * @ingroup tests
 * @brief no doc.
 * @test no doc.
	//!@bug should work for NTL Integers too
 */


#include <linbox/linbox-config.h>
#include "linbox/solutions/smith-form.h"

#include "givaro/zring.h"
#include "linbox/util/commentator.h"
#include "linbox/matrix/dense-matrix.h"
#include "linbox/vector/blas-vector.h"
using namespace LinBox;

#include "test-smith-form.h"

int main(int argc, char** argv)
{

	bool pass = true;
	static size_t m =30;
	static size_t n =20;
	static int iterations = 2;
	static Argument args[] = {
		{ 'm', "-n M", "Set row dim of test matrices to N.", TYPE_INT,  &m },
		{ 'n', "-n N", "Set col dim of test matrices to N.", TYPE_INT,  &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT, &iterations },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);
	//!@bug should be tried on NTZ_LL too
	typedef Givaro::ZRing<Integer> PIR;
	PIR R;

	commentator().start("Smith form test", "Smith");

	size_t k = std::min(m,n);
	DenseMatrix<PIR> A(R,m,n);
	BlasVector<PIR> d(R,k), x(R,k), bumps(R,k), lumps(R,19);
	for (uint64_t i = 0; i <10; ++i) lumps[i] = i;
	for (uint64_t i = 10; i <19; ++i) lumps[i] = i-19;

	makeBumps(bumps, 0);
	makeSNFExample(A,d,bumps,lumps);
	smithForm (x, A);
	pass = pass and checkSNFExample(d,x);

	makeBumps(bumps, 1);
	makeSNFExample(A,d,bumps,lumps);
	smithForm (x, A);
	pass = pass and checkSNFExample(d,x);

	makeBumps(bumps, 2);
	makeSNFExample(A,d,bumps,lumps);
	smithForm (x, A);
	//SmithFormAdaptive::compute_local_long(x, A, 2, 64);
	pass = pass and checkSNFExample(d,x);

	makeBumps(bumps, 3);
	makeSNFExample(A,d,bumps,lumps);
	smithForm (x, A);
	pass = pass and checkSNFExample(d,x);

	commentator().stop("Smith form test");
	return pass ? 0 : -1;

}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
