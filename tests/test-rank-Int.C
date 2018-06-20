/* tests/test-rank-Int.c
 * Time-stamp: <08 Aug 14 07:51:49 Jean-Guillaume.Dumas@imag.fr>
 * -----------------------------------------------------
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


/*! @file  tests/test-rank.C
 * @ingroup tests
 * @brief  no doc
 * @test no doc.
 */

#include "linbox/linbox-config.h"
#include <givaro/modular-integer.h>
#include "test-rank.h"

int main (int argc, char **argv)
{

	bool pass = true;

	static size_t n = 40;
	static integer bigQ("12345678901234567890123456789012345678901234568119");
	static int iterations = 1;
        static double sparsity = 0.05;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 's', "-s S", "Sparse matrices with density S.", TYPE_DOUBLE,     &sparsity },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	srand ((unsigned)time (NULL));
	// srand48 ((unsigned)time (NULL));

	commentator().start("Integer sparse matrix rank test suite", "rank");
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_NORMAL);

	// Givaro::Modular<Integer> Gq(bigQ);
	Givaro::Modular<integer> Gq(bigQ);
	pass = pass && testSparseRank(Gq,n,n+1,(size_t)iterations,sparsity);
	pass = pass && testSparseRank(Gq,LINBOX_USE_BLACKBOX_THRESHOLD+n,LINBOX_USE_BLACKBOX_THRESHOLD+n-1,(size_t)iterations,sparsity);

	commentator().stop("Integer sparse matrix rank TEST suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
