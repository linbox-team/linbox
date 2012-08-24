
/* tests/test-modular.C
 * Copyright (C) 2001, 2002 Bradford Hovinen,
 * Copyright (C) 2002 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Rename from test-large-modular.C to test-modular.C; made other updates in
 * accordance with changes to Modular interace.
 * ------------------------------------
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 *.
 */

/*! @file   tests/test-modular-int.C
 * @ingroup tests
 * @brief run runFieldTests testRandomIterator tests on modular-int32_t
 * @test run runFieldTests testRandomIterator tests on modular-int32_t
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <queue>

#include "linbox/field/modular.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static integer q1("18446744073709551557");
	static integer q2 = 1073741789;
	static integer q3 = 65521U;
	static int q4 = 101;
	static integer q5 = 144115188075855881;
	static size_t n = 10000;
	static int iterations = 1;
	static int trials = 10000;
	static int categories = 1000;
	static int hist_level = 10;

	static Argument args[] = {
		{ 'K', "-K Q", "Operate over the \"field\" GF(Q) [1] for integer modulus.", TYPE_INTEGER, &q1 },
		{ 'k', "-k Q", "Operate over the \"field\" GF(Q) [1] for int64_t modulus.", TYPE_INTEGER, &q5 },
		{ 'Q', "-Q Q", "Operate over the \"field\" GF(Q) [1] for int32_t modulus.", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for int16_t modulus.", TYPE_INTEGER, &q3 },
		{ 'p', "-p P", "Operate over the \"field\" GF(Q) [1] for int8_t modulus.", TYPE_INT, &q4 },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	bool pass = true;

	{
		bool part_pass = true ;
		commentator().start("Modular<integer> field test suite", "Modular<integer>");

		Modular<integer> F_int (q1);
		// integer k = FieldTraits<Modular<integer> >::maxModulus() ;
		// prevprime(k,k);
		Modular<integer> I_int(q5);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<integer>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<integer>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<integer>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<integer>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<integer> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Modular<int32_t> field test suite", "Modular<int32_t>");

		Modular<int32_t> F_int (q2);
		integer k = FieldTraits<Modular<int32_t> >::maxModulus() ;
		prevprime(k,k);
		Modular<int32_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<int32_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<int32_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<int32_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<int32_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<int32_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Modular<int64_t> field test suite", "Modular<int64_t>");
		Modular<int64_t> F_int (q5);
		integer k = FieldTraits<Modular<int64_t> >::maxModulus() ;
		prevprime(k,k);
		Modular<int64_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<int64_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<int64_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<int64_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<int64_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<int64_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Modular<uint32_t> field test suite", "Modular<uint32_t>");

		Modular<uint32_t> F_int (q2);
		integer k = FieldTraits<Modular<uint32_t> >::maxModulus() ;
		prevprime(k,k);
		Modular<uint32_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<uint32_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<uint32_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<uint32_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<uint32_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<uint32_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Modular<uint64_t> field test suite", "Modular<uint64_t>");
		Modular<uint64_t> F_int (q5);
		integer k = FieldTraits<Modular<uint64_t> >::maxModulus() ;
		prevprime(k,k);
		Modular<uint64_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<uint64_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<uint64_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<uint64_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<uint64_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<uint64_t> field test suite");
		pass &= part_pass ;
	}


	{
		bool part_pass = true ;
		commentator().start("Modular<int8_t> field test suite", "Modular<int8_t>");

		Modular<int8_t> F_int (q4);
		integer k = FieldTraits<Modular<int8_t> >::maxModulus() ;
		prevprime(k,k);
		Modular<int8_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<int8_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<int8_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<int8_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<int8_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<int8_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Modular<int16_t> field test suite", "Modular<int16_t>");
		Modular<int16_t> F_int (q3);
		integer k = FieldTraits<Modular<int16_t> >::maxModulus() ;
		prevprime(k,k);
		Modular<int16_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<int16_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<int16_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<int16_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<int16_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<int16_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Modular<uint8_t> field test suite", "Modular<uint8_t>");

		Modular<uint8_t> F_int (q4);
		integer k = FieldTraits<Modular<uint8_t> >::maxModulus() ;
		prevprime(k,k);
		Modular<uint8_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<uint8_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<uint8_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<uint8_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<uint8_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<uint8_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Modular<uint16_t> field test suite", "Modular<uint16_t>");
		Modular<uint16_t> F_int (q3);
		integer k = FieldTraits<Modular<uint16_t> >::maxModulus() ;
		prevprime(k,k);
		Modular<uint16_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Modular<uint16_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Modular<uint16_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Modular<uint16_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Modular<uint16_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Modular<uint16_t> field test suite");
		pass &= part_pass ;
	}


	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

