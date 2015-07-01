
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
 * accordance with changes to Givaro::Modular interace.
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

#include <queue>

#include "linbox/ring/modular.h"

#include "test-common.h"
#include "test-generic.h"

#include "linbox/field/field-traits.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static integer q1("18446744073709551557");
	static integer q2 = 1073741789;
	static integer q3 = 65521U;
	static int q4 = 101;
	static integer q5("144115188075855881");
	static size_t n = 10000;
	static unsigned int iterations = 1;
	static unsigned int trials = 10000;
	static unsigned int categories = 1000;
	static unsigned int hist_level = 10;

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

	/******* blocks follow for these cases in this order
	 Modular<integer>
	 Modular<int32_t>
	 Modular<int64_t>
	 Modular<uint32_t>
	 Modular<uint64_t> 
	 Modular<int8_t>
	 Modular<int16_t>
	 Modular<uint8_t>
	 Modular<uint16_t>
	 ******/
	{ 
		bool part_pass = true ;
		commentator().start("Givaro::Modular<integer> field test suite", "Givaro::Modular<integer>");

		Givaro::Modular<integer> F_int (q1);
		// integer k = FieldTraits<Givaro::Modular<integer> >::maxModulus() ;
		// prevprime(k,k);
		Givaro::Modular<integer> I_int(q5);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<integer>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<integer>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<integer>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<integer>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<integer> field test suite");
		pass &= part_pass ;
	}

	{ 
		bool part_pass = true ;
		commentator().start("Givaro::Modular<int32_t> field test suite", "Givaro::Modular<int32_t>");

		integer qq = std::min(q2, FieldTraits<Givaro::Modular<int32_t> >::maxModulus());
		// if (isOdd(qq)) --qq;
		Givaro::prevprime(qq,qq);
		Givaro::Modular<int32_t> F_int (qq);
		integer k = FieldTraits<Givaro::Modular<int32_t> >::maxModulus() ;
		prevprime(k,k);
		Givaro::Modular<int32_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<int32_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<int32_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<int32_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<int32_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<int32_t> field test suite");
		pass &= part_pass ;
	}

	{ 
		bool part_pass = true ;
		commentator().start("Givaro::Modular<int64_t> field test suite", "Givaro::Modular<int64_t>");
		integer qq = FieldTraits<Givaro::Modular<int64_t> >::maxModulus()/2 ;
		prevprime(qq,qq);
		Givaro::Modular<int64_t> F_int (qq);
		integer k = FieldTraits<Givaro::Modular<int64_t> >::maxModulus() ;
		prevprime(k,k);
		Givaro::Modular<int64_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<int64_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<int64_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<int64_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<int64_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<int64_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Givaro::Modular<uint32_t> field test suite", "Givaro::Modular<uint32_t>");

		Givaro::Modular<uint32_t> F_int (q2);
		integer k = FieldTraits<Givaro::Modular<uint32_t> >::maxModulus() ;
		prevprime(k,k);
		Givaro::Modular<uint32_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<uint32_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<uint32_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<uint32_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<uint32_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<uint32_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Givaro::Modular<uint64_t> field test suite", "Givaro::Modular<uint64_t>");
		integer qq = FieldTraits<Givaro::Modular<uint64_t> >::maxModulus()/2 ;
		prevprime(qq,qq);
		Givaro::Modular<uint64_t> F_int (qq);
		integer k = FieldTraits<Givaro::Modular<uint64_t> >::maxModulus() ;
		prevprime(k,k);
		Givaro::Modular<uint64_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<uint64_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<uint64_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<uint64_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<uint64_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<uint64_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Givaro::Modular<int8_t> field test suite", "Givaro::Modular<int8_t>");

		Givaro::Modular<int8_t> F_int (q4);
		integer k = FieldTraits<Givaro::Modular<int8_t> >::maxModulus() ;
		prevprime(k,k);
		Givaro::Modular<int8_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<int8_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<int8_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<int8_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<int8_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<int8_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Givaro::Modular<int16_t> field test suite", "Givaro::Modular<int16_t>");
		integer qq = FieldTraits<Givaro::Modular<int16_t> >::maxModulus()/2 ;
		prevprime(qq,qq);

		Givaro::Modular<int16_t> F_int (qq);
		integer k = FieldTraits<Givaro::Modular<int16_t> >::maxModulus() ;
		prevprime(k,k);
		Givaro::Modular<int16_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<int16_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<int16_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<int16_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<int16_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<int16_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Givaro::Modular<uint8_t> field test suite", "Givaro::Modular<uint8_t>");

		Givaro::Modular<uint8_t> F_int ((uint32_t)q4);
		integer k = FieldTraits<Givaro::Modular<uint8_t> >::maxModulus() ;
		prevprime(k,k);
		Givaro::Modular<uint8_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<uint8_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<uint8_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<uint8_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<uint8_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<uint8_t> field test suite");
		pass &= part_pass ;
	}

	{
		bool part_pass = true ;
		commentator().start("Givaro::Modular<uint16_t> field test suite", "Givaro::Modular<uint16_t>");
		Givaro::Modular<uint16_t> F_int (q3);
		integer k = FieldTraits<Givaro::Modular<uint16_t> >::maxModulus() ;
		prevprime(k,k);
		Givaro::Modular<uint16_t> I_int(k);


		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		if (!runFieldTests (F_int,  "Givaro::Modular<uint16_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (F_int,  "Givaro::Modular<uint16_t>", trials, categories, hist_level)) part_pass = false;

		if (!runFieldTests (I_int,  "Givaro::Modular<uint16_t>",  iterations, n, false)) part_pass = false;
		if (!testRandomIterator (I_int,  "Givaro::Modular<uint16_t>", trials, categories, hist_level)) part_pass = false;


		commentator().stop(MSG_STATUS(part_pass),"Givaro::Modular<uint16_t> field test suite");
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

