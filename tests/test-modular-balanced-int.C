/* tests/test-modular-balanced-int.C
 * Copyright (C) 2001, 2002 Bradford Hovinen,
 * Copyright (C) 2002, 2015 Dave Saunders
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>,
 *            Dave Saunders <saunders@cis.udel.edu>
 *
 * ------------------------------------
 * 2002-04-10 Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Rename from test-large-modular.C to test-modular.C; made other updates in
 * accordance with changes to Givaro::ModularBalanced interace.
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

/*! @file   tests/test-modular-balanced-int.C
 * @ingroup tests
 * @brief For each integer type T, Givaro::ModularBalanced<T> is tested with a small primm and with a large prime using runFieldTests and testRandomIterator.
 */

#include <linbox/linbox-config.h>

#include "givaro/modular-balanced.h"
#include "givaro/givintprime.h"

#include "test-field.h"
using namespace LinBox;

template<class int_type>
bool launchTests(string int_type_name, integer q, size_t n, uint32_t trials, uint32_t categories, uint32_t hist_level) 
	{ 
		using Field = Givaro::ModularBalanced<int_type>;
		bool pass = true ;

		string field_name = "Givaro::ModularBalanced<" + int_type_name + ">";
		string title = field_name + " field test suite";

		commentator().start(title.c_str(), field_name.c_str());
		// Make sure some more detailed messages get printed
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
		commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

		commentator().report() << field_name << " range: [" << Field::minCardinality() << "," << Field::maxCardinality() << "]" << std::endl;

		Field FS(q);
		pass &= runFieldTests (FS,  field_name.c_str(), 1, n, false);
		pass &= testRandomIterator (FS, field_name.c_str(), trials, categories, hist_level);

		Givaro::IntPrimeDom IPD;
		integer k = Field::maxCardinality();
		if (k > 0) {
			IPD.prevprime(k,k + 1);
		} else {
			static integer qi("18446744073709551557");
			k = qi;
		}
		Field FL(k);
		pass &= runFieldTests (FL,  field_name.c_str(), 1, n, false);
		pass &= testRandomIterator (FL, field_name.c_str(), trials, categories, hist_level);

		commentator().stop(MSG_STATUS(pass), "field test-suite");
		return pass;
	}

int main (int argc, char **argv)
{
	// for field testing
	static integer q = 5; // small prime valid for all int types.

	// for randiter testing
	static size_t n = 10000;
	static unsigned int trials = 10000;
	static unsigned int categories = 1000;
	static unsigned int hist_level = 10;

	static Argument args[] = {
		{ 'K', "-K Q", "to use ModularBalanced<T> F(q). Must be valid for all int types T. (A large prime is also used.)", TYPE_INTEGER, &q },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	bool pass = true;


	pass &= launchTests<int64_t>("int64_t",q,n,trials,categories,hist_level);
	pass &= launchTests<int32_t>("int32_t",q,n,trials,categories,hist_level);

/*  these are not in givaro/modular-balanced.h
	pass &= launchTests<integer>("integer",q,n,trials,categories,hist_level);
	pass &= launchTests<uint64_t>("uint64_t",q,n,trials,categories,hist_level);
	pass &= launchTests<uint32_t>("uint32_t",q,n,trials,categories,hist_level);
	pass &= launchTests<int16_t>("int16_t",q,n,trials,categories,hist_level);
	pass &= launchTests<int8_t>("int8_t",q,n,trials,categories,hist_level);
	pass &= launchTests<uint16_t>("uint16_t",q,n,trials,categories,hist_level);
	pass &= launchTests<uint8_t>("uint8_t",q,n,trials,categories,hist_level);
*/

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
