/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

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
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

template <class Field>
bool runFieldTests (const Field &F, const char *desc, unsigned int iterations, size_t n, bool runCharacteristicTest) 
{
	bool pass = true;
	ostringstream str1, str2;

	str1 << "Testing " << desc << " field" << ends;
	str2 << "Testing " << desc << " FieldAXPY" << ends;

	commentator.start (str1.str ().c_str (), "runFieldTests", runCharacteristicTest ? 11 : 10);

	if (!testField                 (F, str1.str ().c_str ()))                pass = false; commentator.progress ();
	if (!testFieldNegation         (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldInversion        (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldAxioms           (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldAssociativity    (F, desc, iterations))                    pass = false; commentator.progress ();

	if (runCharacteristicTest) {
		if (!testFieldCharacteristic (F, desc, iterations))
			pass = false;

		commentator.progress ();
	}

	if (!testGeometricSummation    (F, desc, iterations, 100))               pass = false; commentator.progress ();
	if (!testFreshmansDream        (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testArithmeticConsistency (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testAxpyConsistency       (F, desc, iterations))                    pass = false; commentator.progress ();
	if (!testFieldAXPY             (F, n, iterations, str2.str ().c_str ())) pass = false; commentator.progress ();

	commentator.stop (MSG_STATUS (pass), (const char *) 0, "runFieldTests");

	return pass;
}

int main (int argc, char **argv)
{
	static integer q1("18446744073709551557");
	static integer q2 = 2147483647U;
	static integer q3 = 65521U;
	static size_t n = 10000;
	static int iterations = 10;

	static Argument args[] = {
		{ 'K', "-K Q", "Operate over the \"field\" GF(Q) [1] for integer modulus (default 18446744073709551557)", TYPE_INTEGER, &q1 },
		{ 'Q', "-Q Q", "Operate over the \"field\" GF(Q) [1] for uint32 modulus (default 2147483647)", TYPE_INTEGER, &q2 },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] for uint16 modulus (default 65521)", TYPE_INTEGER, &q3 },
		{ 'n', "-n N", "Set dimension of test vectors to NxN (default 10000)",      TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",           TYPE_INT,     &iterations },
		{ '\0' }
	};

	parseArguments (argc, argv, args);

	cout << "Modular field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;

	Modular<integer> F_integer (q1);
	Modular<uint32> F_uint32 ((uint32) q2);
	Modular<uint16> F_uint16 ((uint16) q3);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F_integer, "Modular<integer>", n, iterations, false)) pass = false;
	if (!runFieldTests (F_uint32,  "Modular<uint32>",  n, iterations, false)) pass = false;
	if (!runFieldTests (F_uint16,  "Modular<uint16>",  n, iterations, false)) pass = false;

#if 0
	FieldArchetype K(new LargeModular(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of Modular field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
