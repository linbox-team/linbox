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
	Modular<uint32> F_uint32((uint32) q2);
	Modular<uint16> F_uint16 ((uint16) q3);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!testField                 (F_integer, "Testing Modular<integer> field"))      pass = false;
	if (!testFieldNegation         (F_integer, "Modular<integer>", iterations))        pass = false;
	if (!testFieldInversion        (F_integer, "Modular<integer>", iterations))        pass = false;
	if (!testFieldAxioms           (F_integer, "Modular<integer>", iterations))        pass = false;
	if (!testFieldAssociativity    (F_integer, "Modular<integer>", iterations))        pass = false;
	if (!testGeometricSummation    (F_integer, "Modular<integer>", iterations, 100))   pass = false;
	if (!testFreshmansDream        (F_integer, "Modular<integer>", iterations))        pass = false;
	if (!testArithmeticConsistency (F_integer, "Modular<integer>", iterations))        pass = false;
	if (!testAxpyConsistency       (F_integer, "Modular<integer>", iterations))        pass = false;
	if (!testFieldAXPY             (F_integer, n, iterations, "Testing Modular<integer> FieldAXPY")) pass = false;

	if (!testField                 (F_uint32,    "Testing Modular<uint32> field"))         pass = false;
	if (!testFieldNegation         (F_uint32,    "Modular<uint32>", iterations))           pass = false;
	if (!testFieldInversion        (F_uint32,    "Modular<uint32>", iterations))           pass = false;
	if (!testFieldAxioms           (F_uint32,    "Modular<uint32>", iterations))           pass = false;
	if (!testFieldAssociativity    (F_uint32,    "Modular<uint32>", iterations))           pass = false;
	if (!testGeometricSummation    (F_uint32,    "Modular<uint32>", iterations, 100))      pass = false;
	if (!testFreshmansDream        (F_uint32,    "Modular<uint32>", iterations))           pass = false;
	if (!testArithmeticConsistency (F_uint32,    "Modular<uint32>", iterations))           pass = false;
	if (!testAxpyConsistency       (F_uint32,    "Modular<uint32>", iterations))           pass = false;
	if (!testFieldAXPY             (F_uint32,    n, iterations, "Testing Modular<uint32> FieldAXPY")) pass = false;

	if (!testField                 (F_uint16,   "Testing Modular<uint16> field")) pass = false;
	if (!testFieldNegation         (F_uint16,   "Modular<uint16>", iterations)) pass = false;
	if (!testFieldInversion        (F_uint16,   "Modular<uint16>", iterations)) pass = false;
	if (!testFieldAxioms           (F_uint16,   "Modular<uint16>", iterations)) pass = false;
	if (!testFieldAssociativity    (F_uint16,   "Modular<uint16>", iterations)) pass = false;
	if (!testGeometricSummation    (F_uint16,   "Modular<uint16>", iterations, 100)) pass = false;
	if (!testFieldCharacteristic   (F_uint16,   "Modular<uint16>", iterations)) pass = false;
	if (!testFreshmansDream        (F_uint16,   "Modular<uint16>", iterations)) pass = false;
	if (!testArithmeticConsistency (F_uint16,   "Modular<uint16>", iterations)) pass = false;
	if (!testAxpyConsistency       (F_uint16,   "Modular<uint16>", iterations)) pass = false;
	if (!testFieldAXPY             (F_uint16,   n, iterations, "Testing Modular<uint16> FieldAXPY")) pass = false;

#if 0
	FieldArchetype K(new LargeModular(101));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of Modular field"))
		pass = false;
#endif

	return pass ? 0 : -1;
}
