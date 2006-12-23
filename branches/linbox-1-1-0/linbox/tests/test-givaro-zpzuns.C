/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-givaro-zpz.C
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
 * Boston, MA 02111-1307, USA.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/givaro-zpz.h"
#include "linbox/field/givaro-gfq.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
        static integer q = 10733;
	static size_t n = 10000;
	static int iterations = 10;
	static int e;

        static Argument args[] = {
                { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 10733)", TYPE_INTEGER, &q },
                { 'e', "-e E", "Use GF(q^e) for the extension field [1] (default 2)",  TYPE_INT,     &e },
		{ 'n', "-n N", "Set dimension of test vectors to NxN (default 10000)", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations (default 10)",      TYPE_INT,     &iterations },
                { '\0' }
        };

        parseArguments (argc, argv, args);

	cout << endl << "GivaroZpz<Unsigned32> field test suite" << endl;
	cout.flush ();
	bool pass = true;
	
	GivaroZpz<Unsigned32> F1 (2);
	GivaroZpz<Unsigned32> F2 (q);
	GivaroZpz<Unsigned32> F3 (3);
	GivaroZpz<Unsigned32> F4 (32749);
	GivaroZpz<Unsigned32> F5 (65521);

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F1, "2", iterations, n, false)) pass = false;
	if (!runFieldTests (F2, "10733", iterations, n, false)) pass = false;
	if (!runFieldTests (F3, "3", iterations, n, false)) pass = false;
	if (!runFieldTests (F4, "32749", iterations, n, false)) pass = false;
	if (!runFieldTests (F5, "65521", iterations, n, false)) pass = false;

#if 0

	FieldArchetype K1(new GivaroZpz<Std16> (101));
	cout<<"coco\n";
	if (!testField<FieldArchetype> (K1, "Testing archetype with envelope of GivaroZpz<Std16> field"))
		pass = false;
#endif

#if 0
	FieldArchetype K2(new GivaroZpz<Unsigned32>(101));
	if (!testField<FieldArchetype> (K2, "Testing archetype with envelope of GivaroZpz<Unsigned32> field"))
		pass = false;
#endif

#if 0
	FieldArchetype K3(new GivaroZpz<Log16>(101));

	if (!testField<FieldArchetype> (K3, "Testing archetype with envelope of GivaroZpz<Log16> field"))
		pass = false;
#endif

#if 0
	FieldArchetype K4(new GivaroGfq(101,1));

	if (!testField<FieldArchetype> (K4, "Testing archetype with envelope of GivaroGfq prime field"))
		pass = false;
#endif


	return pass ? 0 : -1;
}
