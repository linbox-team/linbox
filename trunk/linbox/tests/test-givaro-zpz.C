
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

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/givaro-zpz.h"
#include "linbox/field/givaro-montg.h"
#include "linbox/field/givaro-gfq.h"
#include "linbox/field/givaro-extension.h"

#include "test-common.h"
#include "test-generic.h"

#define TEST_ARCHETYPES 1

using namespace LinBox;

int main (int argc, char **argv)
{
        static integer q = 10733;
	static size_t n = 10000;
	static int iterations = 10;
	static int e = 3;
	static int trials = 1000;
	static int categories = 100;
	static int hist_level = 1;


        static Argument args[] = {
                { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
                { 'e', "-e E", "Use GF(q^e) for the extension field [1].",  TYPE_INT,     &e },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
        { '\0' }
        };

        parseArguments (argc, argv, args);

	//cout << endl << "GivaroZpz<Std16> field test suite" << endl;
	//cout.flush ();
	bool pass = true;
	
	GivaroZpz<Std16> F1 (q);
	GivaroZpz<Std32> F2 (q);
//	GivaroZpz<Log16> F3 (q);
	GivaroMontg F3 (39989);
	GivaroGfq F4 (q, 1);
	GivaroGfq F5 (11, e);
 	GivaroExtension<GivaroGfq> F6 (F5, e );
 	GivaroExtension<> F7 (103, 4 );

	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F1, "GivaroZpz<Std16>", iterations, n, false)) pass = false;
	if (!runFieldTests (F2, "GivaroZpz<Std32>", iterations, n, false)) pass = false;
	if (!runFieldTests (F3, "GivaroMontg", iterations, n, false)) pass = false;
	if (!runFieldTests (F4, "GivaroGfq (prime)", iterations, n, false)) pass = false;
	if (!runFieldTests (F5, "GivaroGfq (simple extension)", iterations, n, false)) pass = false;
	if (!runFieldTests (F6, "GivaroExtension (small polynomial extension)", iterations, n, false)) pass = false;
	if (!runFieldTests (F7, "GivaroExtension (large polynomial extension)", iterations, n, false)) pass = false;


	if (!testRandomIterator (F1,  "GivaroZpz<Std16>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F2,  "GivaroZpz<Std32>", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F3,  "GivaroMontgomery", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F4,  "GivaroGfq (prime)", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F5,  "GivaroGfq (simple extension)", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F6,  "GivaroExtension (small polynomial extension)", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F7,  "GivaroExtension (large polynomial extension)", trials, categories, hist_level)) pass = false;



#if TEST_ARCHETYPES

	GivaroZpz<Std16> * K1g = new GivaroZpz<Std16> (101);
	FieldArchetype K1(K1g);
	if (!testField<FieldArchetype> (K1, "Testing archetype with envelope of GivaroZpz<Std16> field"))
		pass = false;
	delete K1g;
#endif

#if TEST_ARCHETYPES
	GivaroZpz<Std32> * K2g = new GivaroZpz<Std32>(101);
	FieldArchetype K2(K2g);

	if (!testField<FieldArchetype> (K2, "Testing archetype with envelope of GivaroZpz<Std32> field"))
		pass = false;
	delete K2g;
#endif

#if TEST_ARCHETYPES
	GivaroZpz<Log16> * K3g = new GivaroZpz<Log16>(101);
	FieldArchetype K3(K3g);

	if (!testField<FieldArchetype> (K3, "Testing archetype with envelope of GivaroZpz<Log16> field"))
		pass = false;
	delete K3g;
#endif

#if TEST_ARCHETYPES
	GivaroGfq * K4g = new GivaroGfq (101,1);
	FieldArchetype K4(K4g);

	if (!testField<FieldArchetype> (K4, "Testing archetype with envelope of GivaroGfq prime field"))
		pass = false;
	delete K4g;
#endif


	return pass ? 0 : -1;
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
