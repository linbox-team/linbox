/* tests/test-givaro-zpz.C
 * Copyright (C) 2002 Pascal Giorgi
 *
 * Written by Pascal Giorgi  <pascal.giorgi@ens-lyon.fr>
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


/*! @file  tests/test-givaro-zpz.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>


#include "givaro/modular.h"
#include "givaro/montgomery.h"
#include "givaro/gfq.h"

#include "test-common.h"
#include "test-field.h"

using namespace LinBox;

int main (int argc, char **argv)
{
    static integer q = 10733;
	static size_t n = 10000;
	static int iterations = 10;
	static int e = 3;
	static int trials = 10000;
	static int categories = 1000;
	static int hist_level = 10;


    static Argument args[] = {
        { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
        { 'e', "-e E", "Use GF(q^e) for the extension field [1].",  TYPE_INT,     &e },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
    };
    
    parseArguments (argc, argv, args);
    
        //cout << endl << "Givaro::Modular< int16_t> field test suite" << endl;
        //cout.flush ();
	bool pass = true;
    
        // Does not work with q > 256
    Givaro::Modular< int16_t> F1 ( (q<256?q:integer(101)) ); 
	Givaro::Modular< int32_t> F2 (q);
	Givaro::Montgomery< int32_t > F3 (39989);
	Givaro::GFqDom<int64_t> F4 (q, 1);
	Givaro::GFqDom<int64_t> F5 (11, e);
 	Givaro::Extension<Givaro::GFqDom<int64_t>> F6 (F5, e );
 	Givaro::Extension<> F7 (103, 4 );
        // Does not work with q > 256
    Givaro::Modular< Givaro::Log16> F8 ( (q<256?q:integer(101)) );
	Givaro::Modular< uint32_t> F1u (2);
	Givaro::Modular< uint32_t> F2u (q);
	Givaro::Modular< uint32_t> F3u (3);
	Givaro::Modular< uint32_t> F4u (32749);
	Givaro::Modular< uint32_t> F5u (65521);


	LinBox::commentator().start("Givaro-zpz test suite", "GivZpz");
	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F1, "Givaro::Modular< int16_t>", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F2, "Givaro::Modular< int32_t>", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F3, "Givaro::Montgomery< int32_t >", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F4, "Givaro::GFqDom<int64_t> (prime)", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F5, "Givaro::GFqDom<int64_t> (simple extension)", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F6, "Givaro::Extension (small polynomial extension)", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F7, "Givaro::Extension (large polynomial extension)", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F8, "GivaroLog13", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F1u, "Unsigned32-2",     (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F2u, "Unsigned32-q", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F3u, "Unsigned32-3",     (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F4u, "Unsigned32-32749", (unsigned int)iterations, n, false)) pass = false;
	if (!runFieldTests (F5u, "Unsigned32-65521", (unsigned int)iterations, n, false)) pass = false;


	if (!testRandomIterator (F1,  "Givaro::Modular< int16_t>", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F2,  "Givaro::Modular< int32_t>", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F3,  "Givaro::Montgomery< int32_t >", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F4,  "Givaro::GFqDom<int64_t> (prime)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F5,  "Givaro::GFqDom<int64_t> (simple extension)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F6,  "Givaro::Extension (small polynomial extension)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F7,  "Givaro::Extension (large polynomial extension)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F1u,  "Givaro::Modular< uint32_t>(2)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F2u,  "Givaro::Modular< uint32_t>(q)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F3u,  "Givaro::Modular< uint32_t>(3)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F4u,  "Givaro::Modular< uint32_t>(32749)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;
	if (!testRandomIterator (F5u,  "Givaro::Modular< uint32_t>(65521)", (unsigned int)trials, (unsigned int)categories, (unsigned int)hist_level)) pass = false;

	LinBox::commentator().stop(MSG_STATUS (pass), "Givaro::Modular test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
