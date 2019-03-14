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


/*! @file  tests/test-givaro-zpzuns.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>


#include "givaro/modular.h"

#include "test-common.h"
#include "test-field.h"

using namespace LinBox;

int main (int argc, char **argv)
{
    static integer q = 10733;
	static size_t n = 10000;
	static unsigned int iterations = 10;
	static int e;
	static unsigned int trials = 10000;
	static unsigned int categories = 1000;
	static unsigned int hist_level = 10;


    static Argument args[] = {
        { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
        { 'e', "-e E", "Use GF(q^e) for the extension field [1].", TYPE_INT,     &e },
		{ 'n', "-n N", "Set dimension of test vectors to NxN.", TYPE_INT,     &n },
		{ 'i', "-i I", "Perform each test for I iterations.", TYPE_INT,     &iterations },
		{ 't', "-t T", "Number of trials for the random iterator test.", TYPE_INT, &trials },
		{ 'c', "-c C", "Number of categories for the random iterator test.", TYPE_INT, &categories },
		{ 'H', "-H H", "History level for random iterator test.", TYPE_INT, &hist_level },
		END_OF_ARGUMENTS
    };
    
    parseArguments (argc, argv, args);
    
	//cout << endl << "Givaro::Modular< uint32_t> field test suite" << endl;
	//cout.flush ();
	bool pass = true;

	Givaro::Modular< uint32_t> F1 (2);
	Givaro::Modular< uint32_t> F2 (q);
	Givaro::Modular< uint32_t> F3 (3);
	Givaro::Modular< uint32_t> F4 (32749);
	Givaro::Modular< uint32_t> F5 (65521);

	LinBox::commentator().start("Givaro-zpzuns test suite", "GivZpzu");
	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	if (!runFieldTests (F1, "2",     iterations, n, false)) pass = false;
	if (!runFieldTests (F2, "10733", iterations, n, false)) pass = false;
	if (!runFieldTests (F3, "3",     iterations, n, false)) pass = false;
	if (!runFieldTests (F4, "32749", iterations, n, false)) pass = false;
	if (!runFieldTests (F5, "65521", iterations, n, false)) pass = false;

	if (!testRandomIterator (F1,  "Givaro::Modular< uint32_t>(2)", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F2,  "Givaro::Modular< uint32_t>(10733)", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F3,  "Givaro::Modular< uint32_t>(3)", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F4,  "Givaro::Modular< uint32_t>(32749)", trials, categories, hist_level)) pass = false;
	if (!testRandomIterator (F5,  "Givaro::Modular< uint32_t>(65521)", trials, categories, hist_level)) pass = false;

	LinBox::commentator().stop(MSG_STATUS (pass), "Givaro::Modularuns test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
