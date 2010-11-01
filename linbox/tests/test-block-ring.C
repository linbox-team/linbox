/* Copyright (C) LinBox
 *
 * Copyright (C) 2007 b d saunders
 *
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


#include "linbox/field/modular.h"
#include "linbox/field/block-ring.h"

#include "test-common.h"
#include "test-field.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static size_t n = 10;
    static integer q = 10733;
	static int iterations = 2;

    static Argument args[] = {
	{ 'n', "-n N", "Set dimension of blocks to N.", TYPE_INT,     &n },
    { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
	{ 'i', "-i I", "Perform each test for I iterations.",      TYPE_INT,     &iterations },
    { '\0' }
    };

    parseArguments (argc, argv, args);

	commentator.start("block-ring test suite", "block-ring");
	bool pass = true;
	
	typedef Modular<int> Field1;  
	typedef Modular<double> Field2;  
        
        Field1 F1(q, 1);
	BlockRing<Field1> R1(F1, n);

        Field2 F2(q, 1);
        BlockRing<Field2> R2(F2, n);
                        	
	// Make sure some more detailed messages get printed
	//commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
        commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
        commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);
                

	if (!runBasicRingTests(R1, "BlockRing of Modular<int>", iterations)) pass = false;
	if (!runBasicRingTests(R2, "BlockRing of Modular<double>", iterations)) pass = false;

	commentator.stop("block-ring test suite");
	return pass ? 0 : -1;
        
}
/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s:syntax=cpp.doxygen:foldmethod=syntax
