/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* tests/test-block-ring.C
 * Copyright (C) 2007 b d saunders
 *
 */

#include "linbox-config.h"

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
    static integer q = 10733;
	static size_t n = 10;
	static int iterations = 2;

    static Argument args[] = {
            { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 10733)", TYPE_INTEGER, &q },
	{ 'n', "-n N", "Set dimension of blocks to N (default 10)", TYPE_INT,     &n },
	{ 'i', "-i I", "Perform each test for I iterations (default 2)",      TYPE_INT,     &iterations },
            { '\0' }
    };

    parseArguments (argc, argv, args);

	std::cout << endl << "block-ring test suite" << std::endl;
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
	return pass ? 0 : -1;
        
}
