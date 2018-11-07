/* tests/test-poweroftwomodular.C
 * Copyright (C) 2014 The LinBox group,
 *
 * Written by Gavin Harrison <gmh33@drexel.edu>, Jean-Guillaume.Dumas@imag.fr
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


/*! @file  tests/test-poweroftwomodular.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>


#include "linbox/ring/poweroftwomodular.h"

#include "test-field.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static int p = 13;
	static int e = 1; // currently not used

	static Argument args[] = {
		{ 'p', "-p P", "Set the base field prime.", TYPE_INT, &p },
		{ 'e', "-e E", "Set the base field exponent.", TYPE_INT, &e },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("PowerOfTwoGivaro::Modular test", "2powRing");
	bool pass = true;

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	PowerOfTwoGivaro::Modular R();
	if ( not testRing (R, "PowerOfTwoGivaro::Modular"))
		pass = false;
	
	commentator().start("PowerOfTwoGivaro::Modular test");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
