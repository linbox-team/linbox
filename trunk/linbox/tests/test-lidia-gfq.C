/* -*- mode: c; style: linux -*- */

/* tests/test-lidia-gfq.C
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

#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <iostream>
#include <fstream>
#include <vector>


#include "linbox/field/lidia-gfq.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
        static integer q = 10733;
 

        static Argument args[] = {
                { 'q', "-q Q", "Operate over the \"field\" GF(Q) [1] (default 1073741789)", TYPE_INTEGER, &q },
                { '\0' }
        };

        parseArguments (argc, argv, args);




	cout << " \n\nLiDIAGfq field test suite" << endl << endl;
	cout.flush ();
	bool pass = true;
	
	LidiaGfq F(q,1);
	
	// Make sure some more detailed messages get printed
	commentator.getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);

	if (!testField<LidiaGfq > (F, "Testing LidiaGfq as a prime field"))
	
	  pass = false;

#if 0
	FieldArchetype K(new LidiaGfq(101,1));

	if (!testField<FieldArchetype> (K, "Testing archetype with envelope of LidiaGfq prime field"))
		pass = false;
#endif


	return pass ? 0 : -1;
}
