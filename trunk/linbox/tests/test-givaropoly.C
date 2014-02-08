
/* tests/test-givaropoly.C
 * Copyright (C) 2014 Gavin Harrison,
 *
 * Written by Gavin Harrison <gmh33@drexel.edu>,
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


/*! @file  tests/test-gf2.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <sstream>

#include <queue>

#include "linbox/ring/givaro-poly.h"
#include "givaro/givgfq.h"

#include "test-field.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	static int p = 13;
	static int e = 1;

	static Argument args[] = {
		{ 'p', "-p P", "Set the base field prime.", TYPE_INT, &p },
		{ 'e', "-e E", "Set the base field exponent.", TYPE_INT, &e },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("GivaroPoly field test suite", "GivaroPoly");
	bool pass = true;

	typedef Givaro::GFqDom<int> BaseDom;
	typedef typename Givaro::Poly1Dom<BaseDom, Givaro::Dense> PolyDom;
	
	BaseDom GFq(p, e);
	PolyDom Poly(GFq);
	GivaroPoly<PolyDom> F(Poly);

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator().start ("Testing GivaroPoly", "main", 10);
	
	if ( not testField (F, "GivaroPoly"))
		pass = false;
	commentator().progress ();
	
	commentator().stop("GivaroPoly field test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

