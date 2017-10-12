
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


//#include "linbox/ring/polynomial-ring.h"
#include <givaro/givpoly1factor.h>
#include "linbox/ring/modular.h"

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

	commentator().start("PolynomialRing test suite", "PolynomialRing");
	bool pass = true;

	typedef Givaro::Modular<float> BaseDom;
	typedef Givaro::Poly1Dom<BaseDom, Givaro::Dense> PolyDom;
	typedef Givaro::Poly1FactorDom<BaseDom, Givaro::Dense> PolyFactorDom;
    // PolynomialRing has some ring members missing, eg., init, one, mOne, zero.
	//typedef PolynomialRing<BaseDom> PolyDom;
	
	BaseDom Fp(p);
	PolyDom Poly(Fp);
	PolyFactorDom PolyFac(Fp);

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator().start ("Testing GivaroPoly", "main", 10);
	
	if ( not testRing (Poly, "Poly1Dom<Modular<float>>"))
		pass = false;
	commentator().progress ();
	
	if ( not runPIRTests(Poly, "Poly1Dom<Modular<float>>"))
		pass = false;
	commentator().progress ();
	
	if ( not testRing (PolyFac, "Poly1FactorDom<Modular<float>>"))
		pass = false;
	commentator().progress ();
	
	if ( not runPIRTests(PolyFac, "Poly1FactorDom<Modular<float>>"))
		pass = false;
	commentator().progress ();
	
	commentator().stop("PolynomialRing test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

