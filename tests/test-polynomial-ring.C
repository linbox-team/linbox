/* tests/test-givaropoly.C
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


/*! @file  tests/test-gf2.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"

#include <iostream>


#include "linbox/ring/polynomial-ring.h"
#include "linbox/ring/modular.h"
#include <givaro/givquotientdomain.h>

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
	typedef PolynomialRing<BaseDom> PolyDom;
	typedef PolynomialRing<PolyDom> Bivariate;
	

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (4);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	commentator().start ("Testing LB Polynomial ring", "main", 10);
	
	BaseDom Fp(p);
	PolyDom Poly(Fp,"X");
	if ( not testRing (Poly, "PolynomialRing<Modular<float>>"))
		pass = false;
	commentator().progress ();
	
	Bivariate PPo(Poly,"Y");
	if ( not testRing (PPo, "PolynomialRing<PolynomialRing<Modular<float>>>"))
		pass = false;
	commentator().progress ();
	
    PolyDom::Element Q;
    Poly.init(Q, Givaro::Degree(2)); // X^2

    typedef Givaro::QuotientDom<PolyDom> Quotient;
    Quotient QD(Poly, Q);

    if ( not testRing (QD, "QuotientDom<PolynomialRing<Modular<float>>>"))
		pass = false;
	commentator().progress ();

    Fp.init(Q[0], 1); // X^2+1
    if (p==2) Fp.init(Q[1],1); // X^2+X+1

    if ( not testField(QD, "Irreducible QuotientDom<PolynomialRing<Modular<float>>>"))
		pass = false;
	commentator().progress ();

	commentator().stop("PolynomialRing test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
