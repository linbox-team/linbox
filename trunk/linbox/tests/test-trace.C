
/* tests/test-trace.C
 * Copyright (C) -bds
 *
 * Earlier version by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------------
 *
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
 *
 */

/*! @file  tests/test-trace.C
 * @ingroup tests
 *
 * @brief no doc
 *
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>


#include "test-common.h"

#include "linbox/field/modular.h"
#include "linbox/solutions/trace.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 20;
	static integer q = 101;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

	commentator().start("Trace test suite", "Trace");
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	ostream& report = commentator().report();

	typedef Modular<int> Field;
	typedef Field::Element Element;
	// typedef vector<Element> Vector;
	Field F (q);
	Element t, t1, t2, t3;
	F.assign(t, F.zero); F.assign(t1, F.zero); F.assign(t2, F.zero); F.assign(t3, F.zero);

// scalar matrix
	Element s; F.init(s, 3);
	ScalarMatrix<Field> A(F, n, n, s);
	F.init(t, n*3);

	trace(t1, A);
	if (! F.areEqual(t1, t) ) pass = false;

	SolutionTags::Local L;
	trace(t2, A, L );
	if (! F.areEqual(t2, t) ) pass = false;

	SolutionTags::Generic G;
	trace(t3, A, G );
	if (! F.areEqual(t3, t) ) pass = false;

#ifdef __LINBOX_HAVE_NTL
// Toeplitz matrix

	NTL_zz_p CF( q );
	NTL_zz_p::Element u, u1, u2, u3;
	CF.assign(u, CF.zero); CF.assign(u1, CF.zero); CF.assign(u2, CF.zero); CF.assign(u3, CF.zero);
    NTL_zz_pX PF(CF);

    NTL_zz_p::Element temp;
    NTL_zz_pX::Element poly;
    PF.assign(poly,PF.zero);

    for( int diff = 1 - ((int)n); diff <= ((int)n) - 1; ++diff ) {
		PF.setCoeff(poly,(size_t)((size_t)diff + n - 1), CF.init(temp,diff) );
	}

	Toeplitz<NTL_zz_p,NTL_zz_pX> B( PF, poly, n );

	CF.assign(u, CF.zero);
	trace(u1, B);
	if (! CF.areEqual(u1, u)) {
		pass = false;
		report << "u and u1: " << u << " " << u1 << endl;
	}

	trace(u2, B, SolutionTags::Local() );
	if (! CF.areEqual(u2, u)) {
		pass = false;
		report << "u and u2: " << u << " " << u2 << endl;
	}

	trace(u3, B, SolutionTags::Generic() );
	if (! CF.areEqual(u3, u)) {
		pass = false;
		report << "u and u3: " << u << " " << u3 << endl;
	}
#endif // __LINBOX_HAVE_NTL

	commentator().stop("Trace solution test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

