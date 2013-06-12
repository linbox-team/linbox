/* tests/test-direct-sum.C
 * Copyright (C) LinBox
 * Written by Austin Lobo, David Saunders
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
 *.
 */

/*! @file  tests/test-frobenius.C
 * @ingroup tests
 * @brief no doc.
 * @test NO DOC
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>


#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/scalar-matrix.h"
#include "linbox/blackbox/frobenius.h"

#include "test-generic.h"

using namespace LinBox;

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 101;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	typedef Modular<uint32_t> Field;
	typedef vector<Field::Element> Vector;
	typedef Vector Polynomial;
	typedef vector<Polynomial> Plist;

	parseArguments (argc, argv, args);
	Field F (q);
	Plist plist(3);

	Field::RandIter r(F);

	size_t  pdeg = 10;
	plist[0].resize(pdeg+1);
	for ( size_t i=0; i < pdeg; ++i) r.random(plist[0][i]);
	F.init(plist[0][pdeg],1);

	pdeg = 6;
	plist[1].resize(pdeg+1);
	for ( size_t i=0; i < pdeg; ++i) r.random(plist[1][i]);
	F.init(plist[1][pdeg],1);

	pdeg = 4;
	plist[2].resize(pdeg+1);
	for ( size_t i=0; i < pdeg; ++i) r.random(plist[2][i]);
	F.init(plist[2][pdeg],1);

	commentator().start("Frobenius form black box test suite", "frobenius");
	Frobenius<Field>  A(F, plist.begin(), plist.end());

	pass = pass && testBlackboxNoRW(A);

	commentator().stop("Frobenius form black box test suite");
	return pass ? 0 : -1;
}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

