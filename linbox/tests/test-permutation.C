
/* tests/test-permutation.C
 * bds
 * adapted from test-diagonal by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * Time-stamp: <22 Jun 10 15:59:39 Jean-Guillaume.Dumas@imag.fr>
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


/*! @file  tests/test-diagonal.C
 * @ingroup tests
 * @brief  no doc
 * @test NO DOC
 */



#include "linbox/linbox-config.h"


#include <cstdio>

#include "linbox/blackbox/permutation.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/util/commentator.h"
#include "linbox/field/modular.h"
#include "linbox/vector/stream.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

// Test 1: PP^T should = I.
template <class Blackbox>
static bool testInvEqTrans (Blackbox & P)
{
	typedef typename Blackbox::Field Field;
	const Field & F = P.field();
	typedef BlasVector<Field> Vector;

	Transpose<Blackbox> PT(P);

	RandomDenseStream<Field, Vector> stream1 (F, P.rowdim(), 3);
	commentator().start ("Testing PP^T = I", "testInvEqTrans", stream1.m ());

	Vector u(F), v(F), w(F);
	VectorWrapper::ensureDim (u, stream1.n ());
	VectorWrapper::ensureDim (v, stream1.n ());
	VectorWrapper::ensureDim (w, stream1.n ());
	stream1.next (u);

	PT.apply(v, u);
	P.apply(w, v);

	VectorDomain<Field> VD (F);
	bool pass = VD.areEqual(w, u);
	if (! pass || (pass && P.rowdim() <= 20)) {
		ostream &report = commentator().report (Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
		VD.write (report << "u = ", u) << endl;
		//PT.write(report << "PT = ") << endl; //transpose has no write!
		VD.write (report << "v = P^Tu = ", v) << endl;
		P.write(report << "P = ") << endl;
		VD.write (report << "w = Pv = ", w) << endl;
	}

	commentator().stop (MSG_STATUS (pass), (const char *) 0, "testInvEqTrans");

	return pass;
}

int main (int argc, char **argv)
{
	bool pass = true;

	static size_t n = 10;
	static integer q = 2147483647U;

	static Argument args[] = {
		{ 'n', "-n N", "Set dimension of test matrices to NxN.", TYPE_INT,     &n },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	typedef Modular<uint32_t> Field;

	parseArguments (argc, argv, args);
	Field F (q);

	srand ((unsigned)time (NULL));

	commentator().start("Permutation matrix black box test suite", "permutation");

	// Make sure some more detailed messages get printed
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDepth (3);
	commentator().getMessageClass (INTERNAL_DESCRIPTION).setMaxDetailLevel (Commentator::LEVEL_UNIMPORTANT);

	LinBox::Permutation<Field> P(F);
	P.random(n);
	pass = pass && testBlackboxNoRW(P);
	pass = pass && testInvEqTrans(P);

	commentator().stop (MSG_STATUS (pass));

	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
