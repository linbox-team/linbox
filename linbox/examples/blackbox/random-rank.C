
/*
 * examples/blackbox/random-rank.C
 *
 * Copyright (C) 2001, 2002, 2010 Bradford Hovinen <hovinen@cis.udel.edu>
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
 */
/** @name examples/blackbox/random-rank.C
 *
 * @author Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * @memo rank by several algorithms over large and small fields.
 * @doc FIXME what's it do more precisely?
 */
//@{

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/field/givaro-gfq.h"
#include "linbox/matrix/sparse.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/methods.h"
#include "linbox/vector/stream.h"

using namespace LinBox;
using namespace std;

// Select our field: integers modulo a word-size (max. 31-bit) modulus
typedef Modular<LinBox::uint32_t> Field;
typedef GivaroGfq FieldExtn;

// The Gaussian elimiation code only works with sparse sequence vectors, so we'll use those
typedef Vector<Field>::SparseSeq Row;
typedef Vector<FieldExtn>::SparseSeq RowE;

// Select our black box: a sparse matrix over the above-mentioned field with
// default application vector and row representation types
typedef SparseMatrix<Field, Row> Blackbox;
typedef SparseMatrix<FieldExtn, RowE> BlackboxE;

// Constants: we are working with an n x n matrix over GF(q)
const int n = 1000;
//const int q = 65521U;
const int q = 2U;
const double p = 0.026;

/// no command line args
int main (int argc, char **argv)
{
	unsigned long rank_A_Wiedemann, rank_A_elimination;

	commentator().setMaxDepth (2);
	commentator().setReportStream (cout);

	Field F (q);
	FieldExtn E (q, 10);

	RandomSparseStream<Field, Row> A_stream (F, p, n, n);
	RandomSparseStream<FieldExtn, RowE> B_stream (E, p, n, n);

	Blackbox A (F, A_stream);
	BlackboxE B (E, B_stream);

	rank (rank_A_Wiedemann, A, F, MethodTrait::Wiedemann ());
	rank (rank_A_elimination, A, F, MethodTrait::Elimination ());

	cout << "Rank by Wiedemann is " << rank_A_Wiedemann << endl;
	cout << "Rank by Elimination is " << rank_A_elimination << endl;

	rank (rank_A_Wiedemann, B, E, MethodTrait::Wiedemann ());
	rank (rank_A_elimination, B, E, MethodTrait::Elimination ());

	cout << "Over Extension Field Rank by Wiedemann is " << rank_A_Wiedemann << endl;
	cout << "Over Extension Field Rank by Elimination is " << rank_A_elimination << endl;

	return 0;
}
//@}

// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,:0,t0,+0,=s
// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:

