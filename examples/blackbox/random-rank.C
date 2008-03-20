/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

// Copyright (C) 2001, 2002 Bradford Hovinen
// See COPYING for license information.
/** @name examples/blackbox/random-rank.C
 *
 * @author Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * @memo rank by several algorithms over large and small fields.
 * @doc FIXME what's it do more precisely?
 */
//@{

#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/field/givaro-gfq.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/methods.h"
#include "linbox/vector/stream.h"

using namespace LinBox;
using namespace std;

// Select our field: integers modulo a word-size (max. 31-bit) modulus
typedef Modular<LinBox::uint32> Field;
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

	commentator.setMaxDepth (2);
	commentator.setReportStream (cout);

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
