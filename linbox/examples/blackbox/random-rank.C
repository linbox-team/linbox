/* -*- mode: C++; tab-width: 8; indent-tabs-mode: t; c-basic-offset: 8 -*- */

/* examples/blackbox/random-rank.C
 * Copyright (C) 2001, 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <bghovinen@math.uwaterloo.ca>
 *
 * --------------------------------------------------
 *
 * See COPYING for license information.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/blackbox/sparse.h"
#include "linbox/solutions/rank.h"
#include "linbox/solutions/methods.h"
#include "linbox/vector/stream.h"

using namespace LinBox;
using namespace std;

// Select our field: integers modulo a word-size (max. 31-bit) modulus
typedef Modular<uint32> Field;

// The Gaussian elimiation code only works with sparse sequence vectors, so we'll use those
typedef Vector<Field>::SparseSeq Row;

// Select our black box: a sparse matrix over the above-mentioned field with
// default application vector and row representation types
typedef SparseMatrix<Field, Vector<Field>::Dense, Row> Blackbox;

// Constants: we are working with an n x n matrix over GF(q)
const int n = 1000;
const int q = 65521U;
const double p = 0.005;

int main (int argc, char **argv)
{
	unsigned long rank_A_Wiedemann, rank_A_elimination;

	commentator.setMaxDepth (2);
	commentator.setReportStream (cout);

	Field F (q);

	RandomSparseStream<Field, Row> A_stream (F, p, n, n);

	Blackbox A (F, A_stream);

	rank (rank_A_Wiedemann, A, F, MethodTrait::Wiedemann ());
	rank (rank_A_elimination, A, F, MethodTrait::Elimination ());

	cout << "Rank by Wiedemann is " << rank_A_Wiedemann << endl;
	cout << "Rank by Elimination is " << rank_A_elimination << endl;

	return 0;
}
