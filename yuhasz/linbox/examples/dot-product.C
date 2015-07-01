/* -*- mode: C++; style: linux -*- */

// Copyright (C) 2002 Bradford Hovinen
// See COPYING for license information.
/** @name examples/dot-product.C
 *
 * @author Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * @memo  
 * Use of vectors meeting the LinBox dense and sparse vector archetypes
 * is illustrated and their dot-product times are benchmarked.
 *
 * @doc
 * Constructs random vectors and computes their dot product, giving the
 * required time.
 */
//@{
#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/vector/stream.h"
#include "linbox/util/commentator.h"

using namespace LinBox;
using namespace std;

typedef Modular<uint32> Field;

// Constants: we are working with an n x n matrix over GF(q)
const int n = 10000000;
const double p = .001;
const int q = 32749;

/// no command line args
int main (int argc, char **argv)
{
	commentator.setMaxDepth (1);
	commentator.setReportStream (cout);

	Field F (q);

	RandomDenseStream<Field, Vector<Field>::Dense> factory1 (F, n);
	RandomSparseStream<Field, Vector<Field>::SparseSeq> factory2 (F, p, n);
	RandomSparseStream<Field, Vector<Field>::SparsePar> factory3 (F, p, n);

	Vector<Field>::Dense v1, v2;
	Vector<Field>::SparseSeq v3;
	Vector<Field>::SparsePar v4;

	factory1 >> v1 >> v2;
	factory2 >> v3;
	factory3 >> v4;

	VectorDomain<Field> VD (F);
	Field::Element res;

	commentator.start ("dense/dense dot product (1)");
	for (int i = 0; i < 1; i++)
		VD.dot (res, v1, v2);
	commentator.stop ("done");

	//////  FIXME  sparse seq dot segfaults on stimpy -bds
	commentator.start ("dense/sparse sequence dot product (1000)");
	for (int i = 0; i < 1000; i++)
		VD.dot (res, v1, v3);
	commentator.stop ("done");

	commentator.start ("dense/sparse parallel dot product (1000)");
	for (int i = 0; i < 1000; i++)
		VD.dot (res, v1, v4);
	commentator.stop ("done");

	return 0;
}
//@}
