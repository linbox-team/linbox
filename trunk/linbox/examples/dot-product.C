/* -*- mode: C++; style: linux -*- */

/* examples/dot-product.C
 * Copyright (C) 2002 Bradford Hovinen
 *
 * Written by Bradford Hovinen <hovinen@cis.udel.edu>
 *
 * --------------------------------------------------
 *
 * See COPYING for license information.
 *
 * --------------------------------------------------
 * Constructs random vectors and computes their dot product, giving the
 * required time.
 */

#include "linbox-config.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "linbox/field/modular.h"
#include "linbox/field/vector-domain.h"
#include "linbox/util/vector-factory.h"
#include "linbox/util/commentator.h"

using namespace LinBox;
using namespace std;

typedef Modular<long> Field;
typedef vector<Field::Element> DenseVector;
typedef vector<pair<size_t, Field::Element> > SparseSeqVector;
typedef pair<vector<size_t>, vector<Field::Element> > SparseParVector;

// Constants: we are working with an n x n matrix over GF(q)
const int n = 10000000;
const int k = 10000;
const int q = 32749;

int main (int argc, char **argv)
{
	srand (time (NULL));

	commentator.setMaxDepth (2);
	commentator.setReportStream (cout);

	Field F (q);

	RandomDenseVectorFactory<Field> factory1 (F, n);
	RandomSparseSeqVectorFactory<Field> factory2 (F, n, k);
	RandomSparseParVectorFactory<Field> factory3 (F, n, k);

	DenseVector v1, v2;
	SparseSeqVector v3;
	SparseParVector v4;

	factory1.next (v1);
	factory1.next (v2);
	factory2.next (v3);
	factory3.next (v4);

	VectorDomain<Field> VD (F);
	Field::Element res;

	commentator.start ("dense/dense dot product (1)");
	for (int i = 0; i < 1; i++)
		VD.dot (res, v1, v2);
	commentator.stop ("done");

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
