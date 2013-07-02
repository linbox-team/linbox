/* Copyright (C) LinBox
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.	 See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/*! @file  tests/test-triplesbb-omp.C
 * @ingroup tests
 *
 * @brief no doc
 *
 * @test no doc.
 */


#include "linbox/linbox-config.h"

#include <iostream>
#include <fstream>

#include <omp.h>
#include <set>
#include <utility>

#include "linbox/field/modular.h"
#include "linbox/blackbox/triplesbb-omp.h"
#include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/blas-matrix.h"

#include "test-common.h"
#include "test-generic.h"

using namespace LinBox;

int randRange(int start, int end)
{
        double rval = rand();
        static const double NORMALIZING_CONSTANT = 1.0/(1.0+RAND_MAX);
        double normedRVal = rval*NORMALIZING_CONSTANT;
        double rangeSize = end-start;
        int offset = rangeSize*normedRVal;
        return start+offset;
} 

int main (int argc, char **argv)
{

	bool pass = true;

	static size_t m = 400;
	static size_t n = 200;
        static size_t nnz = m*(n/100);
	static integer q = 101;

	static Argument args[] = {
		{ 'm', "-m M", "Set row dimension of test matrix to M.", TYPE_INT,     &m },
		{ 'n', "-n N", "Set col dimension of test matrix to N.", TYPE_INT,     &n },
		{ 'z', "-z NNZ", "Set number of nonzero entries in test matrix.", TYPE_INT,     &nnz },
		{ 'q', "-q Q", "Operate over the \"field\" GF(Q) [1].", TYPE_INTEGER, &q },
		END_OF_ARGUMENTS
	};

	parseArguments (argc, argv, args);

        if (nnz > m*n) {
                nnz=m*(n/100);
        }

	srand ((unsigned)time (NULL));

	commentator().start("TriplesBBOMP black box test suite", "triplesbbomp");

	typedef Modular<double> Field;
	typedef Field::Element Element;
	typedef vector <Element> Vector;
	typedef TriplesBBOMP<Field> OMPBlackbox;
        typedef TriplesBB<Field> SeqBlackbox;
        typedef BlasMatrix<Field> DenseMat;

	Field F (q);

        //DenseMat MatIn(F,m,n);
	OMPBlackbox TestMat(F, m, n);
        SeqBlackbox RefMat(F, m, n);
        Vector x(n), testY(m), refY(m);
	Element d;

        typedef pair<size_t,size_t> CoordPair;
        typedef set<CoordPair> PairSet;
        PairSet pairs;
	for(int i = 0; i < (int)nnz; ++i)
	{
                size_t row,col;
                do {
                        row = randRange(0,m);
                        col = randRange(0,n);
                } while (pairs.count(CoordPair(row,col))!=0);

                F.init(d, randRange(0,q));
                TestMat.setEntry(row,col,d);
                RefMat.setEntry(row,col,d);
                pairs.insert(CoordPair(row,col));
	}
        TestMat.sortBlock();
        TestMat.sortRow();
        for (int i=0;i<(int)n;++i) {
                F.init(x[i],randRange(0,q));
        }
        omp_set_num_threads(2);
	ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);
        report << "Testing triplesbb-omp against triplesbb" << std::endl;
        TestMat.apply(testY,x);
        RefMat.apply(refY,x);
        if (testY==refY) {
                report << "PASS: Equality test passed" << std::endl;
        } else {
                pass=false;
                report << "FAILURE: Equality test failed" << std::endl;
        }
	commentator().stop("TriplesBBOMP black box test suite");
	return pass ? 0 : -1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s

