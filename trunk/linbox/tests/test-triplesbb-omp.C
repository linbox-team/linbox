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

typedef Modular<double> Field;
typedef Field::Element Element;
typedef std::vector <Element> STLVector;
typedef TriplesBBOMP<MatrixDomain<Field> > OMPBlackbox;
typedef TriplesBB<MatrixDomain<Field> > SeqBlackbox;
typedef BlasMatrix<Field> DenseMat;


bool testRand(size_t m, size_t n, size_t nnz, integer q, bool shouldFail, int numThreads, ostream &report)
{
        bool pass=true;

	Field F (q);

	OMPBlackbox TestMat(F, m, n);
        SeqBlackbox RefMat(F, m, n);
        STLVector x(n), xT(m), testY(m), refY(m), testYT(n), refYT(n);
	Element d;

        for (int i=0;i<(int)n;++i) {
                F.init(x[i],randRange(0,q));
        }

        for (int i=0;i<(int)m;++i) {
                F.init(xT[i],randRange(0,q));
        }

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

                int randVal=randRange(0,q);
                if (shouldFail && (i==((int)nnz/2))) {
                        if (F.isZero(x[col])) {
                                int temp=randRange(1,q);
                                F.init(x[col],temp);
                        }
                        if (F.isZero(xT[row])) {
                                int temp=randRange(1,q);
                                F.init(xT[row],temp);
                        }
                        F.init(d, randVal);
                        RefMat.setEntry(row,col,d);

                        randVal = (randVal+1)%q;
                        F.init(d, randVal);
                        TestMat.setEntry(row,col,d);

                } else {
                        F.init(d, randRange(0,q));
                        TestMat.setEntry(row,col,d);
                        RefMat.setEntry(row,col,d);
                }
                pairs.insert(CoordPair(row,col));
	}

        TestMat.sortBlock();
        TestMat.sortRow();

        omp_set_num_threads(numThreads);

        if (shouldFail) {
                report << "Testing triplesbb-omp apply() against triplesbb apply() on distinct random matrices using: " <<
                        "m=" << m << " n=" << n << " nnz=" << nnz << " numthreads=" << numThreads << " q=" << q << std::endl;
                TestMat.apply(testY,x);
                RefMat.apply(refY,x);
                if (testY!=refY) {
                        report << "PASS: Apply-Inequality test passed" << std::endl;
                } else {
                        pass=false;
                        report << "FAILURE: Apply-Inequality test failed" << std::endl;
                }

                report << "Testing triplesbb-omp applyTranspose() against triplesbb applyTranspose() on distinct random matrices using: " <<
                        "m=" << m << " n=" << n << " nnz=" << nnz << " numthreads=" << numThreads << " q=" << q << std::endl;
                TestMat.applyTranspose(testYT,xT);
                RefMat.applyTranspose(refYT,xT);
                if (testYT!=refYT) {
                        report << "PASS: Transpose-Inequality test passed" << std::endl;
                } else {
                        pass=false;
                        report << "FAILURE: Transpose-Inequality test failed" << std::endl;
                }
        } else {
                report << "Testing triplesbb-omp apply() against triplesbb apply() on random matrix using: " <<
                        "m=" << m << " n=" << n << " nnz=" << nnz << " numthreads=" << numThreads << " q=" << q << std::endl;
                TestMat.apply(testY,x);
                RefMat.apply(refY,x);
                if (testY==refY) {
                        report << "PASS: Apply-Equality test passed" << std::endl;
                } else {
                        pass=false;
                        report << "FAILURE: Apply-Equality test failed" << std::endl;
                }

                report << "Testing triplesbb-omp applyTranspose() against triplesbb applyTranspose() on random matrix using: " <<
                        "m=" << m << " n=" << n << " nnz=" << nnz << " numthreads=" << numThreads << " q=" << q << std::endl;
                TestMat.applyTranspose(testYT,x);
                RefMat.applyTranspose(refYT,x);
                if (testYT==refYT) {
                        report << "PASS: Transpose-Equality test passed" << std::endl;
                } else {
                        pass=false;
                        report << "FAILURE: Transpose-Equality test failed" << std::endl;
                }
        }
        return pass;
}

bool testRandSuite(integer q, ostream &report)
{
        bool pass=true;

        pass = pass && testRand(100,100,100,q,false,1,report);
        pass = pass && testRand(100,100,100,q,true,1,report);
        pass = pass && testRand(100,100,100,q,false,2,report);
        pass = pass && testRand(100,100,100,q,true,2,report);
        pass = pass && testRand(100,100,100,q,false,4,report);
        pass = pass && testRand(100,100,100,q,true,4,report);

        pass = pass && testRand(100000,100000,100000,q,false,1,report);
        pass = pass && testRand(100000,100000,100000,q,true,1,report);
        pass = pass && testRand(100000,100000,100000,q,false,2,report);
        pass = pass && testRand(100000,100000,100000,q,true,2,report);
        pass = pass && testRand(100000,100000,100000,q,false,4,report);
        pass = pass && testRand(100000,100000,100000,q,true,4,report);

        pass = pass && testRand(100000,1000,100000,q,false,1,report);
        pass = pass && testRand(100000,1000,100000,q,true,1,report);
        pass = pass && testRand(100000,1000,100000,q,false,2,report);
        pass = pass && testRand(100000,1000,100000,q,true,2,report);
        pass = pass && testRand(100000,1000,100000,q,false,4,report);
        pass = pass && testRand(100000,1000,100000,q,true,4,report);

        pass = pass && testRand(1000,100000,100000,q,false,1,report);
        pass = pass && testRand(1000,100000,100000,q,true,1,report);
        pass = pass && testRand(1000,100000,100000,q,false,2,report);
        pass = pass && testRand(1000,100000,100000,q,true,2,report);
        pass = pass && testRand(1000,100000,100000,q,false,4,report);
        pass = pass && testRand(1000,100000,100000,q,true,4,report);

        pass = pass && testRand(10000,10000,1000000,q,false,1,report);
        pass = pass && testRand(10000,10000,1000000,q,true,1,report);
        pass = pass && testRand(10000,10000,1000000,q,false,2,report);
        pass = pass && testRand(10000,10000,1000000,q,true,2,report);
        pass = pass && testRand(10000,10000,1000000,q,false,4,report);
        pass = pass && testRand(10000,10000,1000000,q,true,4,report);

        pass = pass && testRand(100000,100000,1000000,q,false,2,report);
        pass = pass && testRand(100000,100000,1000000,q,true,2,report);

        return pass;
}

bool isScaled(STLVector y,STLVector x,int n,int scaleFactor)
{
        for (int i=0;i<n;++i) {
                if (y[i]!=(x[i]*scaleFactor)) {
                        return false;
                }
        }
        return true;
}

/*
bool testDenseRand(size_t n, size_t m, integer q, bool shouldFail, int numThreads, ostream &report)
{
        bool pass=true;

	Field F (q);

        DenseMat MatIn(F,m,n);
        DenseMat TestY(F,m,n);
        DenseMat RefY(F,m,n);
	OMPBlackbox TestMat(F, n, n);
        SeqBlackbox RefMat(F,n,n);
        Element d;

        TestY.zero();
        RefY.zero();
        for (int i=0;i<n;++i) {
                for (int j=0;j<m;++j) {
                        F.init(d,randRange(0,q));
                        MatIn.setEntry(i,j,d);
                }
        }

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

                int randVal=randRange(0,q);
                F.init(d, randRange(0,q));
                TestMat.setEntry(row,col,d);
                RefMat.setEntry(row,col,d);
                pairs.insert(CoordPair(row,col));
        }

        TestMat.sortBlock();
        TestMat.sortRow();

        omp_set_num_threads(numThreads);

        report << "Testing triplesbb-omp apply_right against triplesbb apply_right on random matrix using: " <<
                        "m=" << m << " n=" << n << " nnz=" << nnz << " numthreads=" << numThreads << " q=" << q << std::endl;
        TestMat.apply(TestY,MatIn);
        RefMat.apply(RefY,MatIn);

        if (testY==refY) {
                report << "PASS: apply_right test passed" << std::endl;
        } else {
                pass=false;
                report << "FAILURE: apply_right test failed" << std::endl;
        }

        return pass;
}
*/

bool testIdent(size_t n, integer q, bool shouldFail, int numThreads, ostream &report)
{
        bool pass=true;

	Field F (q);

	OMPBlackbox TestMat(F, n, n);
        STLVector x(n), y(n);
	Element d;

        for (int i=0;i<(int)n;++i) {
                F.init(x[i],randRange(1,q/2));
        }

        typedef pair<size_t,size_t> CoordPair;
        typedef set<CoordPair> PairSet;
        PairSet pairs;
        int scaleFactor=(q==2)?1:2;
	for(int i = 0; i < (int)n; ++i)
	{
                if (shouldFail && (i==(n/3))) {
                        F.init(d,0);
                } else {
                        F.init(d,scaleFactor);
                }
                TestMat.setEntry(i,i,d);
        }

        TestMat.sortBlock();
        TestMat.sortRow();

        omp_set_num_threads(numThreads);
        if (shouldFail) {
                report << "Testing triplesbb-omp apply() on non-diagonal matrix using: "<<
                        "n=" << n << " q=" << q << " numThreads=" << numThreads << std::endl;
                TestMat.apply(y,x);
                if (isScaled(y,x,n,scaleFactor)) {
                        report << "FAILURE: Diagonal-nonscale test failed" << std::endl;
                        pass=false;
                } else {
                        report << "PASS: Diagonal-nonscale test passed" << std::endl;
                }
        } else {
                report << "Testing triplesbb-omp apply() on diagonal matrix using: " <<
                        "n=" << n << " q=" << q << " numThreads=" << numThreads << std::endl;
                TestMat.apply(y,x);
                if (!isScaled(y,x,n,scaleFactor)) {
                        report << "FAILURE: Diagonal-scale test failed" << std::endl;
                        pass=false;
                } else {
                        report << "PASS: Diagonal-scale test passed" << std::endl;
                }
        }

        if (shouldFail) {
                report << "Testing triplesbb-omp applyTranspose() on non-diagonal matrix using: "<<
                        "n=" << n << " q=" << q << " numThreads=" << numThreads << std::endl;
                TestMat.applyTranspose(y,x);
                if (isScaled(y,x,n,scaleFactor)) {
                        report << "FAILURE: Diagonal-nonscale test failed" << std::endl;
                        pass=false;
                } else {
                        report << "PASS: Diagonal-nonscale test passed" << std::endl;
                }
        } else {
                report << "Testing triplesbb-omp applyTranspose() on diagonal matrix using: " <<
                        "n=" << n << " q=" << q << " numThreads=" << numThreads << std::endl;
                TestMat.applyTranspose(y,x);
                if (!isScaled(y,x,n,scaleFactor)) {
                        report << "FAILURE: Diagonal-scale test failed" << std::endl;
                        pass=false;
                } else {
                        report << "PASS: Diagonal-scale test passed" << std::endl;
                }
        }


        return pass;
}

bool testIdentSuite(integer q,ostream &report)
{
        bool pass=true;

        pass = pass && testIdent(100,q,true,1,report);
        pass = pass && testIdent(100,q,false,1,report);
        pass = pass && testIdent(100,q,true,2,report);
        pass = pass && testIdent(100,q,false,2,report);
        pass = pass && testIdent(100,q,true,4,report);
        pass = pass && testIdent(100,q,false,4,report);

        pass = pass && testIdent(100000,q,false,1,report);
        pass = pass && testIdent(100000,q,true,1,report);
        pass = pass && testIdent(100000,q,false,1,report);
        pass = pass && testIdent(100000,q,true,2,report);
        pass = pass && testIdent(100000,q,false,2,report);
        pass = pass && testIdent(100000,q,true,4,report);
        pass = pass && testIdent(100000,q,false,4,report);

        return pass;
}

int main (int argc, char **argv)
{

	bool pass = true;

	static size_t m = 400;
	static size_t n = 200;
        static size_t nnz = m*(n/100);
	static integer q = 101;

        // TODO: Remove this, it's no longer used
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

	ostream &report = LinBox::commentator().report (LinBox::Commentator::LEVEL_IMPORTANT, INTERNAL_DESCRIPTION);

        pass = pass && testIdentSuite(2,report);
        pass = pass && testIdentSuite(2147483629,report);

        pass = pass && testRandSuite(2,report);
        pass = pass && testRandSuite(2147483629,report);

        report << "Testing" << std::endl;

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

