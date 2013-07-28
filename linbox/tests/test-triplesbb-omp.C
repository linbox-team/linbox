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
#include <sstream>

#include "linbox/field/modular.h"
#include "linbox/blackbox/triplesbb-omp.h"
#include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/blas-matrix.h"

#include "examples/map-sparse.h"

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
typedef TriplesBBOMP<Field> OMPBlackbox;
typedef TriplesBB<Field> SeqBlackbox;

void generateRandMat(MapSparse<Field>& mat,int m, int n, int nnz, int q)
{
	Element d;

        typedef pair<size_t,size_t> CoordPair;
        typedef set<CoordPair> PairSet;
        PairSet pairs;

	for(int i = 0; i < (int)nnz; ++i) {
                size_t row,col;
                do {
                        row = randRange(0,m);
                        col = randRange(0,n);
                } while (pairs.count(CoordPair(row,col))!=0);
                
                mat.field().init(d, randRange(0,q));
                mat.setEntry(row,col,d);
                pairs.insert(CoordPair(row,col));
        }
}

void generateRandVec(STLVector& vec,int n, int q)
{
        Field F(q);

        for (int i=0;i<n;++i) {
                F.init(vec[i],randRange(0,q));
        }
}

void varyMatVecPair(MapSparse<Field>& leftMat,STLVector& rightVec, int m, int n, int q)
{
        Field F(q);
        Element d;

        int row=randRange(0,m);
        int col=randRange(0,n);

        if (F.isZero(rightVec[col])) {
                F.init(rightVec[col],randRange(1,q));
                F.init(d,randRange(1,q));
        } else {
                do {
                        F.init(d,randRange(0,q));
                } while (F.areEqual(d,leftMat.getEntry(row,col)));
        }
        leftMat.setEntry(row,col,d);
}

void varyMatVecPair(STLVector& leftVec, MapSparse<Field>& rightMat, int m, int n, int q)
{
        Field F(q);
        Element d;

        int row=randRange(0,m);
        int col=randRange(0,n);

        if (F.isZero(rightMat.getEntry(row,col))) {
                F.init(d,randRange(1,q));
                rightMat.setEntry(row,col,d);
                F.init(leftVec[row],randRange(1,q));
        } else {
                do {
                        F.init(d,randRange(0,q));
                } while (F.areEqual(d,leftVec[row]));
                leftVec[row]=d;
        }
}

void varyMatMatPair(MapSparse<Field>& leftMat, MapSparse<Field>& rightMat, int m, int n, int p, int q)
{
        Field F(q);
        Element d;

        int leftRow=randRange(0,m);
        int leftCol=randRange(0,n);
        int rightCol=randRange(0,p);

        if (F.isZero(rightMat.getEntry(leftCol,rightCol))) {
                F.init(d,randRange(1,q));
                rightMat.setEntry(leftCol,rightCol,d);
                F.init(d,randRange(1,q));
                leftMat.setEntry(leftRow,leftCol,d);
        } else {
                do {
                        F.init(d,randRange(0,q));
                } while (F.areEqual(d,leftMat.getEntry(leftRow,leftCol)));
                leftMat.setEntry(leftRow,leftCol,d);

        }
}

bool outputResult(std::string testName,
                  bool expected,
                  bool actual,
                  int m, int n, int nnz, int numThreads, int q,
                  ostream& report)
{
        report << "Testing " << testName << " using: " <<
                "m=" << m << " n=" << n << " nnz=" << nnz << " numthreads=" << numThreads << " q=" << q << std::endl;

        if (expected==actual) {
                report << "PASS: ";
        } else {
                report << "FAILURE: ";
        }
        report << testName << " expected " << expected << " got " << actual << std::endl;

        return expected==actual;
}

bool testRandVec(size_t m, size_t n, size_t nnz, int q, bool shouldFail, int numThreads, ostream &report)
{
        bool pass=true;

	Field F (q);

	OMPBlackbox TestMat(F, m, n);
        SeqBlackbox RefMat(F, m, n);
        MapSparse<Field> CommonMat(F,m,n);
        STLVector x(n), xT(m), testY(m), refY(m), testYT(n), refYT(n);
        VectorDomain<Field> VD(F);

        generateRandMat(CommonMat,m,n,nnz,q);
        generateRandVec(x,n,q);
        generateRandVec(xT,m,q);

        std::stringstream testMatSS,refMatSS;

        CommonMat.write(testMatSS);
        TestMat.read(testMatSS);
        TestMat.finalize();

        if (shouldFail) {
                varyMatVecPair(CommonMat,x,m,n,q);
                varyMatVecPair(xT,CommonMat,m,n,q);
        }

        CommonMat.write(refMatSS);
        RefMat.read(refMatSS);

        omp_set_num_threads(numThreads);

        TestMat.apply(testY,x);
        RefMat.apply(refY,x);
        pass=pass&&outputResult("apply() random",!shouldFail,VD.areEqual(testY,refY),m,n,nnz,numThreads,q,report);
        
        TestMat.applyTranspose(testYT,xT);
        RefMat.applyTranspose(refYT,xT);
        pass=pass&&outputResult("applyTranspose() random",!shouldFail,VD.areEqual(testYT,refYT),m,n,nnz,numThreads,q,report);

        return pass;
}

bool testRandSuite(int q, ostream &report)
{
        bool pass=true;

        pass = pass && testRandVec(100,100,100,q,false,1,report);
        pass = pass && testRandVec(100,100,100,q,true,1,report);
        pass = pass && testRandVec(100,100,100,q,false,2,report);
        pass = pass && testRandVec(100,100,100,q,true,2,report);
        pass = pass && testRandVec(100,100,100,q,false,4,report);
        pass = pass && testRandVec(100,100,100,q,true,4,report);

        pass = pass && testRandVec(100000,100000,100000,q,false,1,report);
        pass = pass && testRandVec(100000,100000,100000,q,true,1,report);
        pass = pass && testRandVec(100000,100000,100000,q,false,2,report);
        pass = pass && testRandVec(100000,100000,100000,q,true,2,report);
        pass = pass && testRandVec(100000,100000,100000,q,false,4,report);
        pass = pass && testRandVec(100000,100000,100000,q,true,4,report);

        pass = pass && testRandVec(100000,1000,100000,q,false,1,report);
        pass = pass && testRandVec(100000,1000,100000,q,true,1,report);
        pass = pass && testRandVec(100000,1000,100000,q,false,2,report);
        pass = pass && testRandVec(100000,1000,100000,q,true,2,report);
        pass = pass && testRandVec(100000,1000,100000,q,false,4,report);
        pass = pass && testRandVec(100000,1000,100000,q,true,4,report);

        pass = pass && testRandVec(1000,100000,100000,q,false,1,report);
        pass = pass && testRandVec(1000,100000,100000,q,true,1,report);
        pass = pass && testRandVec(1000,100000,100000,q,false,2,report);
        pass = pass && testRandVec(1000,100000,100000,q,true,2,report);
        pass = pass && testRandVec(1000,100000,100000,q,false,4,report);
        pass = pass && testRandVec(1000,100000,100000,q,true,4,report);

        pass = pass && testRandVec(10000,10000,1000000,q,false,1,report);
        pass = pass && testRandVec(10000,10000,1000000,q,true,1,report);
        pass = pass && testRandVec(10000,10000,1000000,q,false,2,report);
        pass = pass && testRandVec(10000,10000,1000000,q,true,2,report);
        pass = pass && testRandVec(10000,10000,1000000,q,false,4,report);
        pass = pass && testRandVec(10000,10000,1000000,q,true,4,report);

        pass = pass && testRandVec(100000,100000,1000000,q,false,2,report);
        pass = pass && testRandVec(100000,100000,1000000,q,true,2,report);

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

bool testDenseRand(size_t m, size_t n, size_t nnz, int q, int numThreads, ostream &report)
{
        bool pass=true;

	Field F (q);

        typedef MatrixDomain<Field>::OwnMatrix OwnMatrix;
        LinBox::MatrixDomain<Field> MD(F);
        OwnMatrix MatIn(F,n,m);
        OwnMatrix TestY(F,n,m);
        OwnMatrix RefY(F,n,m);
	OMPBlackbox TestMat(F, n, n);
        SeqBlackbox RefMat(F,n,n);
        Element d;

        TestY.zero();
        RefY.zero();
        MatIn.random();

        typedef pair<size_t,size_t> CoordPair;
        typedef set<CoordPair> PairSet;
        PairSet pairs;
	for(int i = 0; i < (int)nnz; ++i)
	{
                size_t row,col;
                do {
                        row = randRange(0,n);
                        col = randRange(0,n);
                } while (pairs.count(CoordPair(row,col))!=0);

                F.init(d, randRange(0,q));
                TestMat.setEntry(row,col,d);
                RefMat.setEntry(row,col,d);
                pairs.insert(CoordPair(row,col));
        }

        TestMat.finalize();

        omp_set_num_threads(numThreads);

        report << "Testing triplesbb-omp applyLeft against triplesbb applyLeft on random matrix using: " <<
                        "m=" << m << " n=" << n << " nnz=" << nnz << " numthreads=" << numThreads << " q=" << q << std::endl;
        RefMat.applyLeft(RefY,MatIn);
        TestMat.applyLeft(TestY,MatIn);
        if (MD.areEqual(TestY,RefY)) {
                report << "PASS: applyLeft test passed" << std::endl;
        } else {
                pass=false;
                report << "FAILURE: applyLeft test failed" << std::endl;
        }

        return pass;
}

bool testIdent(size_t n, int q, bool shouldFail, int numThreads, ostream &report)
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
                if (shouldFail && (((size_t)i)==(n/3))) {
                        F.init(d,0);
                } else {
                        F.init(d,scaleFactor);
                }
                TestMat.setEntry(i,i,d);
        }

        TestMat.finalize();

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

bool testIdentSuite(int q,ostream &report)
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

        pass = pass && testDenseRand(5,100,1000,65537,1,report);
        pass = pass && testDenseRand(5,100,1000,65537,4,report);

        pass = pass && testRandSuite(65537,report);
        pass = pass && testRandSuite(2,report);

        pass = pass && testIdentSuite(65537,report);
        pass = pass && testIdentSuite(2,report);
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

