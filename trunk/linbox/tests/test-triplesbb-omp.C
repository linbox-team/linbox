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
#include "linbox/matrix/sparse-matrix.h"
// #include "linbox/blackbox/triplesbb-omp.h"
// #include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/vector/vector-domain.h"
#include "linbox/matrix/dense-matrix.h"

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
        int offset = (int)(rangeSize*normedRVal);
        return start+offset;
}

template<class Field>
void varyMatMatPair(MapSparse<Field>& leftMat, MapSparse<Field>& rightMat, int q)
{
        typedef typename Field::Element Element;

        Field F(q);
        Element d;

        int m=(int)leftMat.rowdim();
        int n=(int)leftMat.coldim();
        int p=(int)rightMat.coldim();

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

template<class Field,class Blackbox>
void computeAnswer(const MapSparse<Field>& A,
                   const MapSparse<Field>& X,
                   MapSparse<Field>& Y,
                   const Field& F,
                   bool isRight,
                   bool useVector)
{
	typedef BlasVector<Field> Vector ;
        typedef typename MatrixDomain<Field>::OwnMatrix OwnMatrix;

        LinBox::MatrixDomain<Field> MD(F);

        if (useVector) {
                if (isRight) {
                        Vector vecX(F,X.coldim()),vecY(F,A.coldim());
                        X.toVector(vecX);
                        Blackbox matA(F);
                        A.copy(matA);
                        matA.applyTranspose(vecY,vecX);
                        Y.fromVector(vecY,1,A.coldim());
                } else {
                        Vector vecX(F,X.rowdim()),vecY(F,A.rowdim());
                        X.toVector(vecX);
                        Blackbox matA(F);
                        A.copy(matA);
                        matA.apply(vecY,vecX);
                        Y.fromVector(vecY,A.rowdim(),1);
                }
        } else {
                OwnMatrix matX(F,X.rowdim(),X.coldim());
                if (isRight) {
                        OwnMatrix matY(F,X.rowdim(),A.coldim());
                        X.copy(matX);
                        Blackbox matA(F);
                        A.copy(matA);
                        matA.applyRight(matY,matX);
                        Y.copyFrom(matY);
                } else {
                        OwnMatrix matY(F,A.rowdim(),X.coldim());
                        X.copy(matX);
                        Blackbox matA(F);
                        A.copy(matA);
                        matA.applyLeft(matY,matX);
                        Y.copyFrom(matY);
                }
        }
}

bool outputResult(std::string testName,
                  bool expected,
                  bool actual,
                  int m, int n, int p, int nnz, int numThreads, int q,
                  ostream& report)
{
        report << "Testing " << testName << " using: " <<
                "m=" << m << " n=" << n << " p=" << p << " nnz=" << nnz << " numthreads=" << numThreads << " q=" << q << std::endl;

        if (expected==actual) {
                report << "PASS: ";
        } else {
                report << "FAILURE: ";
        }
        report << testName << " expected " << expected << " got " << actual << std::endl;

        return expected==actual;
}

template<class Field>
void reshapeAndScale(MapSparse<Field>& matOut,
                     MapSparse<Field>& matIn,
                     int m, int n, int alpha)
{
        typedef typename Field::Element Element;

        matOut.shape(m,n);
        m=(m<(int)matIn.rowdim())?m:(int)matIn.rowdim();
        n=(n<(int)matIn.coldim())?n:(int)matIn.coldim();
        Element d,a,t;
        matIn.field().init(a,alpha);
        for (int i=0;i<m;++i) {
                for (int j=0;j<n;++j) {
                        d=matIn.getEntry(i,j);
                        matIn.field().mul(t,d,a);
                        matOut.setEntry(i,j,t);
                }
        }
}

template <class Field>
bool runIdentTest(int n,
                  int m,
                  int p,
                  int q,
                  int numThreads,
                  bool shouldFail,
                  bool isRight,
                  bool useVector,
                  ostream& report)
{
        typedef SparseMatrix<Field,SparseMatrixFormat::TPL_omp> OMPBlackbox;

        omp_set_num_threads(numThreads);

        Field F(q);
        MapSparse<Field> A(F,n,m),X,YTest,YRef;
        if (isRight) {
                X.init(F,p,n);
                YTest.init(F,p,m);
                YRef.init(F,p,m);
        } else {
                X.init(F,m,p);
                YTest.init(F,n,p);
                YRef.init(F,n,p);
        }
        int alpha=randRange(1,q);
        MapSparse<Field>::generateScaledIdent(A,alpha);
        MapSparse<Field>::generateDenseRandMat(X,q);
        if (isRight) {
                reshapeAndScale<Field>(YRef,X,p,m,alpha);
        } else {
                reshapeAndScale<Field>(YRef,X,n,p,alpha);
        }

        if (isRight && shouldFail) {
                varyMatMatPair<Field>(X,A,q);
        } else if ((!isRight) && shouldFail) {
                varyMatMatPair<Field>(A,X,q);
        }

        computeAnswer<Field,OMPBlackbox>(A,X,YTest,F,isRight,useVector);
        std::string testName;
        if ((!isRight) && useVector) {
                testName="ident apply()";
        } else if (isRight && (!useVector)) {
                testName="ident applyRight()";
        } else if (isRight && useVector) {
                testName="ident applyTranspose()";
        } else {
                testName="ident applyLeft()";
        }
        int nnz=(n<m)?n:m;
        return outputResult(testName,!shouldFail,YRef.areEqual(YTest),m,n,p,nnz,numThreads,q,report);
}

template <class Field>
bool runRandTest(int n,
                 int m,
                 int p,
                 double density,
                 int q,
                 int numThreads,
                 bool shouldFail,
                 bool isRight,
                 bool useVector,
                 ostream& report)
{
        typedef SparseMatrix<Field,SparseMatrixFormat::TPL_omp> OMPBlackbox;
        typedef SparseMatrix<Field,SparseMatrixFormat::TPL> SeqBlackbox;

        omp_set_num_threads(numThreads);

        double matSize=n*m;
        int nnz=(int) (density*matSize);

        Field F(q);
        MapSparse<Field> A(F,n,m),X,YTest,YRef;
        if (isRight) {
                X.init(F,p,n);
                YTest.init(F,p,m);
                YRef.init(F,p,m);
        } else {
                X.init(F,m,p);
                YTest.init(F,n,p);
                YRef.init(F,n,p);
        }
        MapSparse<Field>::generateRandMat(A,nnz,q);
        MapSparse<Field>::generateDenseRandMat(X,q);

        computeAnswer<Field,SeqBlackbox>(A,X,YRef,F,isRight,useVector);

        if (isRight && shouldFail) {
                varyMatMatPair<Field>(X,A,q);
        } else if ((!isRight) && shouldFail) {
                varyMatMatPair<Field>(A,X,q);
        }

        computeAnswer<Field,OMPBlackbox>(A,X,YTest,F,isRight,useVector);

        std::string testName;
        if ((!isRight) && useVector) {
                testName="rand apply()";
        } else if (isRight && (!useVector)) {
                testName="rand applyRight()";
        } else if (isRight && useVector) {
                testName="rand applyTranspose()";
        } else {
                testName="rand applyLeft()";
        }
        return outputResult(testName,!shouldFail,YRef.areEqual(YTest),m,n,p,nnz,numThreads,q,report);
}


template<class Field>
bool runSizeSuite(int n, int m, int p,
                  double density,
                  std::vector<int>& qs,
                  std::vector<int>& numThreads,
                  bool shouldFail, bool isRight,
                  ostream& report)
{
        bool pass=true;
        for (size_t i=0;i<qs.size();++i) {
                for (size_t j=0;j<numThreads.size();++j) {
                        if (p==1) {
                                bool useVector=true;
                                pass = pass && runRandTest<Field>
                                        (n,m,p,density,qs[i],numThreads[j],
                                         shouldFail,isRight,useVector,report);
                                pass = pass && runIdentTest<Field>
                                        (n,m,p,qs[i],numThreads[j],
                                        shouldFail,isRight,useVector,report);
                        }
                        pass = pass && runRandTest<Field>
                                (n,m,p,density,qs[i],numThreads[j],
                                 shouldFail,isRight,false,report);
                        pass = pass && runIdentTest<Field>
                                (n,m,p,qs[i],numThreads[j],
                                 shouldFail,isRight,false,report);
                }
        }
        return pass;
}

template<class Field>
bool testSuite(std::vector<int>& qs, ostream &report)
{
        bool pass=true;

        std::vector<int> numThreads;
        numThreads.push_back(1);
        numThreads.push_back(2);
        numThreads.push_back(3);
        numThreads.push_back(4);
        bool shouldFail=false,isRight=false;
        do {
                pass=pass&&runSizeSuite<Field>(10000,10,1,0.01,qs,numThreads,shouldFail,isRight,report);
                pass=pass&&runSizeSuite<Field>(10000,10,10,0.01,qs,numThreads,shouldFail,isRight,report);
                pass=pass&&runSizeSuite<Field>(10,10000,1,0.01,qs,numThreads,shouldFail,isRight,report);
                pass=pass&&runSizeSuite<Field>(10,10000,10,0.01,qs,numThreads,shouldFail,isRight,report);
                pass=pass&&runSizeSuite<Field>(1,10000,1,0.1,qs,numThreads,shouldFail,isRight,report);

                pass=pass&&runSizeSuite<Field>(10000,10000,1,0.0,qs,numThreads,shouldFail,isRight,report);

                pass=pass&&runSizeSuite<Field>(2,2,1,1.0,qs,numThreads,shouldFail,isRight,report);
                pass=pass&&runSizeSuite<Field>(2,2,10000,1.0,qs,numThreads,shouldFail,isRight,report);

                pass=pass&&runSizeSuite<Field>(10000,10000,1,0.0001,qs,numThreads,shouldFail,isRight,report);
                pass=pass&&runSizeSuite<Field>(10000,10000,10,0.0001,qs,numThreads,shouldFail,isRight,report);

                shouldFail^=isRight;
                isRight=!isRight;
        } while (shouldFail || isRight);

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

        std::vector<int> qs;
        qs.push_back(2);
        qs.push_back(65537);

        pass = testSuite<Modular<double> >(qs,report);

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

