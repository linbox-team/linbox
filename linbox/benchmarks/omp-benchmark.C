/* Copyright (C) 2013 LinBox
 * Written by AJS <stachnik@udel.edu>
 *
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

/*! @file   benchmarks/omp-benchmark.C
 * @ingroup benchmarks
 * @brief f
 */

#include "linbox/linbox-config.h"

#include <stdlib.h>
#include <fstream>
#include <time.h>
#include <omp.h>

#include "benchmarks/CSValue.h"
#include "benchmarks/BenchmarkFile.h"

#include "linbox/vector/blas-vector.h"
#include "linbox/matrix/blas-matrix.h"
#include "linbox/util/timer.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/triplesbb-omp.h"
#include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/vector/vector-domain.h"
#include "examples/map-sparse.h"

#include "linbox/algorithms/blackbox-block-container.h"
#include "linbox/algorithms/block-massey-domain.h"

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

template<class Blackbox>
void blockWiedemannBenchmark(BenchmarkFile& of,
                             typename Blackbox::Field F,
                             MapSparse<typename Blackbox::Field> matA,
                             MapSparse<typename Blackbox::Field> matU,
                             MapSparse<typename Blackbox::Field> matV,
                             int numThreads)
{
        typedef typename Blackbox::Field Field;
        typedef typename Field::Element Element;
        typedef BlackboxBlockContainer<Field,Blackbox> BlockContainer;
        typedef BlasMatrix<Field> Block;

        int m=matA.rowdim();
        int n=matA.coldim();
        int p=matU.rowdim();
        int q=matV.coldim();

        Blackbox A(F,m,n);
        Block U(F,p,m), V(F,n,q);
        matU.copy(U);
        matV.copy(V);
        matA.copy(A);
        A.finalize();

        omp_set_num_threads(numThreads);

        BlockContainer Sequence(&A,F,U,V);
        BlockMasseyDomain<Field,BlockContainer> MBD(&Sequence);
        std::vector<Block> minpoly;
        std::vector<size_t> degree;

        double start = omp_get_wtime();
        MBD.left_minpoly_rec(minpoly,degree);
        double time=omp_get_wtime()-start;

        of.addDataField("time",CSDouble(time));
}

void randBlockWiedemannBenchmark(int numThreads,
                                 int m,
                                 int n,
                                 int p,
                                 int q,
                                 int nnz,
                                 int fieldPrime)
{
        typedef Modular<double> Field;
        typedef TriplesBBOMP<Field> OMPBB;
        typedef TriplesBB<Field> SeqBB;

        BenchmarkFile benchmarkFile;
        benchmarkFile.addMetadata("problem",CSString("Min-poly"));
        benchmarkFile.addMetadata("date",BenchmarkFile::getDateStamp());
        benchmarkFile.addMetadata("GF(q)",CSInt(fieldPrime));
        benchmarkFile.addMetadata("left-matrix-rows",CSInt(n));
        benchmarkFile.addMetadata("left-matrix-columns",CSInt(n));
        benchmarkFile.addMetadata("right-matrix-columns",CSInt(p));
        benchmarkFile.addMetadata("nnz",CSInt(nnz));
        benchmarkFile.setType("date",BenchmarkFile::getDateFormat());
        benchmarkFile.setType("time","seconds");

        Field F(fieldPrime);

        MapSparse<Field> A(F,m,n);
        MapSparse<Field> U(F,p,m);
        MapSparse<Field> V(F,n,q);
        MapSparse<Field>::generateRandMat(A,nnz,fieldPrime);
        MapSparse<Field>::generateDenseRandMat(U,fieldPrime);
        MapSparse<Field>::generateDenseRandMat(V,fieldPrime);

        benchmarkFile.addDataField("algorithm",CSString("Block-Wiedemann Parallel"));
        benchmarkFile.addDataField("num_threads",CSInt(numThreads));
        blockWiedemannBenchmark<OMPBB>(benchmarkFile,F,A,U,V,numThreads);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("Block-Wiedemann Seq"));
        benchmarkFile.addDataField("num_threads",CSInt(numThreads));
        blockWiedemannBenchmark<SeqBB>(benchmarkFile,F,A,U,V,numThreads);
        benchmarkFile.pushBackTest();

        benchmarkFile.write(std::cout);
}
                                 

template<class Blackbox>
void randMMPTest(BenchmarkFile& of,
                 typename Blackbox::Field F,
                 MapSparse<typename Blackbox::Field> leftMat,
                 MapSparse<typename Blackbox::Field> rightMat,
                 int numThreads, int iters)
{
        typedef typename Blackbox::Field Field;
        typedef typename Field::Element Element;
        typedef typename std::vector<Element> Vector;
        typedef typename MatrixDomain<Field>::OwnMatrix OwnMatrix;

        int m=leftMat.rowdim();
        int n=leftMat.coldim();
        int p=rightMat.coldim();

        Blackbox A(F,m,n);
        OwnMatrix x(F,n,p), y(F,n,p);
        linbox_check(n==p);

        omp_set_num_threads(numThreads);

        leftMat.copy(A);
        rightMat.copy(x);
        A.finalize();

        double start = omp_get_wtime();
        for (int i=0;i<iters;++i) {
                A.applyLeft(y,x);
                A.applyLeft(x,y);
        }
        double time=omp_get_wtime()-start;

        of.addDataField("time",CSDouble(time));
}

template<class Blackbox>
void randMVPTest(BenchmarkFile& of,
                 typename Blackbox::Field F,
                 MapSparse<typename Blackbox::Field> mat,
                 MapSparse<typename Blackbox::Field> vec,
                 int numThreads, int iters)
{
        typedef typename Blackbox::Field Field;
        typedef typename Field::Element Element;
        typedef typename std::vector<Element> Vector;

        int m=mat.rowdim();
        int n=mat.coldim();

        Blackbox A(F,m,n);
        Vector x(n), y(n);

        omp_set_num_threads(numThreads);

        mat.copy(A);
        vec.toVector(x);
        A.finalize();

        double start = omp_get_wtime();
        for (int i=0;i<iters;++i) {
                A.apply(y,x);
                A.apply(x,y);
        }
        double time=omp_get_wtime()-start;

        of.addDataField("time",CSDouble(time));
}

void MMPBenchmark(int numThreads,
                  int n,
                  int p,
                  int nnz,
                  int iters,
                  int q)
{
        typedef Modular<double> Field;
        typedef TriplesBBOMP<Field> OMPBB;
        typedef TriplesBB<Field> SeqBB;

        BenchmarkFile benchmarkFile;
        benchmarkFile.addMetadata("problem",CSString("Matrix-Matrix Product"));
        benchmarkFile.addMetadata("date",BenchmarkFile::getDateStamp());
        benchmarkFile.addMetadata("GF(q)",CSInt(q));
        benchmarkFile.addMetadata("left-matrix-rows",CSInt(n));
        benchmarkFile.addMetadata("left-matrix-columns",CSInt(n));
        benchmarkFile.addMetadata("right-matrix-columns",CSInt(p));
        benchmarkFile.addMetadata("nnz",CSInt(nnz));
        benchmarkFile.addMetadata("iterations",CSInt(iters));
        benchmarkFile.setType("date",BenchmarkFile::getDateFormat());
        benchmarkFile.setType("time","seconds");

        Field F(q);

        MapSparse<Field> A(F,n,n);
        MapSparse<Field> X(F,n,p);
        MapSparse<Field>::generateRandMat(A,nnz,q);
        MapSparse<Field>::generateDenseRandMat(X,q);

        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-applyLeft"));
        benchmarkFile.addDataField("num_threads",CSInt(numThreads));
        randMMPTest<OMPBB>(benchmarkFile,F,A,X,numThreads,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBB-seq-applyLeft"));
        randMMPTest<SeqBB>(benchmarkFile,F,A,X,1,iters);
        benchmarkFile.pushBackTest();
        benchmarkFile.write(std::cout);
}


void MVPBenchmarkSuite(int numThreads,
                       int n,
                       int nnz,
                       int iters,
                       int q)
{
        typedef Modular<double> Field;
        typedef TriplesBBOMP<Field> OMPBB;
        typedef TriplesBB<Field> SeqBB;

        int m=n;
        Field F(q);

        BenchmarkFile benchmarkFile;
        benchmarkFile.addMetadata("problem",CSString("Matrix-Vector Product"));
        benchmarkFile.addMetadata("date",BenchmarkFile::getDateStamp());
        benchmarkFile.addMetadata("GF(q)",CSInt(q));
        benchmarkFile.addMetadata("rows",CSInt(m));
        benchmarkFile.addMetadata("columns",CSInt(n));
        benchmarkFile.addMetadata("nnz",CSInt(nnz));
        benchmarkFile.addMetadata("iterations",CSInt(iters));
        benchmarkFile.setType("date",BenchmarkFile::getDateFormat());
        benchmarkFile.setType("time","seconds");

        MapSparse<Field> A(F,m,n);
        MapSparse<Field> X(F,n,1);
        MapSparse<Field>::generateRandMat(A,nnz,q);
        MapSparse<Field>::generateDenseRandMat(X,q);
        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-apply"));
        benchmarkFile.addDataField("num_threads",CSInt(numThreads));
        randMVPTest<OMPBB>(benchmarkFile,F,A,X,numThreads,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-applyLeft"));
        benchmarkFile.addDataField("num_threads",CSInt(numThreads));
        randMMPTest<OMPBB>(benchmarkFile,F,A,X,numThreads,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBB-seq-applyLeft"));
        benchmarkFile.addDataField("num_threads",CSInt(numThreads));
        randMMPTest<SeqBB>(benchmarkFile,F,A,X,numThreads,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBB-seq-apply"));
        randMVPTest<SeqBB>(benchmarkFile,F,A,X,1,iters);
        benchmarkFile.pushBackTest();
        benchmarkFile.write(std::cout);
}

int main(int argc, char **argv)
{
        int numThreads=1,n=50000,p=100,iters=30,q=65521,m=100,r=100,nnz=50000;
        bool runMMPBenchmarks=false,runMVPBenchmarks=true,runWiedemann=false;

	static Argument args[] = {
		{ 't', "-t THREADS", "Number of threads", TYPE_INT, &numThreads },
                { 'm', "-m M", "Dimension", TYPE_INT, &m },
                { 'r', "-r R", "Dimension", TYPE_INT, &r },
                { 'n', "-n N", "Dimension of N*N matrix", TYPE_INT, &n },
                { 'p', "-p P", "Dimension of N*P matrix", TYPE_INT, &p },
                { 'i', "-i ITERS", "Number of iterations", TYPE_INT, &iters},
                { 'q', "-q PRIME", "Use field GF(Q) for prime Q", TYPE_INT, &q},
                { 'z', "-z NNZ", "Number of non-zero entries", TYPE_INT, &nnz},
                { 'a', NULL, "Run Matrix-Matrix Benchmarks", TYPE_BOOL,&runMMPBenchmarks},
                { 'b', NULL, "Run Matrix-Vector Benchmarks", TYPE_BOOL,&runMVPBenchmarks},
                { 'c', NULL, "Run Block-Wiedemann Benchmarks", TYPE_BOOL,&runWiedemann},
                END_OF_ARGUMENTS};

	parseArguments (argc, argv, args);

        if (runMMPBenchmarks) {
                MMPBenchmark(numThreads,n,p,nnz,iters,q);
        }
        if (runMVPBenchmarks) {
                MVPBenchmarkSuite(numThreads,n,nnz,iters,q);
        }
        if (runWiedemann) {
                randBlockWiedemannBenchmark(numThreads,m,n,p,r,nnz,q);
        }
        return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
