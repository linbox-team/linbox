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

#include "linbox/matrix/sparse-matrix.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/vector/blas-vector.h"
#include "linbox/matrix/dense-matrix.h"
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

typedef Modular<double> Field;
typedef TriplesBBOMP<Field> OMPBB;
typedef TriplesBB<Field> SeqBB;
typedef SparseMatrix<Field,SparseMatrixFormat::VPV> VPVBB;

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
void runBlockWiedemann(BenchmarkFile& of,
                       typename Blackbox::Field F,
                       MapSparse<typename Blackbox::Field> matA,
                       MapSparse<typename Blackbox::Field> matU,
                       MapSparse<typename Blackbox::Field> matV,
                       int iters)
{
        typedef typename Blackbox::Field Field;
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

        BlockContainer Sequence(&A,F,U,V);
        BlockMasseyDomain<Field,BlockContainer> MBD(&Sequence);
        std::vector<Block> minpoly;
        std::vector<size_t> degree;

        double start = omp_get_wtime();
        for (int i=0;i<iters;++i) { 
                MBD.left_minpoly_rec(minpoly,degree);
        }
        double time=omp_get_wtime()-start;

        of.addDataField("time",CSDouble(time));
}

template<class Blackbox>
void runMMP(BenchmarkFile& of,
            typename Blackbox::Field F,
            MapSparse<typename Blackbox::Field> leftMat,
            MapSparse<typename Blackbox::Field> rightMat,
            int iters)
{
        typedef typename MatrixDomain<typename Blackbox::Field>::OwnMatrix OwnMatrix;

        int m=leftMat.rowdim();
        int n=leftMat.coldim();
        int p=rightMat.coldim();

        Blackbox A(F,m,n);
        OwnMatrix x(F,n,p), y(F,n,p);
        linbox_check(n==p);

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
void runMVP(BenchmarkFile& of,
            typename Blackbox::Field F,
            MapSparse<typename Blackbox::Field> mat,
            MapSparse<typename Blackbox::Field> vec,
            int iters)
{
        typedef typename Blackbox::Field Field;
        typedef typename Field::Element Element;
        typedef typename std::vector<Element> Vector;

        int m=mat.rowdim();
        int n=mat.coldim();

        Blackbox A(F,m,n);
        Vector x(n), y(n);

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


void blockWiedemannBenchmark(BenchmarkFile& benchmarkFile,
                             int n,
                             int m,
                             int nnz,
                             int iters,
                             int q)
{
        benchmarkFile.addMetadata("problem",CSString("Min-poly"));

        Field F(q);

        MapSparse<Field> A(F,n,n);
        MapSparse<Field> U(F,n,m);
        MapSparse<Field> V(F,m,n);
        MapSparse<Field>::generateRandMat(A,nnz,q);
        MapSparse<Field>::generateDenseRandMat(U,q);
        MapSparse<Field>::generateDenseRandMat(V,q);

        benchmarkFile.addDataField("algorithm",CSString("Block-Wiedemann Parallel"));
        runBlockWiedemann<OMPBB>(benchmarkFile,F,A,U,V,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("Block-Wiedemann Seq"));
        runBlockWiedemann<SeqBB>(benchmarkFile,F,A,U,V,iters);
        benchmarkFile.pushBackTest();
}


void MMPBenchmark(BenchmarkFile& benchmarkFile,
                  int n,
                  int m,
                  int nnz,
                  int iters,
                  int q)
{
        benchmarkFile.addMetadata("problem",CSString("Matrix-Matrix Product"));

        Field F(q);

        MapSparse<Field> A(F,n,n);
        MapSparse<Field> X(F,n,m);
        MapSparse<Field>::generateRandMat(A,nnz,q);
        MapSparse<Field>::generateDenseRandMat(X,q);

        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-applyLeft"));
        runMMP<OMPBB>(benchmarkFile,F,A,X,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBB-seq-applyLeft"));
        runMMP<SeqBB>(benchmarkFile,F,A,X,iters);
        benchmarkFile.pushBackTest();
}


void MVPBenchmark(BenchmarkFile& benchmarkFile,
                  int n,
                  int nnz,
                  int iters,
                  int q)
{
        benchmarkFile.addMetadata("problem",CSString("Matrix-Vector Product"));

        Field F(q);

        MapSparse<Field> A(F,n,n);
        MapSparse<Field> X(F,n,1);
        MapSparse<Field>::generateRandMat(A,nnz,q);
        MapSparse<Field>::generateDenseRandMat(X,q);

        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-apply"));
        runMVP<OMPBB>(benchmarkFile,F,A,X,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-applyLeft"));
        runMVP<OMPBB>(benchmarkFile,F,A,X,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBB-seq-applyLeft"));
        runMVP<SeqBB>(benchmarkFile,F,A,X,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBB-seq-apply"));
        runMVP<SeqBB>(benchmarkFile,F,A,X,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("SparseMatrix::VPV"));
        runMVP<VPVBB>(benchmarkFile,F,A,X,iters);
        benchmarkFile.pushBackTest();
}

int main(int argc, char **argv)
{
        int numThreads=1,n=50000,m=100,nnz=50000,q=65521,iters=1;
        bool runMMPBenchmarks=false,runMVPBenchmarks=false,runWiedemann=false;

	static Argument args[] = {
		{ 't', "-t THREADS", "Number of threads", TYPE_INT, &numThreads },
                { 'n', "-n N", "Dimension of N*N matrix A", TYPE_INT, &n },
                { 'm', "-m M", "Dimension of N*M matrix X", TYPE_INT, &m },
                { 'z', "-z NNZ", "Number of non-zero entries", TYPE_INT, &nnz},
                { 'q', "-q PRIME", "Use field GF(Q) for prime Q", TYPE_INT, &q},
                { 'i', "-i ITERS", "Number of iterations", TYPE_INT, &iters},

                { 'a', NULL, "Run Matrix-Matrix Benchmarks", TYPE_BOOL,&runMMPBenchmarks},
                { 'b', NULL, "Run Matrix-Vector Benchmarks", TYPE_BOOL,&runMVPBenchmarks},
                { 'c', NULL, "Run Block-Wiedemann Benchmarks", TYPE_BOOL,&runWiedemann},
                END_OF_ARGUMENTS};

	parseArguments (argc, argv, args);

        BenchmarkFile benchmarkFile;

        benchmarkFile.addMetadata("Field-Implementation",CSString("Modular<double>"));
        benchmarkFile.addMetadata("num_threads",CSInt(numThreads));
        benchmarkFile.addMetadata("N*N Blackbox Dimensions",CSInt(n));
        benchmarkFile.addMetadata("Fat Vector Dimension",CSInt(m));
        benchmarkFile.addMetadata("nnz",CSInt(nnz));
        benchmarkFile.addMetadata("iterations",CSInt(iters));
        benchmarkFile.addMetadata("date",BenchmarkFile::getDateStamp());
        benchmarkFile.addMetadata("GF(q)",CSInt(q));

        benchmarkFile.setType("date",BenchmarkFile::getDateFormat());
        benchmarkFile.setType("time","seconds");

        omp_set_num_threads(numThreads);

        if (runMMPBenchmarks) {
                MMPBenchmark(benchmarkFile,n,m,nnz,iters,q);
                benchmarkFile.write(std::cout);
        }
        if (runMVPBenchmarks) {
                MVPBenchmark(benchmarkFile,n,nnz,iters,q);
                benchmarkFile.write(std::cout);
        }
        if (runWiedemann) {
                blockWiedemannBenchmark(benchmarkFile,n,m,nnz,iters,q);
                benchmarkFile.write(std::cout);
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
