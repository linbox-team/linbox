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
#include "linbox/util/timer.h"
#include "linbox/field/modular.h"
#include "linbox/blackbox/triplesbb-omp.h"
#include "linbox/blackbox/triplesbb.h"
#include "linbox/blackbox/transpose.h"
#include "linbox/vector/vector-domain.h"

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
void randMVPTest(BenchmarkFile& of, int p, int numThreads, int m, int n, int nnz, int iters)
{
        typedef typename Blackbox::Field Field;
        typedef typename Field::Element Element;
        typedef std::vector<Element> Vector;

        Field F(p);

        Blackbox A(F,m,n);
        Vector x(n), y(n);
        Element d;

        omp_set_num_threads(numThreads);

	for(int i = 0; i < (int)nnz; ++i)
	{
                size_t row,col;
                row = randRange(0,m);
                col = randRange(0,n);
                F.init(d, randRange(0,p));
                A.setEntry(row,col,d);
	}
        A.finalize();

        for (int i=0;i<(int)n;++i) {
                F.init(x[i],randRange(0,p));
        }

        double start = omp_get_wtime();
        for (int i=0;i<iters;++i) {
                A.apply(y,x);
                A.apply(x,y);
        }
        double time=omp_get_wtime()-start;

        of.addDataField("time",CSDouble(time));
        of.addDataField("num_threads",CSInt(numThreads));
        of.addDataField("GF(p)",CSInt(p));
        of.addDataField("rows",CSInt(m));
        of.addDataField("columns",CSInt(n));
        of.addDataField("nnz",CSInt(nnz));
        of.addDataField("iterations",CSInt(iters));
}

int main(int argc, char **argv)
{
        typedef Modular<double> Field;
        typedef TriplesBBOMP<Field> OMPBB;
        typedef TriplesBB<Field> SeqBB;

        time_t rawTime;
        struct tm *timeInfo;
	srand ((unsigned)time (&rawTime));
        timeInfo=localtime(&rawTime);

        BenchmarkFile benchmarkFile;
        benchmarkFile.addMetadata("problem",CSString("Matrix-Vector Product"));
        benchmarkFile.addMetadata("date",CSDate(*timeInfo));

        benchmarkFile.setType("date", "%a %m/%d %H/%M/%S %Y");
        benchmarkFile.setType("time", "seconds");

        int p=65521;
        int n=50000;
        int m=n;
        int nnz=500000;
        int iters=30;

        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-apply"));
        randMVPTest<OMPBB>(benchmarkFile,p,4,m,n,nnz,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-apply"));
        randMVPTest<OMPBB>(benchmarkFile,p,2,m,n,nnz,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBBOMP-apply"));
        randMVPTest<OMPBB>(benchmarkFile,p,1,m,n,nnz,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.addDataField("algorithm",CSString("TriplesBB-apply"));
        randMVPTest<SeqBB>(benchmarkFile,p,1,m,n,nnz,iters);
        benchmarkFile.pushBackTest();

        benchmarkFile.write(std::cout);

        return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
