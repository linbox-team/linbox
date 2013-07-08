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

#include <omp.h>

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

void runOMPTest(std::ostream &out,int numThreads, integer q, int n, int m, int nnz, int iters)
{
        typedef Modular<double> Field;
	typedef Field::Element Element;
	typedef std::vector <Element> Vector;
	typedef TriplesBBOMP<MatrixDomain<Field> > Blackbox;

	Field F (q);

	Blackbox A(F, m, n);
        Vector x(n), y(n);
	Element d;

	for(int i = 0; i < (int)nnz; ++i)
	{
                size_t row,col;
                row = randRange(0,m);
                col = randRange(0,n);
                F.init(d, randRange(0,q));
                A.setEntry(row,col,d);
        }

        for (int i=0;i<(int)n;++i) {
                F.init(x[i],randRange(0,q));
        }
        omp_set_num_threads(numThreads);
        A.sortBlock();
        A.sortRow();
        double start = omp_get_wtime();
        for (int i=0;i<iters;++i) {
                A.apply(y,x);
                A.apply(x,y);
        }
        double time=omp_get_wtime()-start;
        out.precision(10);
        out << time << "," << numThreads << "," << q << "," << n << "," << m << "," << nnz << "," << iters << "," << "omp" << std::endl;
}

void runSeqTest(std::ostream &out, integer q, int n, int m, int nnz, int iters)
{
        typedef Modular<double> Field;
	typedef Field::Element Element;
	typedef std::vector <Element> Vector;
	typedef TriplesBB<MatrixDomain<Field> > Blackbox;

	Field F (q);

	Blackbox A(F, m, n);
        Vector x(n), y(n);
	Element d;


	for(int i = 0; i < (int)nnz; ++i)
	{
                size_t row,col;
                row = randRange(0,m);
                col = randRange(0,n);
                F.init(d, randRange(0,q));
                A.setEntry(row,col,d);
	}

        for (int i=0;i<(int)n;++i) {
                F.init(x[i],randRange(0,q));
        }
        double start = omp_get_wtime();
        for (int i=0;i<iters;++i) {
                A.apply(y,x);
                A.apply(x,y);
        }
        double time=omp_get_wtime()-start;
        out.precision(10);
        out << time << "," << 1 << "," << q << "," << n << "," << m << "," << nnz << "," << iters << "," << "seq" << std::endl;
}


void printHeader(std::ostream &out)
{
        out << "time,num_threads,field_size,n,m,nnz,code,iters" << std::endl;
}

int main(int argc, char **argv)
{
	srand ((unsigned)time (NULL));

        printHeader(std::cout);
        runOMPTest(std::cout,4,2147483629,500000,500000,5000000,3);
        runOMPTest(std::cout,2,2147483629,500000,500000,5000000,3);
        runOMPTest(std::cout,1,2147483629,500000,500000,5000000,3);
        runSeqTest(std::cout,2147483629,500000,500000,5000000,3);
        return 0;
}

// Local Variables:
// mode: C++
// tab-width: 8
// indent-tabs-mode: nil
// c-basic-offset: 8
// End:
// vim:sts=8:sw=8:ts=8:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
