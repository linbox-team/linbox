/*
 * benchmarks/benchmark-dense-solve.C
 *
 * Copyright (C) 2019 The LinBox group
 * Author: J-G Dumas
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

/**\file benchmarks/benchmark-dense-solve.C
\brief Solving dense linear system over Q or Zp.
\ingroup benchmarks
*/

#include "linbox/linbox-config.h"
#include <iostream>

#include <givaro/modular.h>
#include "linbox/util/args-parser.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/matrix-stream.h"
#include "linbox/solutions/methods.h"

#ifdef _DEBUG
#define _BENCHMARKS_DEBUG_
#endif

using namespace LinBox;
typedef Givaro::ZRing<Givaro::Integer> Ints;

int main (int argc, char **argv)
{
    Givaro::Integer q = -1 ;
    size_t n = 500 ;
    size_t p = 0;
    size_t bits = 10;

    Argument as[] = {
        { 'q', "-q Q", "Set the field characteristic (-1 for rationals).", TYPE_INTEGER , &q },
        { 'n', "-n N", "Set the matrix dimension.",      TYPE_INT , &n },
        { 'p', "-p P", "0 for sequential, 1 for 2D iterative, 2 for 2D rec, 3 for 2D rec adaptive, 4 for 3D rec in-place, 5 for 3D rec, 6 for 3D rec adaptive.", TYPE_INT , &p },
        { 'b', "-b B", "bit size", TYPE_INT , &bits },
        END_OF_ARGUMENTS
    };

    LinBox::parseArguments(argc,argv,as);

    bool ModComp = false;
    if (q > 0) ModComp = true;

    
    Timer chrono;

    if (ModComp) {
            
        typedef Givaro::Modular<double> Field;
        Field F(q);
	Field::RandIter G(F);
	
#ifdef _BENCHMARKS_DEBUG_
	std::clog << "Setting A ... " << std::endl;
#endif
        chrono.start();		
        DenseMatrix<Field> A(F,n,n);
	PAR_BLOCK { FFLAS::pfrand(F,G, n,n,A.getPointer(),n); }   
        chrono.stop();
#ifdef _BENCHMARKS_DEBUG_
        std::clog << "... A is " << A.rowdim() << " by " << A.coldim() << ", " << chrono << std::endl;
#endif            
        DenseVector<Field> X(F, A.coldim()),B(F, A.rowdim());
	for(auto it=B.begin(); it != B.end(); ++it)
		if (drand48() <0.5)
			F.assign(*it,F.mOne);
                else
			F.assign(*it,F.one);
#ifdef _BENCHMARKS_DEBUG_
        std::clog << "B is [";
        for(DenseVector<Field>::const_iterator it=B.begin();it != B.end(); ++it)
            F.write(std::clog, *it) << " ";
        std::clog << "]" << std::endl;
#endif

            // BlasElimination
        chrono.start();		
        PAR_BLOCK { solve (X, A, B, Method::DenseElimination()); }
        chrono.stop();

#ifdef _BENCHMARKS_DEBUG_
        std::clog << "(BlasElimination) Solution is [";
        for(DenseVector<Field>::const_iterator it=X.begin();it != X.end(); ++it)
            F.write(std::cout, *it) << " ";
        std::clog << "]" << std::endl;		
#endif
        std::cout << "Time: " << chrono.usertime() ;
	FFLAS::writeCommandString(std::cout, as) << std::endl;
    } else { 

        typedef Ints Integers;
        Integers ZZ;
        std::cerr << "Reading A ... " << std::endl;
        chrono.start();		
        size_t n=atoi(argv[1]);
        DenseMatrix<Integers> A(ZZ, n, n);
        DenseMatrix<Integers>::Iterator ait = A.Begin();
        for( ; ait != A.End(); ++ait) {
            double d = drand48()*3;
            if (d < 1)
                *ait = -1;
            else if (d>2)
                *ait = 1;
        }
        chrono.stop();
        std::cerr << "... A is " << A.rowdim() << " by " << A.coldim() << ", " << chrono << std::endl;
            
        Givaro::IntegerDom::Element d;

        DenseVector<Integers> X(ZZ, A.coldim()),B(ZZ, A.rowdim());

        for(DenseVector<Integers>::iterator it=B.begin();
            it != B.end(); ++it)
            if (drand48() <0.5)
                *it = -1;
            else
                *it = 1;

        std::cout << "B is ["; for(DenseVector<Integers>::const_iterator it=B.begin(); it != B.end(); ++it) ZZ.write(std::cout, *it) << " ]" << std::endl;
                
	
            // BlasElimination
        chrono.start();
        solve (X, d, A, B, RingCategories::IntegerTag(), Method::DenseElimination());
        chrono.stop();

        std::cerr << "(BlasElimination) Solution is [";
        for(DenseVector<Integers>::const_iterator it=X.begin();it != X.end(); ++it) ZZ.write(std::cout, *it) << " ";
        ZZ.write(std::cerr << "] / ", d)<< std::endl;
        std::cerr << "CPU time (seconds): " << chrono.usertime() << std::endl;

        std::cerr << "Size in bits: " << Givaro::logtwo(d) << std::endl;
    }

    return 0;
}
