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
    size_t nbiter = 3 ;
    size_t n = 500 ;
    size_t bits = 10;
//     size_t p = 0;

    Argument as[] = {
        { 'i', "-i R", "Set number of repetitions.",       TYPE_INT , &nbiter },
        { 'q', "-q Q", "Set the field characteristic (-1 for rationals).", TYPE_INTEGER , &q },
        { 'n', "-n N", "Set the matrix dimension.",      TYPE_INT , &n },
        { 'b', "-b B", "bit size", TYPE_INT , &bits },
//         { 'p', "-p P", "0 for sequential, 1 for 2D iterative, 2 for 2D rec, 3 for 2D rec adaptive, 4 for 3D rec in-place, 5 for 3D rec, 6 for 3D rec adaptive.", TYPE_INT , &p },
        END_OF_ARGUMENTS
    };

    LinBox::parseArguments(argc,argv,as);

    bool ModComp = false;
    if (q > 0) ModComp = true;

    
    Timer chrono;
    std::vector<std::pair<double,double>> timebits(nbiter);
    for(size_t iter=0; iter<nbiter; ++iter) {

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
            if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::clog <<"A:=",Tag::FileFormat::Maple) << ';' << std::endl;
#endif            
            DenseVector<Field> X(F, A.coldim()),B(F, A.rowdim());
            PAR_BLOCK { FFLAS::pfrand(F,G, n,1,B.getPointer(),1); }
#ifdef _BENCHMARKS_DEBUG_
            std::clog << "B is [";
            for(const auto& it: B)
                F.write(std::clog, it) << ' ';
            std::clog << ']' << std::endl;
#endif

                // DenseElimination
            chrono.start();
            solve (X, A, B, Method::DenseElimination());
            chrono.stop();

#ifdef _BENCHMARKS_DEBUG_
            std::clog << "(DenseElimination) Solution is [";
            for(const auto& it: X)
                F.write(std::clog, it) << ' ';
            std::clog << ']' << std::endl;
#endif
            timebits[iter].second=Givaro::logtwo(GIVMAX(X.front(), 1));
        } else {

            typedef Ints Integers;
            Integers ZZ;
            Integers::RandIter G(ZZ,bits);
#ifdef _BENCHMARKS_DEBUG_
            std::clog << "Reading A ... " << std::endl;
            chrono.start();
#endif
            DenseMatrix<Integers> A(ZZ, n, n);
            PAR_BLOCK { FFLAS::pfrand(ZZ,G, n,n,A.getPointer(),n); }
#ifdef _BENCHMARKS_DEBUG_
            chrono.stop();
            std::clog << "... A is " << A.rowdim() << " by " << A.coldim() << ", " << chrono << std::endl;
            if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::clog <<"A:=",Tag::FileFormat::Maple) << ';' << std::endl;
#endif    
            Givaro::IntegerDom::Element d;

            DenseVector<Integers> X(ZZ, A.coldim()),B(ZZ, A.rowdim());
            PAR_BLOCK { FFLAS::pfrand(ZZ,G, n,1,B.getPointer(),1); }

#ifdef _BENCHMARKS_DEBUG_
            std::clog << "B is [";
            for(const auto& it: B)
                ZZ.write(std::clog, it) << ' ';
            std::clog << ']' << std::endl;
#endif
	
                // DenseElimination
            chrono.start();
            solve (X, d, A, B, RingCategories::IntegerTag(), Method::DenseElimination());
            chrono.stop();

#ifdef _BENCHMARKS_DEBUG_
            std::clog << "(DenseElimination) Solution is [";
            for(const auto& it: X) ZZ.write(std::clog, it) << ' ';
            ZZ.write(std::clog << "] / ", d)<< std::endl;
#endif

            timebits[iter].second=Givaro::logtwo(d);
        }

        timebits[iter].first=chrono.usertime();
    }

#ifdef _BENCHMARKS_DEBUG_
    for(const auto& it: timebits)
        std::clog << it.first << "s, " << it.second << " bits" << std::endl;
#endif

    std::sort(timebits.begin(), timebits.end(),
              [](const std::pair<double,double> & a,
                 const std::pair<double,double> & b) -> bool {
        return a.first > b.first; });

    std::cout << "Time: " << timebits[nbiter/2].first
              << " Bitsize: " << timebits[nbiter/2].second;

    FFLAS::writeCommandString(std::cout, as) << std::endl;

    return 0;
}
