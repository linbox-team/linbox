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
#include "linbox/algorithms/vector-fraction.h"

#ifdef _DEBUG
#define _BENCHMARKS_DEBUG_
#endif

using namespace LinBox;

template<typename Field, typename Vector_t>
std::ostream& printVector(std::ostream& out,
                          const Field& F, const Vector_t& v) {
    out << '[';
    for(const auto& it: v) F.write(std::clog, it) << ',';
    return out << ']';
}
template<typename Vector_t>
double& setBitsize(double& size, const Givaro::Integer& q, const Vector_t& v) {
    return size=Givaro::logtwo(q);
}

typedef Givaro::ZRing<Givaro::Integer> Ints;
typedef VectorFraction<Ints> VectorFractionInts;

template<>
std::ostream& printVector(std::ostream& out, const Ints& Z, const VectorFractionInts& v) {
    return Z.write(printVector(out, Z, v.numer) << " / ", v.denom);
}
template<>
double& setBitsize(double& size, const Givaro::Integer& q, const VectorFractionInts& v) {
    return size=Givaro::logtwo(v.denom);
}


template<typename T1, typename T2, typename T3>
void solve(VectorFractionInts& X, const T1& A, const T2& B, const T3& M) {
    solve(X.numer, X.denom, A, B, M);
}

template<typename Field, typename Vector_t=DenseVector<Field>>
void tmain (std::pair<double,double>& timebits, size_t n,
            const Givaro::Integer& q, size_t bits) {
    Field F(q);						// q is ignored for Integers
    typename Field::RandIter G(F,bits);	// bits is ignored for ModularRandIter

#ifdef _BENCHMARKS_DEBUG_
    std::clog << "Setting A ... " << std::endl;
#endif

    Timer chrono;

    chrono.start();
    DenseMatrix<Field> A(F,n,n);
    PAR_BLOCK { FFLAS::pfrand(F,G, n,n,A.getPointer(),n); }
    chrono.stop();

#ifdef _BENCHMARKS_DEBUG_
    std::clog << "... A is " << A.rowdim() << " by " << A.coldim() << ", " << chrono << std::endl;
    if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::clog <<"A:=",Tag::FileFormat::Maple) << ';' << std::endl;
#endif

    DenseVector<Field> B(F, A.rowdim());
    PAR_BLOCK { FFLAS::pfrand(F,G, n,1,B.getPointer(),1); }

#ifdef _BENCHMARKS_DEBUG_
    printVector(std::clog << "B is ", F, B) << std::endl;
#endif

        // DenseElimination
    Vector_t X(F, A.coldim());
    chrono.start();
    solve (X, A, B, Method::DenseElimination());
    chrono.stop();

#ifdef _BENCHMARKS_DEBUG_
    printVector(std::clog << "(DenseElimination) Solution is ", F, X) << std::endl;
#endif

    setBitsize(timebits.second, q, X);
    timebits.first=chrono.usertime();
}



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

    std::vector<std::pair<double,double>> timebits(nbiter);
    for(size_t iter=0; iter<nbiter; ++iter) {
        if (ModComp) {
            tmain<Givaro::Modular<double>>(timebits[iter],n,q,bits);
        } else {
            tmain<Ints,VectorFractionInts>(timebits[iter],n,q,bits);
        }
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
