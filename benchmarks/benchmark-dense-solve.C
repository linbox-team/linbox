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

#include "linbox/algorithms/vector-fraction.h"
#include "linbox/matrix/sparse-matrix.h"
#include "linbox/solutions/methods.h"
#include "linbox/solutions/solve.h"
#include "linbox/util/args-parser.h"
#include "linbox/util/matrix-stream.h"
#include <givaro/modular.h>

#ifdef _DEBUG
#define _BENCHMARKS_DEBUG_
#endif

using namespace LinBox;

using Ints = Givaro::ZRing<Givaro::Integer>;
using Mods = Givaro::Modular<double>;
using VectorFractionInts = VectorFraction<Ints>;

namespace {
    struct Arguments {
        Givaro::Integer q = -1;
        int nbiter = 3;
        int n = 500;
        int bits = 10;
        int seed = -1;
        std::string dispatchString = "Auto";
        std::string methodString = "Auto";
    };

    template <typename Vector>
    double& setBitsize(double& size, const Givaro::Integer& q, const Vector& v)
    {
        return size = Givaro::logtwo(q);
    }

    template <>
    double& setBitsize(double& size, const Givaro::Integer& q, const VectorFractionInts& v)
    {
        return size = Givaro::logtwo(v.denom);
    }
}

template <typename Field, typename Vector>
void benchmark(typename Field::RandIter &randIter, std::array<double, 3>& timebits, Arguments& args, MethodBase& method)
{
    const Field &F = randIter.ring();

#ifdef _BENCHMARKS_DEBUG_
    std::clog << "Setting A ... " << std::endl;
#endif

    DenseMatrix<Field> A(F, args.n, args.n);
    DenseVector<Field> B(F, A.rowdim());
    Timer chrono;

    if (method.master()) {
        chrono.start();
        PAR_BLOCK { FFLAS::pfrand(F, randIter, args.n, args.n, A.getPointer(), args.n); }
        chrono.stop();

#ifdef _BENCHMARKS_DEBUG_
        std::clog << "... A is " << A.rowdim() << " by " << A.coldim() << ", " << chrono << std::endl;
        if (A.rowdim() <= 20 && A.coldim() <= 20) A.write(std::clog << "A:=", Tag::FileFormat::Maple) << ';' << std::endl;
#endif

        PAR_BLOCK { FFLAS::pfrand(F, randIter, args.n, 1, B.getPointer(), 1); }

#ifdef _BENCHMARKS_DEBUG_
        std::clog << "B is " << B << std::endl;
#endif
    }

    method.pCommunicator->bcast(A, 0);
    method.pCommunicator->bcast(B, 0);

    Vector X(F, A.coldim());

    if (method.master()) {
        chrono.start();
    }

    if (args.methodString == "Elimination")                 solve(X, A, B, Method::Elimination(method));
    else if (args.methodString == "DenseElimination")       solve(X, A, B, Method::DenseElimination(method));
    else if (args.methodString == "SparseElimination")      solve(X, A, B, Method::SparseElimination(method));
    else if (args.methodString == "Dixon")                  solve(X, A, B, Method::Dixon(method));
    else if (args.methodString == "CRA")                    solve(X, A, B, Method::CRAAuto(method));
    else if (args.methodString == "SymbolicNumericOverlap") solve(X, A, B, Method::SymbolicNumericOverlap(method));
    else if (args.methodString == "SymbolicNumericNorm")    solve(X, A, B, Method::SymbolicNumericNorm(method));
    else if (args.methodString == "Blackbox")               solve(X, A, B, Method::Blackbox(method));
    else if (args.methodString == "Wiedemann")              solve(X, A, B, Method::Wiedemann(method));
    else if (args.methodString == "Lanczos")                solve(X, A, B, Method::Lanczos(method));
    // @fixme Won't compile
    // else if (args.methodString == "BlockLanczos")           solve(X, A, B, Method::BlockLanczos(method));
    else                                                    solve(X, A, B, Method::Auto(method));

    if (method.master()) {
        chrono.stop();

        timebits[0] = chrono.usertime();
        timebits[1] = chrono.realtime();
        setBitsize(timebits[2], args.q, X);
    }
}

int main(int argc, char** argv)
{
    Arguments args;
    Argument as[] = {{'i', "-i", "Set number of repetitions.", TYPE_INT, &args.nbiter},
                     {'q', "-q", "Set the field characteristic (-1 for rationals).", TYPE_INTEGER, &args.q},
                     {'n', "-n", "Set the matrix dimension.", TYPE_INT, &args.n},
                     {'b', "-b", "bit size", TYPE_INT, &args.bits},
                     {'s', "-s", "Seed for randomness.", TYPE_INT, &args.seed},
                     {'d', "-d", "Dispatch mode (any of: Auto, Sequential, SMP, Distributed).", TYPE_STR, &args.dispatchString},
                     {'M', "-M",
                      "Choose the solve method (any of: Auto, Elimination, DenseElimination, SparseElimination, "
                      "Dixon, CRA, SymbolicNumericOverlap, SymbolicNumericNorm, "
                      "Blackbox, Wiedemann, Lanczos).",
                      TYPE_STR, &args.methodString},
                     END_OF_ARGUMENTS};
    LinBox::parseArguments(argc, argv, as);

    if (args.seed < 0) {
        args.seed = time(nullptr);
    }

    // Setting up context

    Communicator communicator(&argc, &argv);
    if (communicator.master()) {
        std::clog << "Communicator size: " << communicator.size() << std::endl;
    }

    MethodBase method;
    method.pCommunicator = &communicator;
    if (args.dispatchString == "Sequential")        method.dispatch = Dispatch::Sequential;
    else if (args.dispatchString == "SMP")          method.dispatch = Dispatch::SMP;
    else if (args.dispatchString == "Distributed")  method.dispatch = Dispatch::Distributed;
    else                                            method.dispatch = Dispatch::Auto;

    // Real benchmark

    bool isModular = false;
    if (args.q > 0) isModular = true;

    using Timing = std::array<double, 3>;
    std::vector<Timing> timebits(args.nbiter);
    for (int iter = 0; iter < args.nbiter; ++iter) {
        if (isModular) {
            Mods Fq(args.q);
            Mods::RandIter randIterFq(Fq, args.seed);
            benchmark<Mods, DenseVector<Mods>>(randIterFq, timebits[iter], args, method);
        }
        else {
            Ints ZZ;
            Ints::RandIter randIterZZ(ZZ, args.seed);
            randIterZZ.setBitsize(args.bits);
            benchmark<Ints, VectorFractionInts>(randIterZZ, timebits[iter], args, method);
        }
    }

#ifdef _BENCHMARKS_DEBUG_
    for (const auto& it : timebits) std::clog << it[0] << "s, " << it[2] << " bits" << std::endl;
#endif

    if (method.master()) {
        std::sort(timebits.begin(), timebits.end(), [](const Timing& a, const Timing& b) -> bool { return a[0] > b[0]; });

        std::cout << "UserTime: " << timebits[args.nbiter / 2][0];
        std::cout << " RealTime: " << timebits[args.nbiter / 2][1];
        std::cout << " Bitsize: " << timebits[args.nbiter / 2][2];

        FFLAS::writeCommandString(std::cout, as) << std::endl;
    }

    return 0;
}
