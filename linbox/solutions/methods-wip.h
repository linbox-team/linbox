/*
 * Copyright(C) LinBox
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
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 * ========LICENCE========
 */

#pragma once

#include <linbox/util/mpicpp.h>
#include <string>

#define DEFINE_METHOD_CONTENT(MethodName)                                                                                        \
    MethodName() = default;                                                                                                      \
    MethodName(const MethodName&) = default;                                                                                     \
    MethodName(const MethodBase& methodBase)                                                                                     \
        : MethodBase(methodBase)                                                                                                 \
    {                                                                                                                            \
    }

#define DEFINE_METHOD(_MethodName, _CategoryTag)                                                                                 \
    struct _MethodName : public MethodBase {                                                                                     \
        using CategoryTag = _CategoryTag;                                                                                        \
        static std::string name() { return std::string("Method::") + #_MethodName; }                                             \
        DEFINE_METHOD_CONTENT(_MethodName)                                                                                       \
    };

#define DEFINE_COMPOUND_METHOD(_MethodName, _CategoryTag)                                                                        \
    template <class IterationMethod>                                                                                             \
    struct _MethodName : public MethodBase {                                                                                     \
        using CategoryTag = _CategoryTag;                                                                                        \
        IterationMethod iterationMethod;                                                                                         \
        static std::string name() { return std::string("Method::") + #_MethodName "<" + IterationMethod::name() + ">"; }         \
        DEFINE_METHOD_CONTENT(_MethodName)                                                                                       \
    };                                                                                                                           \
    using _MethodName##Auto = _MethodName<MethodWIP::Auto>;

namespace LinBox {

    /**
     * Singularity of the system.
     * Only meaningful if the matrix is square.
     */
    enum class Singularity {
        Unknown,     //!< We don't know yet, or the matrix is not square.
        Singular,    //!< The matrix is non-invertible.
        NonSingular, //!< The matrix is invertible.
    };

    /**
     * For integer-based methods that evaluate multiple
     * times the system at different moduli,
     * decides how to dispatch each sub-computations.
     */
    enum class Dispatch {
        Auto,        //!< Let implementation decide what to use.
        Sequential,  //!< All sub-computations are done sequentially.
        Smp,         //!< Use symmetric multiprocessing (Paladin) to do sub-computations.
        Distributed, //!< Use MPI to distribute sub-computations accross nodes.
        Combined,    //!< Use MPI then Paladin on each node.
    };

    /**
     * For Dixon method, which solution type to get when the system is singular.
     */
    enum class SingularSolutionType {
        Determinist,    //!< The solution should be the easiest to compute and always the same.
        Random,         //!< The solution should be random and different at each call.
        Diophantine,    //!< The solution is given over the integers.
    };

    /**
     * Holds everything a method needs to know about the problem.
     */
    struct MethodBase {
        // Generic system information.
        Singularity singularity = Singularity::Unknown;
        size_t rank = 0; //!< Rank of the system. 0 means unknown.

        // For Integer-based systems.
        Dispatch dispatch = Dispatch::Auto;
        Communicator* pCommunicator = nullptr;
        bool master() const { return (pCommunicator == nullptr) || pCommunicator->master(); }

        // For Dixon method.
        SingularSolutionType singularSolutionType = SingularSolutionType::Determinist;

        // For random-based systems.
        size_t trialsBeforeThrowing = 100;        //!< Maximum number of trials before giving up.
        bool findInconsistencyCertificate = true; //!< Whether the solver should attempt to find a certificate of inconsistency if
                                                  //!  it suspects the system to be inconsistent.

        // For block-based methods.
        size_t blockingFactor = 16; //!< @fixme CHECK Size of blocks.
    };

    /**
     * Define which method to use when working on a system.
     */
    struct MethodWIP {
        // MethodWIP::Auto chooses between all following methods,
        // given the types and the dimensions of the problem.
        DEFINE_METHOD(Auto, void);

        //
        // Elimination methods
        //

        // MethodWIP::Elimination forwards to MethodWIP::DenseElimination or MethodWIP::SparseElimination.
        DEFINE_METHOD(Elimination, void);

        // MethodWIP::DenseElimination converts the matrix to dense one if necessary
        // and is using a pivot to solve the system.
        // (Gauss algorithm - Dumas, Giorgi, Pernet ISSAC 2004)
        DEFINE_METHOD(DenseElimination, void);

        // MethodWIP::SparseElimination converts the matrix to sparse one if necessary
        // and is using a pivot to solve the system.
        // (Sparse gauss algorithm - Dumas, Villard CASC 2002)
        DEFINE_METHOD(SparseElimination, void);

        //
        // Integer-based methods
        //

        // MethodWIP::Dixon uses Dixon's p-adic lifting.
        // (Numerische Mathematik - Dixon 1982)
        DEFINE_METHOD(Dixon, RingCategories::IntegerTag);

        // MethodWIP::Cra uses the chinese remainder algorithm
        // to solve the problem on multiple modular domains,
        // and finally reconstruct the solution.
        DEFINE_COMPOUND_METHOD(Cra, RingCategories::IntegerTag);

        // MethodWIP::NumericSymbolicOverlap uses Youse's overlap-based numeric/symbolic iteration.
        // (Numeric symbolic overlap iteration - Saunders, Wood, Youse ISSAC 2011)
        DEFINE_METHOD(NumericSymbolicOverlap, RingCategories::IntegerTag);

        // MethodWIP::NumericSymbolicNorm uses Wan's (older) norm-based numeric/symbolic iteration.
        // (Numeric symbolic norm iteration - Saunders, Wan ISSAC 2004)
        DEFINE_METHOD(NumericSymbolicNorm, RingCategories::IntegerTag);

        //
        // Blackbox methods
        //

        // MethodsWIP::Blackbox uses the best available Blackbox method below.
        DEFINE_METHOD(Blackbox, void);

        // MethodWIP::Wiedemann uses a blackbox algorithm
        // that projects random vectors on the both sides of the matrix
        // to find out the minpoly.
        // (IEEE Transactions on Information Theory - Wiedemann 1986)
        DEFINE_METHOD(Wiedemann, void);

        // (On Randomized Lanczos Algorithms - Kaltofel Eberly ISAAC 1997)
        DEFINE_METHOD(Lanczos, void);

        // (Linear Algebra and its Applications - Coppersmith 1993)
        DEFINE_METHOD(BlockLanczos, void);

        //
        // @deprecated Blackbox methods, kept but not tested.
        //

        // (Mathematics of Computations - Coppersmith 1994)
        DEFINE_METHOD(BlockWiedemann, void);

        // (Mathematics of Computations - Coppersmith 1994)
        DEFINE_METHOD(Coppersmith, void);
    };
}

#undef DEFINE_METHOD
