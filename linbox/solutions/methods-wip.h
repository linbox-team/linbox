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

#define DEFINE_METHOD_CONTENT(MethodName, MethodNameString)                                                                      \
    static std::string name() { return std::string("Method::") + MethodNameString; }                                             \
    MethodName() = default;                                                                                                      \
    MethodName(const MethodName&) = default;                                                                                     \
    MethodName(const MethodBase& methodBase)                                                                                     \
        : MethodBase(methodBase)                                                                                                 \
    {                                                                                                                            \
    }
#define DEFINE_METHOD(MethodName)                                                                                                \
    struct MethodName : public MethodBase {                                                                                      \
        DEFINE_METHOD_CONTENT(MethodName, #MethodName)                                                                           \
    };

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
     * For Dixon method, which solution type to use.
     * @fixme Replace with HeuristicTag and such?
     */
    enum class SolutionType {
        Determinist,
        Random,
        Diophantine,
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
        SolutionType solutionType = SolutionType::Determinist;

        // For random-based systems.
        size_t trialsBeforeThrowing = 100;        //!< Maximum number of trials before giving up.
        bool findInconsistencyCertificate = true; //!< Whether the solver should attempt to find a certificate of inconsistency if
                                                  //!  it suspects the system to be inconsistent.

        // For block-based methods.
        size_t blockingFactor = 16; //!< @fixme CHECK Size of blocks.
    };

    //
    // Integer-based methods.
    //

    template <class IterationMethod>
    struct CraCustomMethod : public MethodBase {
        IterationMethod iterationMethod;

        DEFINE_METHOD_CONTENT(CraCustomMethod<IterationMethod>, "Cra<" + IterationMethod::name() + ">");
    };

    template <class IterationMethod>
    struct DixonCustomMethod : public MethodBase {
        IterationMethod iterationMethod;

        DEFINE_METHOD_CONTENT(DixonCustomMethod<IterationMethod>, "Dixon<" + IterationMethod::name() + ">");
    };

    /**
     * Define which method to use when working on a system.
     */
    struct MethodWIP {
        DEFINE_METHOD(Auto);

        // Elimination methods
        DEFINE_METHOD(Elimination);
        DEFINE_METHOD(DenseElimination);
        DEFINE_METHOD(SparseElimination);

        // Integer-based methods
        template <class IterationMethod>
        using DixonCustom = DixonCustomMethod<IterationMethod>;
        using Dixon = DixonCustom<MethodWIP::Auto>;
        template <class IterationMethod>
        using CraCustom = CraCustomMethod<IterationMethod>;
        using Cra = CraCustom<MethodWIP::Auto>;

        // Blackbox methods
        DEFINE_METHOD(Blackbox);
        DEFINE_METHOD(Wiedemann);
        DEFINE_METHOD(BlockWiedemann);
        DEFINE_METHOD(Coppersmith);
    };
}

#undef DEFINE_METHOD
