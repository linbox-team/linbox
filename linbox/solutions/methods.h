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

#include <linbox/field/field-traits.h>
#include <linbox/matrix/dense-matrix.h> // Only for useBlackboxMethod
#include <linbox/solutions/constants.h>
#include <linbox/util/mpicpp.h>
#include <string>

/**
 * These macros are used to define methods.
 *
 * The CategoryTag is used to differenciate what the method is used for.
 * For instance, Method::Dixon has its CategoryTag set to RingCategories::IntegerTag,
 * expressing it is used only for non-modular computations.
 * Setting CategoryTag to void means that we don't care.
 *
 * A compound method is something like Method::CRA<Method::Auto>
 * telling the method CRA which method to use.
 */
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
    using _MethodName##Auto = _MethodName<Method::Auto>;

namespace LinBox {

    // Used to decide which method to use when using Method::Auto on a Blackbox or Sparse matrix.
    template <class Matrix>
    bool useBlackboxMethod(const Matrix& A)
    {
        return (A.coldim() > LINBOX_USE_BLACKBOX_THRESHOLD)
            && (A.rowdim() > LINBOX_USE_BLACKBOX_THRESHOLD);
    }

    template <class Field>
    bool useBlackboxMethod(const LinBox::DenseMatrix<Field>& A)
    {
        return false;
    }

    /**
     * Rank of the system, if known.
     */
    namespace Rank {
        enum Value : int16_t {
            Unknown = -1,
        };
    }

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
        SMP,         //!< Use symmetric multiprocessing (Paladin) to do sub-computations.
        Distributed, //!< Use MPI to distribute sub-computations accross nodes.
        Combined,    //!< Use MPI then Paladin on each node.
    };

    /**
     * For Dixon method, which solution type to get when the system is singular.
     */
    enum class SingularSolutionType {
        Deterministic, //!< The solution should be the easiest to compute and always the same.
        Random,        //!< The solution should be random and different at each call.
        Diophantine,   //!< The solution is given over the integers.
    };

    /**
     * Preconditioner to ensure generic rank profile.
     */
    enum class Preconditioner {
        None,                      //!< Do not use any preconditioner.
        Butterfly,                 //!< Use a butterfly network, see @ref Butterfly.
        Sparse,                    //!< Use a sparse preconditioner, c.f. (Mulders 2000).
        Toeplitz,                  //!< Use a Toeplitz preconditioner, c.f. (Kaltofen and Saunders 1991).
        Symmetrize,                //!< Use At A (Eberly and Kaltofen 1997).
        PartialDiagonal,           //!< Use A D, where D is a random non-singular diagonal matrix (Eberly and Kaltofen 1997).
        PartialDiagonalSymmetrize, //!< Use At D A (Eberly and Kaltofen 1997).
        FullDiagonal,              //!< Use D1 At D2 A D1 (Eberly and Kaltofen 1997).
        Dense,                     //!< Multiply (@fixme or add?) by a random dense matrix (used by Dixon).
    };

    /**
     * Flags decribing the shape of the matrix.
     *
     * @note The namespace is here to avoid collisions.
     */
    namespace Shape {
        enum Value : uint16_t {
            Unknown = 0x00,
            Symmetric = 0x01,       //!< Matrix has its main diagonal as a reflection axis.
            Diagonal = 0x02,        //!< Outside main diagonal are zeroes.
            Toeplitz = 0x04,        //!< Diagonals are constant.
            Hankel = 0x08,          //!< Anti-diagonals are constant.
            Unimodular = 0x10,      //!< Square integer matrix having determinant +1 or −1.
            UpperTriangular = 0x20, //!< Only upper right part of the matrix is non-zero.
            LowerTriangular = 0x40, //!< Only lower left part of the matrix is non-zero.
        };
    };
    struct ShapeFlags {
        uint16_t flags;
        ShapeFlags(const uint16_t& f = Shape::Unknown) : flags(f) {}

        // @note We overload operator== so that one can use
        // method.shapeFlags == Shape::Symmetric seemlessly.
        inline ShapeFlags& operator=(Shape::Value shape)
        {
            (shape == Shape::Unknown) ? (flags = Shape::Unknown) : (flags |= shape);
            return *this;
        }
        inline bool operator==(Shape::Value shape) const
        {
            return (shape == Shape::Unknown) ? (flags == Shape::Unknown) : (flags & shape);
        }
        inline bool operator!=(Shape::Value shape) const { return !(flags == shape); }
    };

    /**
     * Pivoting strategy for elimination-based methods.
     */
    enum class PivotStrategy {
        None,
        Linear,
    };

    /**
     * Holds everything a method needs to know about the problem.
     *
     * @note This "put everything in it" design is used so that
     * Method::Auto can specify any parameter of whatever method
     * is going to be used.
     * Fact is this structure is never copied between method switches.
     * It does not matter if it is too big, but for clarity, we try
     * to regroup elements used depending on the method.
     */
    struct MethodBase {
        MethodBase() = default;
        MethodBase(Singularity _singularity) : singularity(_singularity) {}
        MethodBase(Shape::Value _shapeFlag) : shapeFlags(_shapeFlag) {}
        MethodBase(ShapeFlags _shapeFlags) : shapeFlags(_shapeFlags) {}
        MethodBase(Preconditioner _preconditioner) : preconditioner(_preconditioner) {}
        MethodBase(Dispatch _dispatch) : dispatch(_dispatch) {}
        MethodBase(Communicator* _pCommunicator) : pCommunicator(_pCommunicator) {}
        MethodBase(PivotStrategy _pivotStrategy) : pivotStrategy(_pivotStrategy) {}
        MethodBase(SingularSolutionType _singularSolutionType) : singularSolutionType(_singularSolutionType) {}

        // ----- Generic system information.
        Singularity singularity = Singularity::Unknown;
        Rank::Value rank = Rank::Unknown; //!< Rank of the system. -1 means unknown.
        ShapeFlags shapeFlags;            //!< Shape of the system.

        // ----- Generic solve options.
        Preconditioner preconditioner = Preconditioner::None;
        bool checkResult = false; //!< Ensure that solving worked by checking Ax = b (might not be implemented by all methods).

        // ----- For Integer-based systems.
        Dispatch dispatch = Dispatch::Auto;
        Communicator* pCommunicator = nullptr;
        bool master() const { return (pCommunicator == nullptr) || pCommunicator->master(); }

        // ----- For Elimination-based methods.
        PivotStrategy pivotStrategy = PivotStrategy::Linear;

        // ----- For Dixon method.
        // @fixme SingularSolutionType::Deterministic fails with Dense Dixon
        SingularSolutionType singularSolutionType = SingularSolutionType::Random;
        bool certifyMinimalDenominator = false; //!< Whether the solver should try to find a certificate
                                                //!  that the provided denominator is minimal.

        // ----- For random-based systems.
        size_t trialsBeforeFailure = LINBOX_DEFAULT_TRIALS_BEFORE_FAILURE; //!< Maximum number of trials before giving up.
        bool certifyInconsistency = false; //!< Whether the solver should attempt to find a certificate of inconsistency if
                                           //!  it suspects the system to be inconsistent.

        // ----- For block-based methods.
        size_t blockingFactor = LINBOX_DEFAULT_BLOCKING_FACTOR; //!< Size of blocks.

        // ----- For Wiedemann (Berlekamp Massey) methods.
        size_t earlyTerminationThreshold = LINBOX_DEFAULT_EARLY_TERMINATION_THRESHOLD;
    };

    /**
     * Define which method to use when working on a system.
     */
    struct Method {
        // Method::Auto chooses between all following methods,
        // given the types and the dimensions of the problem.
        DEFINE_METHOD(Auto, void);

        //
        // Elimination methods
        //

        // Method::Elimination forwards to Method::DenseElimination or Method::SparseElimination.
        DEFINE_METHOD(Elimination, void);

        // Method::DenseElimination converts the matrix to dense one if necessary
        // and is using a pivot to solve the system.
        // (Gauss algorithm - Dumas, Giorgi, Pernet ISSAC 2004)
        DEFINE_METHOD(DenseElimination, void);

        // Method::SparseElimination converts the matrix to sparse one if necessary
        // and is using a pivot to solve the system.
        // (Sparse gauss algorithm - Dumas, Villard CASC 2002)
        DEFINE_METHOD(SparseElimination, void);

        //
        // Integer-based methods
        //

        // Method::Dixon uses Dixon's p-adic lifting.
        // (Numerische Mathematik - Dixon 1982)
        DEFINE_METHOD(Dixon, RingCategories::IntegerTag);

        // Method::ChineseRemainder uses the chinese remainder algorithm
        // to solve the problem on multiple modular domains,
        // and finally reconstruct the solution.
        DEFINE_COMPOUND_METHOD(CRA, RingCategories::IntegerTag);

        // Method::SymbolicNumericOverlap uses Youse's overlap-based numeric/symbolic iteration.
        // (Numeric symbolic overlap iteration - Saunders, Wood, Youse ISSAC 2011)
        DEFINE_METHOD(SymbolicNumericOverlap, RingCategories::IntegerTag);

        // Method::SymbolicNumericNorm uses Wan's (older) norm-based numeric/symbolic iteration.
        // (Numeric symbolic norm iteration - Saunders, Wan ISSAC 2004)
        DEFINE_METHOD(SymbolicNumericNorm, RingCategories::IntegerTag);

        //
        // Blackbox methods
        //

        // MethodsWIP::Blackbox uses the best available Blackbox method below.
        DEFINE_METHOD(Blackbox, void);

        // Method::Wiedemann uses a blackbox algorithm
        // that projects random vectors on the both sides of the matrix
        // to find out the minpoly.
        // (IEEE Transactions on Information Theory - Wiedemann 1986)
        DEFINE_METHOD(Wiedemann, void);

        // (On Randomized Lanczos Algorithms - Kaltofel Eberly ISAAC 1997)
        DEFINE_METHOD(Lanczos, void);

        // (Linear Algebra and its Applications - Coppersmith 1993)
        DEFINE_METHOD(BlockLanczos, void);

        //
        // Minpoly methods
        //

        // If the field is too short...
        // @todo This might not be needed as a method per se.
        // An method attribute useExtension might do the work here.
        // But it can be fixed while working on https://github.com/linbox-team/linbox/issues/190
        DEFINE_METHOD(WiedemannExtension, void);

        //
        // @deprecated Blackbox methods, kept but not tested.
        //

        // (Mathematics of Computations - Coppersmith 1994)
        DEFINE_METHOD(BlockWiedemann, void);

        // (Mathematics of Computations - Coppersmith 1994)
        DEFINE_METHOD(Coppersmith, void);

        // Use Krylov sequence to solve.
        DEFINE_METHOD(BlockHankel, RingCategories::IntegerTag);
    };
}

#undef DEFINE_METHOD
