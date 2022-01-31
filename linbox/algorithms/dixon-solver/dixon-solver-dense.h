/*
 * Copyright (C) LinBox Team
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

#pragma once

#include "../rational-solver.h"

namespace LinBox {
#ifdef RSTIMING
    class DixonTimer {
    public:
        mutable Timer ttSetup, ttRecon, ttGetDigit, ttGetDigitConvert, ttRingApply, ttRingOther;
        mutable int rec_elt;
        void clear()
        {
            ttSetup.clear();
            ttRecon.clear();
            ttGetDigit.clear();
            ttGetDigitConvert.clear();
            ttRingOther.clear();
            ttRingApply.clear();
            rec_elt = 0;
        }

        template <class RR, class LC>
        void update(RR& rr, LC& lc)
        {
            ttSetup += lc.ttSetup;
            ttRecon += rr.ttRecon;
            rec_elt += rr._num_rec;
            ttGetDigit += lc.ttGetDigit;
            ttGetDigitConvert += lc.ttGetDigitConvert;
            ttRingOther += lc.ttRingOther;
            ttRingApply += lc.ttRingApply;
        }
    };
#endif

    /** \brief partial specialization of p-adic based solver with Dixon algorithm.
     *
     *   See the following reference for details on this algorithm:
     *   @bib
     *  - John D. Dixon <i>Exact Solution of linear equations using p-adic
     *  expansions</i>. Numerische Mathematik, volume 40, pages 137-141,
     *  1982.
     *
     */
    template <class Ring, class Field, class RandomPrime>
    class DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination> {

    public:
        typedef Ring RingType;
        typedef typename Ring::Element Integer;
        typedef typename Field::Element Element;
        typedef typename RandomPrime::Prime_Type Prime;

        // polymorphic 'certificate' generated when level >= SL_CERTIFIED
        // certificate of inconsistency when any solving routine returns SS_INCONSISTENT
        // certificate of minimal denominator when findRandomSolutionAndCertificate is called &
        // return is SS_OK
        mutable VectorFraction<Ring> lastCertificate;

        // next 2 fields generated only by findRandomSolutionAndCertificate, when return is SS_OK
        mutable Integer lastZBNumer;            // filled in if level >= SL_CERTIFIED
        mutable Integer lastCertifiedDenFactor; // filled in if level >= SL_LASVEGAS
        // note: lastCertificate * b = lastZBNumer / lastCertifiedDenFactor, in lowest form

    protected:
        mutable RandomPrime _genprime;
        mutable Prime _prime;
        Ring _ring;
        Field _field;

        BlasMatrixDomain<Field> _bmdf;

#ifdef RSTIMING
        mutable Timer tSetup, ttSetup, tFastInvert,
            ttFastInvert,                                               // only done in deterministic or inconsistent
            tCheckConsistency, ttCheckConsistency,                      // includes lifting the certificate
            tMakeConditioner, ttMakeConditioner, tInvertBP, ttInvertBP, // only done in random
            tCheckAnswer, ttCheckAnswer, tCertSetup,
            ttCertSetup, // remaining 3 only done when makeMinDenomCert = true
            tCertMaking, ttCertMaking,

            tNonsingularSetup, ttNonsingularSetup, tNonsingularInv, ttNonsingularInv,

            totalTimer;

        mutable DixonTimer ttConsistencySolve, ttSystemSolve, ttCertSolve, ttNonsingularSolve;
#endif

    public:
        /** Constructor
         * @param r   a Ring, set by default
         * @param rp  a RandomPrime generator, set by default
         */
        DixonSolver(const Ring& r = Ring(), const RandomPrime& rp = RandomPrime())
            : lastCertificate(r, 0)
            , _genprime(rp)
            , _ring(r)
        {
            _genprime.setBits(FieldTraits<Field>::bestBitSize());
            _prime = *_genprime;
            ++_genprime;
#ifdef RSTIMING
            clearTimers();
#endif
        }

        /** Constructor, trying the prime p first
         * @param p a Prime
         * @param r a Ring, set by default
         * @param rp a RandomPrime generator, set by default
         */
        DixonSolver(const Prime& p, const Ring& r = Ring(), const RandomPrime& rp = RandomPrime())
            : lastCertificate(r, 0)
            , _genprime(rp)
            , _prime(p)
            , _ring(r)
        {
            _genprime.setBits(FieldTraits<Field>::bestBitSize());
#ifdef RSTIMING
            clearTimers();
#endif
        }

        /** Solve a linear system \c Ax=b over quotient field of a ring.
         *
         * @param num Vector of numerators of the solution
         * @param den  The common denominator. 1/den * num is the rational solution of \c Ax=b.
         * @param A        Matrix of linear system
         * @param b        Right-hand side of system
         * @param s
         * @param maxPrimes maximum number of moduli to try
         * @param level    level of certification to be used
         *
         * @return status of solution. if \c (return != SS_FAILED), and \c (level >= SL_LASVEGAS),
         * solution is guaranteed correct. \c  SS_FAILED - all primes used were bad \c SS_OK -
         * solution found. \c  SS_INCONSISTENT - system appreared inconsistent. certificate is in \p
         * lastCertificate if \c (level >= SL_CERTIFIED)
         */
        template <class IMatrix, class Vector1, class Vector2>
        SolverReturnStatus solve(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, const bool s = false,
                                 const int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT);

        // @fixme Can we remove that?
        /** overload so that the bool 'oldMatrix' argument is not accidentally set to true */
        template <class IMatrix, class Vector1, class Vector2>
        SolverReturnStatus solve(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, const int maxPrimes,
                                 const SolverLevel level = SL_DEFAULT)
        {
            return solve(num, den, A, b, false, maxPrimes, level);
        }

        /** Solve a nonsingular, square linear system \c Ax=b over quotient field of a ring.
         *
         * @param num       Vector of numerators of the solution
         * @param den       The common denominator. <code>1/den * num</code> is the rational
         * solution of <code>Ax = b</code>
         * @param A         Matrix of linear system (it must be square)
         * @param b         Right-hand side of system
         * @param s         unused
         * @param maxPrimes maximum number of moduli to try
         *
         * @return status of solution :
         *   - \c SS_FAILED   all primes used were bad;
         *   - \c SS_OK       solution found, guaranteed correct;
         *   - \c SS_SINGULAR system appreared singular mod all primes.
         *   .
         */
        template <class IMatrix, class Vector1, class Vector2>
        SolverReturnStatus solveNonsingular(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, bool s = false,
                                            int maxPrimes = DEFAULT_MAXPRIMES);

        /** Solve a general rectangular linear system \c Ax=b over quotient field of a ring.
         *  If A is known to be square and nonsingular, calling solveNonsingular is more efficient.
         *
         * @param num       Vector of numerators of the solution
         * @param den       The common denominator. <code>1/den * num</code> is the rational
         * solution of <code>Ax = b</code>
         * @param A         Matrix of linear system
         * @param b         Right-hand side of system
         * @param maxPrimes maximum number of moduli to try
         * @param level     level of certification to be used
         *
         * @return status of solution. if <code>(return != SS_FAILED)</code>, and <code>(level >=
         * SL_LASVEGAS)</code>, solution is guaranteed correct.
         *   - \c SS_FAILED        all primes used were bad
         *   - \c SS_OK            solution found.
         *   - \c SS_INCONSISTENT  system appeared inconsistent. certificate is in \p
         * lastCertificate if <code>(level >= SL_CERTIFIED)</code>
         */
        template <class IMatrix, class Vector1, class Vector2>
        SolverReturnStatus solveSingular(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b,
                                         int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT);

        /** Find a random solution of the general linear system  \c Ax=b over quotient field of a
         * ring.
         *
         * @param num   Vector of numerators of the solution
         * @param den   The common denominator. <code>1/den * num</code> is the rational solution of
         * <code>Ax = b</code>.
         * @param A         Matrix of linear system
         * @param b         Right-hand side of system
         * @param maxPrimes maximum number of moduli to try
         * @param level     level of certification to be used
         *
         * @return status of solution. if <code>(return != SS_FAILED)</code>, and <code>(level >=
         * SL_LASVEGAS)</code>, solution is guaranteed correct.
         *  - \c SS_FAILED  all primes used were bad
         *  - \c SS_OK  solution found.
         *  - \c SS_INCONSISTENT  system appreared inconsistent. certificate is in lastCertificate
         * if <code>(level >= SL_CERTIFIED)</code>
         */
        template <class IMatrix, class Vector1, class Vector2>
        SolverReturnStatus findRandomSolution(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b,
                                              int maxPrimes = DEFAULT_MAXPRIMES, const SolverLevel level = SL_DEFAULT);

        /** Big solving routine to perform random solving and certificate generation.
         * Same arguments and return as findRandomSolution, except
         *
         * @param num  Vector of numerators of the solution
         * @param den  The common denominator. <code>1/den * num</code> is the rational solution of
         * <code>Ax = b</code>
         * @param A
         * @param b
         * @param randomSolution  parameter to determine whether to randomize or not (since
         * solveSingular calls this function as well)
         * @param makeMinDenomCert  determines whether a partial certificate for the minimal
         * denominator of a rational solution is made
         * @param maxPrimes
         * @param level
         *
         * When <code>(randomSolution == true && makeMinDenomCert == true)</code>,
         *  - If <code>(level == SL_MONTECARLO)</code> this function has the same effect as calling
         * findRandomSolution.
         *  - If <code>(level >= SL_LASVEGAS && return == SS_OK)</code>, \c lastCertifiedDenFactor
         * contains a certified factor of the min-solution's denominator.
         *  - If <code>(level >= SL_CERTIFIED && return == SS_OK)</code>, \c lastZBNumer and \c
         * lastCertificate are updated as well.
         *
         */
        template <class IMatrix, class Vector1, class Vector2>
        SolverReturnStatus monolithicSolve(Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, Method::Dixon method);

        Ring getRing() { return _ring; }

        void chooseNewPrime()
        {
            ++_genprime;
            _prime = *_genprime;
        }

#ifdef RSTIMING
        void clearTimers()
        {
            ttSetup.clear();
            ttFastInvert.clear();
            ttCheckConsistency.clear();
            ttMakeConditioner.clear();
            ttInvertBP.clear();
            ttCheckAnswer.clear();
            ttCertSetup.clear();
            ttCertMaking.clear();
            ttNonsingularSetup.clear();
            ttNonsingularInv.clear();

            ttConsistencySolve.clear();
            ttSystemSolve.clear();
            ttCertSolve.clear();
            ttNonsingularSolve.clear();
        }

    public:
        inline std::ostream& printTime(const Timer& timer, const char* title, std::ostream& os, const char* pref = "")
        {
            if (&timer != &totalTimer) totalTimer += timer;
            if (timer.count() > 0) {
                os << pref << title;
                for (int i = strlen(title) + strlen(pref); i < 28; i++) os << ' ';
                return os << timer << std::endl;
            }
            else
                return os;
        }

        inline std::ostream& printDixonTime(const DixonTimer& timer, const char* title, std::ostream& os)
        {
            if (timer.ttSetup.count() > 0) {
                printTime(timer.ttSetup, "Setup", os, title);
                printTime(timer.ttGetDigit, "Field Apply", os, title);
                printTime(timer.ttGetDigitConvert, "Ring-Field-Ring Convert", os, title);
                printTime(timer.ttRingApply, "Ring Apply", os, title);
                printTime(timer.ttRingOther, "Ring Other", os, title);
                printTime(timer.ttRecon, "Reconstruction", os, title);
                os << " number of elt recontructed: " << timer.rec_elt << std::endl;
            }
            return os;
        }

        std::ostream& reportTimes(std::ostream& os)
        {
            totalTimer.clear();
            printTime(ttNonsingularSetup, "NonsingularSetup", os);
            printTime(ttNonsingularInv, "NonsingularInv", os);
            printDixonTime(ttNonsingularSolve, "NS ", os);
            printTime(ttSetup, "Setup", os);
            printTime(ttFastInvert, "FastInvert", os);
            printTime(ttCheckConsistency, "CheckConsistency", os);
            printDixonTime(ttConsistencySolve, "INC ", os);
            printTime(ttMakeConditioner, "MakeConditioner", os);
            printTime(ttInvertBP, "InvertBP", os);
            printDixonTime(ttSystemSolve, "SYS ", os);
            printTime(ttCheckAnswer, "CheckAnswer", os);
            printTime(ttCertSetup, "CertSetup", os);
            printDixonTime(ttCertSolve, "CER ", os);
            printTime(ttCertMaking, "CertMaking", os);
            printTime(totalTimer, "TOTAL", os);
            return os;
        }
#endif

    private:
        /// Internal usage
        template <class TAS>
        SolverReturnStatus solveApparentlyInconsistent(const BlasMatrix<Ring>& A, TAS& tas, BlasMatrix<Field>* Atp_minor_inv,
                                                       size_t rank, const MethodBase& method);

        /// Internal usage
        /// @note P seems to be a preconditioner: a random matrix filled with 0 or 1
        /// @note Please note that Atp_minor_inv content *might* be changed after this function
        /// call.
        template <class TAS>
        void makeConditioner(BlasMatrix<Ring>& A_minor, BlasMatrix<Field>*& Ap_minor_inv, BlasMatrix<Ring>*& B,
                             BlasMatrix<Ring>*& P, const BlasMatrix<Ring>& A, TAS& tas, BlasMatrix<Field>* Atp_minor_inv,
                             size_t rank, const MethodBase& method);

        template <class Vector1, class Vector2, class TAS>
        void certifyMinimalDenominator(const BlasMatrix<Ring>& A, const Vector2& b, const TAS& tas, const BlasMatrix<Ring>& B,
                                       BlasMatrix<Ring>& A_minor, BlasMatrix<Field>& Ap_minor_inv, size_t rank);
    };
}

#include "./dixon-solver-dense.inl"
