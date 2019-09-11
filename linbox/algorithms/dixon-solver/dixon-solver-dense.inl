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

#include "linbox/linbox-config.h"
#include "linbox/util/debug.h"

#include "linbox/algorithms/lifting-container.h"
#include "linbox/algorithms/matrix-inverse.h"
#include "linbox/algorithms/rational-reconstruction.h"

namespace LinBox {

    template <class Ring, class Field, class RandomPrime>
    template <class IMatrix, class Vector1, class Vector2>
    SolverReturnStatus DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination>::solve(Vector1& num, Integer& den,
                                                                                              const IMatrix& A, const Vector2& b,
                                                                                              const bool old, int maxP,
                                                                                              const SolverLevel level)
    {
        SolverReturnStatus status;
        int maxPrimes = maxP;
        while (maxPrimes > 0) {
            auto nonSingularResult = (A.rowdim() == A.coldim()) ? solveNonsingular(num, den, A, b, old, maxPrimes) : SS_SINGULAR;
            switch (nonSingularResult) {
            case SS_OK: return SS_OK; break;

            case SS_SINGULAR:
                status = solveSingular(num, den, A, b, maxPrimes, level);
                if (status != SS_FAILED) return status;
                break;

            case SS_FAILED: break;

            default: throw LinboxError("Bad return value from solveNonsingular");
            }
            maxPrimes--;
            if (maxPrimes > 0) chooseNewPrime();
        }
        return SS_FAILED;
    }

    template <class Ring, class Field, class RandomPrime>
    template <class IMatrix, class Vector1, class Vector2>
    SolverReturnStatus DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination>::solveNonsingular(
        Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, bool oldMatrix, int maxPrimes)
    {

        int trials = 0, notfr;

        // history sensitive data for optimal reason
        // static const IMatrix* IMP;

        BlasMatrix<Field>* FMP = NULL;
        Field* F = NULL;

        do {
#ifdef RSTIMING
            tNonsingularSetup.start();
#endif
            // typedef typename Field::Element Element;
            // typedef typename Ring::Element Integer;

            // checking size of system
            linbox_check(A.rowdim() == A.coldim());
            linbox_check(A.rowdim() == b.size());

            LinBox::integer tmp;

            // if input matrix A is different one.
            if (!oldMatrix) {
                if (trials == maxPrimes) return SS_SINGULAR;
                if (trials != 0) chooseNewPrime();
                ++trials;

                // Could delete a non allocated matrix -> segfault
                if (FMP != NULL) delete FMP;

                // IMP = &A;

                if (F != NULL) delete F;

                F = new Field(_prime);
                FMP = new BlasMatrix<Field>(*F, A.rowdim(), A.coldim());
                MatrixHom::map(*FMP, A); // use MatrixHom to reduce matrix PG 2005-06-16

                BlasMatrix<Field>* invA = new BlasMatrix<Field>(*F, A.rowdim(), A.coldim());
                BlasMatrixDomain<Field> BMDF(*F);
#ifdef RSTIMING
                tNonsingularSetup.stop();
                ttNonsingularSetup += tNonsingularSetup;
                tNonsingularInv.start();
#endif
                assert(FMP != NULL);
                BMDF.invin(*invA, *FMP, notfr); // notfr <- nullity
                delete FMP;
                FMP = invA;

#ifdef RSTIMING
                tNonsingularInv.stop();
                ttNonsingularInv += tNonsingularInv;
#endif
            }
            else {
#ifdef RSTIMING
                tNonsingularSetup.stop();
                ttNonsingularSetup += tNonsingularSetup;
#endif
                notfr = 0;
            }
        } while (notfr);

        typedef DixonLiftingContainer<Ring, Field, IMatrix, BlasMatrix<Field>> LiftingContainer;
        LiftingContainer lc(_ring, *F, A, *FMP, b, _prime);
        RationalReconstruction<LiftingContainer> re(lc);
        if (!re.getRational(num, den, 0)) {
            delete FMP;
            return SS_FAILED;
        }
#ifdef RSTIMING
        ttNonsingularSolve.update(re, lc);
#endif
        if (F != NULL) delete F;
        if (FMP != NULL) delete FMP;
        return SS_OK;
    }

    template <class Ring, class Field, class RandomPrime>
    template <class IMatrix, class Vector1, class Vector2>
    SolverReturnStatus DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination>::solveSingular(
        Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, int maxPrimes, const SolverLevel level)
    {
        Method::Dixon m;
        m.certifyMinimalDenominator = false;
        m.certifyInconsistency = false;
        m.trialsBeforeFailure = maxPrimes;
        return monolithicSolve(num, den, A, b, m);
    }

    template <class Ring, class Field, class RandomPrime>
    template <class IMatrix, class Vector1, class Vector2>
    SolverReturnStatus DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination>::findRandomSolution(
        Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, int maxPrimes, const SolverLevel level)
    {
        Method::Dixon m;
        m.singularSolutionType = SingularSolutionType::Random;
        m.certifyMinimalDenominator = false;
        m.certifyInconsistency = true;
        m.trialsBeforeFailure = maxPrimes;
        return monolithicSolve(num, den, A, b, m);
    }

    // TAS stands for Transpose Augmented System (A|b)t
    // this provides a factorization (A|b) = Pt . Lt . Ut . Qt
    // such that
    // - Q . (A|b) . P   has nonzero principal minors up to TAS.rank()
    // - P permutes b to the (TAS.rank())th column of A iff the system is inconsistent mod p

    template <class Field>
    class TransposeAugmentedSystem {
    public:
        BlasMatrix<Field>* factors = nullptr;
        PLUQMatrix<Field>* PLUQ = nullptr;

        BlasPermutation<size_t> P;
        BlasPermutation<size_t> Q;

        std::vector<size_t> srcRow;
        std::vector<size_t> srcCol;

    public:
        template <class Ring, class IMatrix, class IVector>
        TransposeAugmentedSystem(Ring& R, Field& _field, const IMatrix& A, const IVector& b)
            : P(A.coldim() + 1)
            , Q(A.rowdim())
        {
            factors = new BlasMatrix<Field>(_field, A.coldim() + 1, A.rowdim());

            BlasMatrix<Field> Ap(A, _field); // Getting into the field

            // Setting factors = [Ap|0]t
            for (size_t i = 0; i < A.rowdim(); ++i)
                for (size_t j = 0; j < A.coldim(); ++j) factors->setEntry(j, i, Ap.getEntry(i, j));

            // And then factors = [Ap|bp]t
            Integer tmpInteger;
            for (size_t i = 0; i < A.rowdim(); ++i) {
                typename Field::Element tmpElement;
                _field.init(tmpElement, R.convert(tmpInteger, b[i]));
                factors->setEntry(A.coldim(), i, tmpElement);
            }

            // Getting factors -> P L U Q
            PLUQ = new PLUQMatrix<Field>(*factors, P, Q);

            // Used to check consistency
            srcRow.resize(A.rowdim());
            srcCol.resize(A.coldim() + 1);
            FFPACK::LAPACKPerm2MathPerm(srcRow.data(), Q.getStorage().data(), srcRow.size());
            FFPACK::LAPACKPerm2MathPerm(srcCol.data(), P.getStorage().data(), srcCol.size());

            // @note srcCol/srcRow now hold permutations (indices of (A|b)t)
        }

        ~TransposeAugmentedSystem()
        {
            delete factors;
            delete PLUQ;
        }

        size_t rank() { return PLUQ->getRank(); }
    };

    // Returns:
    // - FAILED if A != 0
    // - INCONSISTENT if A == 0 and b != 0 (also sets certificate)
    // - OK if A == 0 and b == 0
    template <class Matrix, class Vector, class Ring = typename Matrix::Field>
    SolverReturnStatus certifyEmpty(const Matrix& A, const Vector& b, const MethodBase& method, Integer& certifiedDenFactor)
    {
        const Ring& R = A.field();

        // In monte carlo, we assume A is actually empty.
        bool aEmpty = true;
        if (method.certifyInconsistency) {
            MatrixDomain<Ring> MD(R);
            aEmpty = MD.isZero(A);
        }

        if (!aEmpty) {
            return SS_FAILED;
        }

        // @fixme Use VectorDomain::isZero
        for (size_t i = 0; i < b.size(); ++i) {
            if (!R.areEqual(b[i], R.zero)) {
                return SS_INCONSISTENT;
            }
        }

        // System is consistent
        if (method.certifyInconsistency) {
            R.assign(certifiedDenFactor, R.one);
        }

        return SS_OK;
    }

    // ---------------------------------------------
    // SS_FAILED means that we will need to try a new prime
    template <class Ring, class Field, class RandomPrime>
    template <class TAS>
    SolverReturnStatus DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination>::solveApparentlyInconsistent(
        const BlasMatrix<Ring>& A, TAS& tas, BlasMatrix<Field>* Atp_minor_inv, size_t rank, const MethodBase& method)
    {
        using LiftingContainer = DixonLiftingContainer<Ring, Field, BlasMatrix<Ring>, BlasMatrix<Field>>;

        if (!method.certifyInconsistency) return SS_INCONSISTENT;

        // @fixme Put these as class members!
        BlasMatrixDomain<Ring> BMDI(_ring);
        BlasApply<Ring> BAR(_ring);

#ifdef RSTIMING
        tCheckConsistency.start();
#endif

        BlasVector<Ring> zt(_ring, rank);
        for (size_t i = 0; i < rank; ++i) _ring.assign(zt[i], A.getEntry(tas.srcRow[rank], tas.srcCol[i]));

        BlasMatrix<Ring> At_minor(_ring, rank, rank);
        for (size_t i = 0; i < rank; ++i)
            for (size_t j = 0; j < rank; ++j) _ring.assign(At_minor.refEntry(j, i), A.getEntry(tas.srcRow[i], tas.srcCol[j]));

#ifdef RSTIMING
        tCheckConsistency.stop();
        ttCheckConsistency += tCheckConsistency;
#endif

        LiftingContainer lc(_ring, _field, At_minor, *Atp_minor_inv, zt, _prime);
        RationalReconstruction<LiftingContainer> re(lc);

        BlasVector<Ring> shortNum(A.field(), rank);
        Integer shortDen;

        // Dirty, but should not be called under normal circumstances
        if (!re.getRational(shortNum, shortDen, 0)) {
            return SS_FAILED;
        }

#ifdef RSTIMING
        ttConsistencySolve.update(re, lc);
        tCheckConsistency.start();
#endif

        // Build up certificate
        VectorFraction<Ring> cert(_ring, shortNum.size());
        cert.numer = shortNum;
        cert.denom = shortDen;
        cert.numer.resize(A.rowdim());
        _ring.subin(cert.numer[rank], cert.denom);
        _ring.assign(cert.denom, _ring.one);

        BMDI.mulin_left(cert.numer, tas.Q);

#ifdef DEBUG_INC
        cert.write(std::cout << "cert:") << std::endl;
#endif

        bool certifies = true; // check certificate
        BlasVector<Ring> certnumer_A(_ring, A.coldim());
        BAR.applyVTrans(certnumer_A, A, cert.numer);
        typename BlasVector<Ring>::iterator cai = certnumer_A.begin();
        for (size_t i = 0; certifies && i < A.coldim(); ++i, ++cai) {
            certifies = certifies && _ring.isZero(*cai);
        }

#ifdef RSTIMING
        tCheckConsistency.stop();
        ttCheckConsistency += tCheckConsistency;
#endif

        if (certifies) {
            if (method.certifyInconsistency) lastCertificate.copy(cert);
            return SS_INCONSISTENT;
        }
        commentator().report(Commentator::LEVEL_IMPORTANT, PARTIAL_RESULT)
            << "system is suspected to be inconsistent but it was only a bad prime" << std::endl;

        // @fixme Try new prime, is SS_FAILED the right API?
        // Analogous to u.A12 != A22 in Muld.+Storj.
        return SS_FAILED;
    }

    template <class Ring, class Field, class RandomPrime>
    template <class TAS>
    void DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination>::makeConditioner(
        BlasMatrix<Ring>& A_minor, BlasMatrix<Field>*& Ap_minor_inv, BlasMatrix<Ring>*& B, BlasMatrix<Ring>*& P,
        const BlasMatrix<Ring>& A, TAS& tas, BlasMatrix<Field>* Atp_minor_inv, size_t rank, const MethodBase& method)
    {
#ifdef RSTIMING
        tMakeConditioner.start();
#endif

        if (method.singularSolutionType != SingularSolutionType::Random) {
            // Transpose Atp_minor_inv to get Ap_minor_inv
            // @note minor inv = L1\U1
            Element _rtmp;
            Ap_minor_inv = Atp_minor_inv;
            for (size_t i = 0; i < rank; ++i) {
                for (size_t j = 0; j < i; ++j) {
                    Ap_minor_inv->getEntry(_rtmp, i, j);
                    Ap_minor_inv->setEntry(i, j, Ap_minor_inv->refEntry(j, i));
                    Ap_minor_inv->setEntry(j, i, _rtmp);
                }
            }

            // permute original entries into A_minor
            // @note A_minor = Pt A Qt
            for (size_t i = 0; i < rank; ++i)
                for (size_t j = 0; j < rank; ++j) _ring.assign(A_minor.refEntry(i, j), A.getEntry(tas.srcRow[i], tas.srcCol[j]));
#ifdef RSTIMING
            tMakeConditioner.stop();
            ttMakeConditioner += tMakeConditioner;
#endif

            if (method.certifyMinimalDenominator) {
                B = new BlasMatrix<Ring>(_ring, rank, A.coldim());
                for (size_t i = 0; i < rank; ++i)
                    for (size_t j = 0; j < A.coldim(); ++j) _ring.assign(B->refEntry(i, j), A.getEntry(tas.srcRow[i], j));
            }
            // @note B = Pt A
        }
        else { // if randomSolution == true
            P = new BlasMatrix<Ring>(_ring, A.coldim(), rank);
            B = new BlasMatrix<Ring>(_ring, rank, A.coldim());
            BlasMatrix<Field> Ap_minor(_field, rank, rank);
            Ap_minor_inv = new BlasMatrix<Field>(_field, rank, rank);

            LinBox::integer tmp2 = 0;
            size_t maxBitSize = 0;
            for (size_t i = 0; i < rank; ++i)
                for (size_t j = 0; j < A.coldim(); ++j) {
                    _ring.assign(B->refEntry(i, j), A.getEntry(tas.srcRow[i], j));
                    _ring.convert(tmp2, A.getEntry(tas.srcRow[i], j));
                    maxBitSize = std::max(maxBitSize, tmp2.bitsize());
                }
                // @note B = Pt A
#ifdef RSTIMING
            bool firstLoop = true;
#endif
            // prepare B to be preconditionned through BLAS matrix mul
            MatrixApplyDomain<Ring, BlasMatrix<Ring>> MAD(_ring, *B);
            MAD.setup(2); // @fixme Useless?

            int nullity;
            do { // O(1) loops of this preconditioner expected
#ifdef RSTIMING
                if (firstLoop)
                    firstLoop = false;
                else
                    tMakeConditioner.start();
#endif
                // compute P a n*r random matrix of entry in [0,1]
                typename BlasMatrix<Ring>::Iterator iter;
                for (iter = P->Begin(); iter != P->End(); ++iter) {
                    if (rand() % 2 == 1)
                        _ring.assign(*iter, _ring.one);
                    else
                        _ring.assign(*iter, _ring.zero);
                }

                // compute A_minor = B.P
                MAD.applyM(A_minor, *P);

                // set Ap_minor = A_minor mod p, try to compute inverse
                for (size_t i = 0; i < rank; ++i)
                    for (size_t j = 0; j < rank; ++j)
                        _field.init(Ap_minor.refEntry(i, j), _ring.convert(tmp2, A_minor.getEntry(i, j)));
#ifdef RSTIMING
                tMakeConditioner.stop();
                ttMakeConditioner += tMakeConditioner;
                tInvertBP.start();
#endif

                // @fixme Seems sad to be forced to specify these BlasMatrix<Field>& casts
                _bmdf.inv((BlasMatrix<Field>&)*Ap_minor_inv, (BlasMatrix<Field>&)Ap_minor, nullity);

#ifdef RSTIMING
                tInvertBP.stop();
                ttInvertBP += tInvertBP;
#endif
            } while (nullity > 0);
        }
    }

    template <class Ring, class Field, class RandomPrime>
    template <class Vector1, class Vector2, class TAS>
    void DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination>::certifyMinimalDenominator(
        const BlasMatrix<Ring>& A, const Vector2& b, const TAS& tas, const BlasMatrix<Ring>& B, BlasMatrix<Ring>& A_minor,
        BlasMatrix<Field>& Ap_minor_inv, size_t rank)
    {
        // To make this certificate we solve with the same matrix as to get the
        // solution, except transposed.
#ifdef RSTIMING
        tCertSetup.start();
#endif

        // @note We transpose Ap and A minors in-place because it won't be used anymore
        Integer _rtmp;
        Element _ftmp;
        for (size_t i = 0; i < rank; ++i)
            for (size_t j = 0; j < i; ++j) {
                Ap_minor_inv.getEntry(_ftmp, i, j);
                Ap_minor_inv.setEntry(i, j, Ap_minor_inv.refEntry(j, i));
                Ap_minor_inv.setEntry(j, i, _ftmp);
            }

        for (size_t i = 0; i < rank; ++i)
            for (size_t j = 0; j < i; ++j) {
                A_minor.getEntry(_rtmp, i, j);
                A_minor.setEntry(i, j, A_minor.refEntry(j, i));
                A_minor.setEntry(j, i, _rtmp);
            }

        // we then try to create a partial certificate
        // the correspondance with Algorithm MinimalSolution from Mulders/Storjohann:
        // paper | here
        // P     | TAS_Q
        // Q     | transpose of TAS_P
        // B     | *B (== TAS_Q . A,  but only top #rank rows)
        // c     | newb (== TAS_Q . b,   but only top #rank rows)
        // P     | P
        // q     | q
        // U     | {0, 1}
        // u     | u
        // z-hat | lastCertificate

        // we multiply the certificate by TAS_Qt at the end
        // so it corresponds to b instead of newb

        // q in {0, 1}^rank
        Givaro::ZRing<Integer> Z;
        BlasVector<Givaro::ZRing<Integer>> q(Z, rank);
        typename BlasVector<Givaro::ZRing<Integer>>::iterator q_iter;

        bool allzero;
        do {
            allzero = true;
            for (q_iter = q.begin(); q_iter != q.end(); ++q_iter) {
                if (rand() > RAND_MAX / 2) {
                    _ring.assign((*q_iter), _ring.one);
                    allzero = false;
                }
                else
                    (*q_iter) = _ring.zero;
            }
        } while (allzero);

#ifdef RSTIMING
        tCertSetup.stop();
        ttCertSetup += tCertSetup;
#endif

        using LiftingContainer = DixonLiftingContainer<Ring, Field, BlasMatrix<Ring>, BlasMatrix<Field>>;
        LiftingContainer lc2(_ring, _field, A_minor, Ap_minor_inv, q, _prime);

        RationalReconstruction<LiftingContainer> rere(lc2);
        Vector1 u_num(_ring, rank);
        Integer u_den;

        // Failure
        if (!rere.getRational(u_num, u_den, 0)) return;

#ifdef RSTIMING
        ttCertSolve.update(rere, lc2);
        tCertMaking.start();
#endif

        // remainder of code does   z <- denom(partial_cert . Mr) * partial_cert * Qt
        BlasApply<Ring> BAR(_ring);
        VectorFraction<Ring> u_to_vf(_ring, u_num.size());
        u_to_vf.numer = u_num;
        u_to_vf.denom = u_den;
        BlasVector<Ring> uB(_ring, A.coldim());
        BAR.applyVTrans(uB, B, u_to_vf.numer);

        Integer numergcd = _ring.zero;
        vectorGcdIn(numergcd, _ring, uB);

        // denom(partial_cert . Mr) = partial_cert_to_vf.denom / numergcd
        VectorFraction<Ring> z(_ring, b.size()); // new constructor
        u_to_vf.numer.resize(A.rowdim());

        BlasMatrixDomain<Ring> BMDI(_ring);
        BMDI.mul(z.numer, u_to_vf.numer, tas.Q);

        z.denom = numergcd;

        lastCertificate.copy(z);

        // output new certified denom factor
        Integer znumer_b, zbgcd;
        VectorDomain<Ring> VDR(_ring);
        VDR.dotprod(znumer_b, z.numer, b);
        _ring.gcd(zbgcd, znumer_b, z.denom);
        _ring.div(lastCertifiedDenFactor, z.denom, zbgcd);

        _ring.div(lastZBNumer, znumer_b, zbgcd);
#ifdef RSTIMING
        tCertMaking.stop();
        ttCertMaking += tCertMaking;
#endif
    }

    // Most solving is done by the routine below.
    // There used to be one for random and one for deterministic, but they have been merged to ease
    // with
    //  repeated code (certifying inconsistency, optimization are 2 examples)

    template <class Ring, class Field, class RandomPrime>
    template <class IMatrix, class Vector1, class Vector2>
    SolverReturnStatus DixonSolver<Ring, Field, RandomPrime, Method::DenseElimination>::monolithicSolve(
        Vector1& num, Integer& den, const IMatrix& A, const Vector2& b, Method::Dixon method)
    {
        using LiftingContainer = DixonLiftingContainer<Ring, Field, BlasMatrix<Ring>, BlasMatrix<Field>>;
        if (method.certifyMinimalDenominator && !method.certifyInconsistency) {
            method.certifyInconsistency = true;
            std::cerr << "WARNING: forcing certifyInconsistency due to certifyMinimalDenominator" << std::endl;
        }

        size_t trials = 0;
        while (trials < method.trialsBeforeFailure) {
            if (trials != 0) chooseNewPrime();
            ++trials;

#ifdef RSTIMING
            tSetup.start();
#endif
            // ----- Build Transposed Augmented System (TAS)

            // checking size of system
            linbox_check(A.rowdim() == b.size());

            _field = Field(_prime);
            _bmdf = BlasMatrixDomain<Field>(_field);

            LinBox::integer tmp;
            BlasMatrixDomain<Ring> BMDI(_ring);
            BlasApply<Ring> BAR(_ring);
            MatrixDomain<Ring> MD(_ring);

            BlasMatrix<Ring> A_check(A); // used to check answer later

            // TAS stands for Transpose Augmented System (A|b)t
            TransposeAugmentedSystem<Field> tas(_ring, _field, A, b);

#ifdef RSTIMING
            tSetup.stop();
            ttSetup += tSetup;
#endif

            // @note If permutation shows that b was needed, means b is not in the columns' span of
            // A (=> Ax=b inconsistent)
            bool appearsInconsistent = (tas.srcCol[tas.rank() - 1] == A.coldim());
            size_t rank = tas.rank() - (appearsInconsistent ? 1 : 0);

            // ----- Handle A == 0 mod p

            // Special case when A = 0, mod p. Deal with it to avoid later deadlock.
            if (rank == 0) {
                SolverReturnStatus status = certifyEmpty(A, b, method, lastCertifiedDenFactor);

                if (status == SS_FAILED) {
                    // A was empty mod p but not over Z, we try new prime.
                    continue;
                }
                else if (status == SS_INCONSISTENT) {
                    return SS_INCONSISTENT;
                }

                // It is consistent, set 0 as a valid answer
                _ring.assign(den, _ring.one);
                for (typename Vector1::iterator p = num.begin(); p != num.end(); ++p) {
                    _ring.assign(*p, _ring.zero);
                }

                return SS_OK;
            }

            // ----- Certifying inconsistency

            std::unique_ptr<BlasMatrix<Field>> Atp_minor_inv = nullptr;
            if ((appearsInconsistent && method.certifyInconsistency)
                || method.singularSolutionType != SingularSolutionType::Random) {
                // take advantage of the (PLUQ)t factorization to compute
                // an inverse to the leading minor of (TAS_P . (A|b) . TAS_Q)
#ifdef RSTIMING
                tFastInvert.start();
#endif

                // @note std::make_unique is only C++14
                Atp_minor_inv = std::unique_ptr<BlasMatrix<Field>>(new BlasMatrix<Field>(_field, rank, rank));

                FFLAS::fassign(_field, rank, rank, tas.factors->getPointer(), tas.factors->getStride(), Atp_minor_inv->getPointer(), Atp_minor_inv->getStride());
                FFPACK::ftrtri (_field, FFLAS::FflasUpper, FFLAS::FflasNonUnit, rank, Atp_minor_inv->getPointer(), Atp_minor_inv->getStride());
                FFPACK::ftrtri (_field, FFLAS::FflasLower, FFLAS::FflasUnit, rank, Atp_minor_inv->getPointer(), Atp_minor_inv->getStride());
                FFPACK::ftrtrm (_field, FFLAS::FflasLeft, FFLAS::FflasNonUnit, rank, Atp_minor_inv->getPointer(), Atp_minor_inv->getStride());

#ifdef RSTIMING
                tFastInvert.stop();
                ttFastInvert += tFastInvert;
#endif
            }

            // ----- Confirm inconsistency if it looks like it

            // If the system appears inconsistent, we either try a new prime,
            // a validate the inconsistency of (A,b).
            if (appearsInconsistent) {
                auto status = solveApparentlyInconsistent(A_check, tas, Atp_minor_inv.get(), rank, method);

                // Failing means we should try a new prime
                if (status == SS_FAILED) continue;
                return status;
            }

            //
            // Starting from here, we know that the system is consistent mod p.
            //

            // Build A_minor and Ap_minor_inv
            BlasMatrix<Ring>* B = nullptr;               // @fixme Make that a std::unique_ptr
            BlasMatrix<Ring>* P = nullptr;               // @fixme Make that a std::unique_ptr
            BlasMatrix<Ring> A_minor(_ring, rank, rank); // -- will have the full rank minor of A
            BlasMatrix<Field>* Ap_minor_inv = nullptr;   // -- will have inverse mod p of A_minor
            makeConditioner(A_minor, Ap_minor_inv, B, P, A_check, tas, Atp_minor_inv.get(), rank, method);

            // Compute newb = (TAS_P.b)[0..(rank-1)]
            BlasVector<Ring> newb(b);
            BMDI.mulin_right(tas.Q, newb);
            newb.resize(rank);

            // ----- Do lifting on sub matrix

            BlasMatrix<Ring> BBA_minor(A_minor);
            LiftingContainer lc(_ring, _field, BBA_minor, *Ap_minor_inv, newb, _prime);

            // ----- Reconstruct rational

            RationalReconstruction<LiftingContainer> re(lc);
            VectorFraction<Ring> resultVF(_ring, rank);
            if (!re.getRational(resultVF.numer, resultVF.denom, 0)) {
                // dirty, but should not be called
                return SS_FAILED;
            }

#ifdef RSTIMING
            ttSystemSolve.update(re, lc);
            tCheckAnswer.start();
#endif

            // ----- Build effective solution from sub matrix

            if (method.singularSolutionType != SingularSolutionType::Random) {
                // short_answer = TAS_Q * short_answer
                resultVF.numer.resize(A.coldim() + 1, _ring.zero);
                BMDI.mulin_left(resultVF.numer, tas.P);
                resultVF.numer.resize(A.coldim());
            }
            else {
                // short_answer = P * short_answer
                BlasVector<Ring> newNumer(_ring, A.coldim());
                BAR.applyV(newNumer, *P, resultVF.numer);
                resultVF.numer = newNumer;
            }

            // ----- Check consistency

            if (method.checkResult) { // check consistency
                BlasVector<Ring> A_times_xnumer(_ring, b.size());
                BAR.applyV(A_times_xnumer, A_check, resultVF.numer);

                Integer tmpi;

                typename Vector2::const_iterator ib = b.begin();
                typename BlasVector<Ring>::iterator iAx = A_times_xnumer.begin();
                int thisrow = 0;
                bool needNewPrime = false;

                for (; !needNewPrime && ib != b.end(); ++iAx, ++ib, ++thisrow)
                    if (!_ring.areEqual(_ring.mul(tmpi, *ib, resultVF.denom), *iAx)) {
                        // should attempt to certify inconsistency now
                        // as in "if [A31 | A32]y != b3" of step (4)
                        needNewPrime = true;
                    }

                if (needNewPrime) {
                    if (Ap_minor_inv != Atp_minor_inv.get()) {
                        delete Ap_minor_inv;
                    }
                    if (method.singularSolutionType == SingularSolutionType::Random) {
                        delete P;
                    }
#ifdef RSTIMING
                    tCheckAnswer.stop();
                    ttCheckAnswer += tCheckAnswer;
#endif
                    continue; // go to start of main loop
                }
            }

#ifdef RSTIMING
            tCheckAnswer.stop();
            ttCheckAnswer += tCheckAnswer;
#endif

            // ----- We have the result values!
            num = resultVF.numer;
            den = resultVF.denom;


            // ----- Checking minimal denominator

            if (method.certifyMinimalDenominator) {
                certifyMinimalDenominator<Vector1>(A, b, tas, *B, A_minor, *Ap_minor_inv, rank);
            }

            if (Ap_minor_inv != Atp_minor_inv.get()) {
                delete Ap_minor_inv;
            }

            if (method.singularSolutionType == SingularSolutionType::Random) {
                delete P;
            }

            // done making certificate, lets blow this popstand
            return SS_OK;
        }

        // All primes were bad
        return SS_FAILED;
    }
}
