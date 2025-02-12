/* Copyright (C) 2007 LinBox
 * Updated by Hongguang ZHU
 * Written by bds and zw
 * author: B. David Saunders and Zhendong Wan
 * parallelized for BOINC computing by Bryan Youse
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

#pragma once

#include <unordered_set>
#include <utility>
#include <vector>

#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/rational-cra.h"
#include "linbox/algorithms/rational-cra-var-prec.h"
#include "linbox/integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/mpicpp.h"
#include "linbox/util/timer.h"

#if defined(__LINBOX_HAVE_MPI)

namespace LinBox {

    template <class CRABase>
    struct ChineseRemainderDistributed {
        using Domain = typename CRABase::Domain;

        // @note This causes the algorithm to give wrong answer sometimes...
        // using MaskedPrimeGenerator = MaskedPrimeIterator<IteratorCategories::HeuristicTag>;
        using MaskedPrimeGenerator = MaskedPrimeIterator<IteratorCategories::DeterministicTag>;

    protected:
        CRABase Builder_;
        Communicator* _pCommunicator;
        double _hadamardLogBound;
        double _workerHadamardLogBound = 0.0; //!< Each worker will compute primes until this is hit.

    public:
        ChineseRemainderDistributed(double b, Communicator* c)
            : Builder_(b)
            , _pCommunicator(c)
            , _hadamardLogBound(b)
        {
            if (c && c->size() > 1) {
                _workerHadamardLogBound = _hadamardLogBound / (c->size() - 1);
            }
        }

        /** \brief The CRA loop.
         *
         * \param Iteration  Function object of two arguments, \c
         * Iteration(r, p), given prime \c p it outputs residue(s) \c
         * r.  This loop may be parallelized.  \p Iteration must be
         * reentrant, thread safe.  For example, \p Iteration may be
         * returning the coefficients of the minimal polynomial of a
         * matrix \c mod \p p.
         *
         * @warning  we won't detect bad primes.
         *
         * \param primeGenerator  RandIter object for generating primes.
         * \param[out] res an integer
         */
        template <class Vect, class Function, class PrimeIterator>
        Vect& operator()(Vect& num, Integer& den, Function& Iteration, PrimeIterator& primeGenerator)
        {
            // Defer to standard CRA loop if no parallel usage is desired
            if (_pCommunicator == 0 || _pCommunicator->size() == 1) {
                RationalChineseRemainder<CRABase> sequential(Builder_);
                return sequential(num, den, Iteration, primeGenerator);
            }

            para_compute(num, Iteration, primeGenerator);

            if (_pCommunicator->master()) {
                return Builder_.result(num, den);
            }
            else {
                return num;
            }
        }

        template <class Any, class Function, class PrimeIterator>
        Any& operator()(Any& res, Function& Iteration, PrimeIterator& primeGenerator)
        {
            // Defer to standard CRA loop if no parallel usage is desired
            if (_pCommunicator == 0 || _pCommunicator->size() == 1) {
                ChineseRemainder<CRABase> sequential(Builder_);
                return sequential(res, Iteration, primeGenerator);
            }

            para_compute(res, Iteration, primeGenerator);

            if (_pCommunicator->master()) {
                return Builder_.result(res);
            }
            else {
                return res;
            }
        }

        template <class Any, class Function, class PrimeIterator>
        void para_compute(Any& res, Function& Iteration, PrimeIterator& primeGenerator) {
            Domain D(*primeGenerator);
            typename Domain::Element r;

            if (_pCommunicator->master()) {
                master_process_task(Iteration, D, r);
            }
            else {
                worker_process_task(Iteration, r);
            }
        }

        template <class Ring, class Function, class PrimeIterator>
        void para_compute(BlasVector<Ring>& num, Function& Iteration, PrimeIterator& primeGenerator)
        {
            Domain D(*primeGenerator);
            BlasVector<Domain> r(D);

            if (_pCommunicator->master()) {
                master_process_task(Iteration, D, r);
            }
            else {
                worker_process_task(Iteration, r);
            }
        }

        template <class Any, class PrimeIterator, class Function>
        void worker_compute(PrimeIterator& gen, Function& Iteration, Any& r)
        {
            // Process mutual independent prime number generation
            ++gen;
            while (Builder_.noncoprime(*gen)) {
                ++gen;
            }

            Domain D(*gen);
            Iteration(r, D);
        }

        template <class Any, class Function>
        void worker_process_task(Function& Iteration, Any& r)
        {
            MaskedPrimeGenerator gen(_pCommunicator->rank() - 1, _pCommunicator->size() - 1);

            // Each worker will work until _workerHadamardLogBound is hit
            double primesLogSum = 0.0;
            while (primesLogSum < _workerHadamardLogBound) {
                worker_compute(gen, Iteration, r);

                uint64_t p = *gen;
                primesLogSum += Givaro::logtwo(p);
                _pCommunicator->send(p, 0);
                _pCommunicator->send(r, 0);
            }

            uint64_t poisonPill = 0;
            _pCommunicator->send(poisonPill, 0);
        }

        template <class Any, class Function>
        void master_process_task(Function& Iteration, Domain& D, Any& r)
        {
            Iteration(r, D);
            Builder_.initialize(D, r);

            uint32_t workersDone = _pCommunicator->size() - 1;
            while (workersDone > 0) {
                // Receive the prime
                uint64_t p;
                _pCommunicator->recv(p, MPI_ANY_SOURCE);
                if (p == 0) {
                    workersDone -= 1;
                    continue;
                }

                // Receive result vector and update builder
                _pCommunicator->recv(r, _pCommunicator->status().MPI_SOURCE);

                Domain D(p);
                Builder_.progress(D, r);
            }
        }
    };
}

#endif

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
