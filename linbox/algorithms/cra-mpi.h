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
#include "linbox/algorithms/rational-cra2.h"
#include "linbox/integer.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/mpicpp.h"
#include "linbox/util/timer.h"

namespace LinBox {

    template <class CRABase>
    struct MPIChineseRemainder {
        using Domain = typename CRABase::Domain;

        // @note This causes the algorithm to give wrong answer sometimes...
        // using MaskedPrimeGenerator = MaskedPrimeIterator<IteratorCategories::HeuristicTag>;
        using MaskedPrimeGenerator = MaskedPrimeIterator<IteratorCategories::DeterministicTag>;

    protected:
        CRABase Builder_;
        Communicator* _communicator;
        double _hadamardLogBound;

    public:
        MPIChineseRemainder(double b, Communicator* c)
            : Builder_(b)
            , _communicator(c)
            , _hadamardLogBound(b)
        {
        }

        int iterationCount()
        {
            auto bits = MaskedPrimeGenerator(0, _communicator->size()).getBits();
            return std::ceil(_hadamardLogBound / (double)(bits - 1));
        }

        /** \brief The CRA loop.
         *
         * termination condition.
         *
         * \param Iteration  Function object of two arguments, \c
         * Iteration(r, p), given prime \c p it outputs residue(s) \c
         * r.  This loop may be parallelized.  \p Iteration must be
         * reentrant, thread safe.  For example, \p Iteration may be
         * returning the coefficients of the minimal polynomial of a
         * matrix \c mod \p p.
         @warning  we won't detect bad primes.
         *
         * \param primeGenerator  RandIter object for generating primes.
         * \param[out] res an integer
         */
        template <class Vect, class Function, class PrimeIterator>
        Vect& operator()(Vect& num, Integer& den, Function& Iteration, PrimeIterator& primeGenerator)
        {
            // Defer to standard CRA loop if no parallel usage is desired
            if (_communicator == 0 || _communicator->size() == 1) {
                RationalRemainder<CRABase> sequential(Builder_);
                return sequential(num, den, Iteration, primeGenerator);
            }

            para_compute(num, Iteration, primeGenerator);

            if (_communicator->master()) {
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
            if (_communicator == 0 || _communicator->size() == 1) {
                ChineseRemainder<CRABase> sequential(Builder_);
                return sequential(res, Iteration, primeGenerator);
            }

            para_compute(res, Iteration, primeGenerator);

            if (_communicator->master()) {
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

            if (_communicator->master()) {
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

            if (_communicator->master()) {
                master_process_task(Iteration, D, r);
            }
            else {
                worker_process_task(Iteration, r);
            }
        }

        template <class Vect, class PrimeIterator, class Function>
        void worker_compute(PrimeIterator& gen, Function& Iteration, Vect& r)
        {
            // Process mutual independent prime number generation
            ++gen;
            while (Builder_.noncoprime(*gen)) {
                ++gen;
            }

            Domain D(*gen);
            Iteration(r, D);
        }

        template <class Vect, class Function>
        void worker_process_task(Function& Iteration, Vect& r)
        {
            // Identify how many tasks we need
            int iterations = iterationCount();
            int procs = _communicator->size() - 1; // We don't count master

            int Ntask = 0;
            if ((iterations % procs) >= _communicator->rank()) {
                Ntask = (int)std::ceil((double)iterations / (double)procs);
            }
            else {
                Ntask = (int)std::floor((double)iterations / (double)procs);
            }

            // Ok, just go compute them
            MaskedPrimeGenerator gen(_communicator->rank(), _communicator->size());

            for (long i = 0; i < Ntask; i++) {
                worker_compute(gen, Iteration, r);

                uint64_t p = *gen;
                _communicator->send(p, 0);
                _communicator->send(r, 0);
            }
        }

        template <class Vect, class Function>
        void master_process_task(Function& Iteration, Domain& D, Vect& r)
        {
            Iteration(r, D);
            Builder_.initialize(D, r);

            int Nrecv = iterationCount();
            while (Nrecv > 0) {
                // Receive the prime and a result vector
                uint64_t p;
                _communicator->recv(p, MPI_ANY_SOURCE);
                _communicator->recv(r, _communicator->status().MPI_SOURCE);

                // Update the number of iterations for the next step
                Nrecv--;

                Domain D(p);
                Builder_.progress(D, r);
            }
        }
    };
}
