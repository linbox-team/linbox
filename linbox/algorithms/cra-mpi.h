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

#include <stdlib.h>
#include <unordered_set>
#include <utility>
#include <vector>

#include "linbox/algorithms/cra-domain-omp.h"
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
        typedef typename CRABase::Domain Domain;
        typedef typename CRABase::DomainElement DomainElement;

    protected:
        CRABase Builder_;
        Communicator* _commPtr;
        unsigned int _numprocs;
        double HB; // hadamard bound

    public:
        template <class Param>
        MPIChineseRemainder(const Param& b, Communicator* c)
            : Builder_(b)
            , _commPtr(c)
            , _numprocs(c->size())
            , HB(b) // Init with hadamard bound
        {
        }

        int getNiter()
        {
            auto bits = LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>(0, _commPtr->size()).getBits();
            return std::ceil(1.442695040889 * HB / (double)(bits - 1));
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
         * \param primeg  RandIter object for generating primes.
         * \param[out] res an integer
         */
        template <class Function, class PrimeIterator>
        Integer& operator()(Integer& res, Function& Iteration, PrimeIterator& primeg)
        {
            //  defer to standard CRA loop if no parallel usage is desired
            if (_commPtr == 0 || _commPtr->size() == 1) {
                ChineseRemainder<CRABase> sequential(Builder_);
                return sequential(res, Iteration, primeg);
            }

            para_compute(res, Iteration, primeg);
            if (_commPtr->rank() == 0) {
                return Builder_.result(res);
            }
            else {
                return res;
            }
        }

        template <class Function, class PrimeIterator>
        Integer& operator()(Integer& num, Integer& den, Function& Iteration, PrimeIterator& primeg)
        {

            //  defer to standard CRA loop if no parallel usage is desired
            if (_commPtr == 0 || _commPtr->size() == 1) {
                RationalRemainder<CRABase> sequential(Builder_);
                return sequential(num, den, Iteration, primeg);
            }
            para_compute(num, Iteration, primeg);
            if (_commPtr->rank() == 0) {
                return Builder_.result(num, den);
            }
            else {
                return num;
            }
        }

        template <class Vect, class Function, class PrimeIterator>
        Vect& operator()(Vect& num, Integer& den, Function& Iteration, PrimeIterator& primeg)
        {
            //  if there is no communicator or if there is only one process,
            //  then proceed normally (without parallel)
            if (_commPtr == 0 || _commPtr->size() == 1) {

                RationalRemainder<CRABase> sequential(Builder_);
                return sequential(num, den, Iteration, primeg);
            }
            para_compute(num, Iteration, primeg);
            if (_commPtr->rank() == 0) {
                return Builder_.result(num, den);
            }
            else {
                return num;
            }
        }

        template <class Vect, class Function, class PrimeIterator>
        Vect& operator()(Vect& num, Function& Iteration, PrimeIterator& primeg)
        {
            //  if there is no communicator or if there is only one process,
            //  then proceed normally (without parallel)
            if (_commPtr == 0 || _commPtr->size() == 1) {

                ChineseRemainder<CRABase> sequential(Builder_);
                return sequential(num, Iteration, primeg);
            }
            para_compute(num, Iteration, primeg);
            if (_commPtr->rank() == 0) {
                return Builder_.result(num);
            }
            else {
                return num;
            }
        }

        template <class Vect, class Function, class PrimeIterator>
        void para_compute(Vect& num, Function& Iteration, PrimeIterator& primeg)
        {

            Domain D(*primeg);
            BlasVector<Domain> r(D);
            Timer chrono;

            //  parent propcess
            if (_commPtr->rank() == 0) {

                master_process_task(Iteration, D, r);
            }
            //  child process
            else {

                worker_process_task(Iteration, r);
            }
        }

        template <class Vect, class PrimeIterator, class Function>
        void worker_compute(std::unordered_set<int>& prime_used, PrimeIterator& gen, Function& Iteration, Vect& r)
        {
            // Process mutual independent prime number generation
            ++gen;
            while (Builder_.noncoprime(*gen)) ++gen;
            prime_used.insert(*gen);
            Domain D(*gen);
            Iteration(r, D);
        }

        template <class Vect, class Function>
        void worker_process_task(Function& Iteration, Vect& r)
        {

            int Ntask = 0;
            LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag> gen(_commPtr->rank(), _commPtr->size());
            // LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag>   gen(_commPtr->rank(),_commPtr->size());

            _commPtr->recv(Ntask, 0);

            if (Ntask != 0) {
                std::unordered_set<int> prime_used;

                for (long i = 0; i < Ntask; i++) {
                    worker_compute(prime_used, gen, Iteration, r);

                    // Add corresponding prime number as the last element in the result vector
                    r.push_back(*gen);

                    _commPtr->send(r.begin(), r.end(), 0, 0);
                }
            };
        }

        template <class Vect>
        void compute_state_comm(int* vTaskDist, Vect& r, int& recvprime, int& Nrecv)
        {

            r.resize(r.size() + 1);

            // receive the beginnin and end of a vector in heapspace
            _commPtr->recv(r.begin(), r.end(), MPI_ANY_SOURCE, 0);

            // Update the number of iterations for the next step
            Nrecv--;

            // Store the corresponding prime number
            recvprime = r[r.size() - 1];

            // Restructure the vector without added prime number
            r.resize(r.size() - 1);
        }
        template <class Vect>
        void master_compute(int* vTaskDist, Vect& r, int Niter)
        {

            int Nrecv = Niter;
            int recvprime;

            while (Nrecv > 0) {

                compute_state_comm(vTaskDist, r, recvprime, Nrecv);

                Domain D(recvprime);

                Builder_.progress(D, r);
            }
        }

        template <class Vect, class Function>
        void master_init(int* vTaskDist, Function& Iteration, Domain& D, Vect& r)
        {
            int procs = _commPtr->size();

            int Niter = getNiter();

            // Compute nb of tasks ought to be realized for each process

            if (Niter < (procs - 1)) {

                for (long i = 1; i < Niter + 1; i++) {
                    vTaskDist[i - 1] = 1;
                    _commPtr->send(vTaskDist[i - 1], i);
                }
                for (long i = Niter + 1; i < procs; i++) {
                    vTaskDist[i - 1] = 0;
                    _commPtr->send(vTaskDist[i - 1], i);
                }
            }
            else {
                for (long i = 1; i < Niter % (procs - 1) + 1; i++) {
                    vTaskDist[i - 1] = Niter / (procs - 1) + 1;
                    _commPtr->send(vTaskDist[i - 1], i);
                }
                for (long i = Niter % (procs - 1) + 1; i < procs; i++) {
                    vTaskDist[i - 1] = Niter / (procs - 1);
                    _commPtr->send(vTaskDist[i - 1], i);
                }
            }

            // Initialize the buider and the receiver vector r
            Builder_.initialize(D, Iteration(r, D));
        }

        template <class Vect, class Function>
        void master_process_task(Function& Iteration, Domain& D, Vect& r)
        {
            int vTaskDist[_commPtr->size() - 1];
            int Niter = getNiter();
            master_init(vTaskDist, Iteration, D, r);

            master_compute(vTaskDist, r, Niter);
        }
    };
}
