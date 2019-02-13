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

#include "linbox/algorithms/cra-domain.h"
#include "linbox/algorithms/rational-cra.h"
#include "linbox/algorithms/rational-cra2.h"
#include "linbox/integer.h"
#include "linbox/solutions/methods.h"
#include "linbox/util/mpicpp.h"
#include "linbox/util/timer.h"
#include <stdlib.h>
#include <utility>
#include <vector>

#include "linbox/randiter/random-prime.h"
#include <unordered_set>

#include "linbox/algorithms/cra-domain-omp.h"

namespace LinBox {

    template <class CRABase>
    struct HybridChineseRemainder {
        typedef typename CRABase::Domain Domain;
        typedef typename CRABase::DomainElement DomainElement;
        using MaskedPrimedGenerator = LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag>;

    protected:
        CRABase _builder;
        Communicator* _pCommunicator;
        double _hadamardBound;       // hadamard bound
        double _workerHadamardBound; // Local Hadamard bound
        bool _builderInitialized = false;
        int _threadsCount = 1;

    public:
        template <class Param>
        HybridChineseRemainder(const Param& b, Communicator* c)
            : _builder(b)
            , _pCommunicator(c)
            , _hadamardBound(b) // Init with hadamard bound
        {
            // @fixme _threadsCount should be available through some parameters...
            #pragma omp parallel
            #pragma omp single
            _threadsCount = omp_get_num_threads();

            // Each (MPI) worker will sum up primes until this is hit
            _workerHadamardBound = _hadamardBound / (_pCommunicator->size() - 1);
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
        template <class Vect, class Function, class PrimeIterator>
        Vect& operator()(Vect& num, Integer& den, Function& Iteration, PrimeIterator& primeg)
        {
            //  if there is no communicator or if there is only one process,
            //  then proceed normally (without parallel)
            if (_pCommunicator == 0 || _pCommunicator->size() == 1) {
                RationalRemainder<CRABase> sequential(_builder);
                return sequential(num, den, Iteration, primeg);
            }

            para_compute(num, Iteration, primeg);

            if (_pCommunicator->rank() == 0) {
                return _builder.result(num, den);
            }
            else {
                return num;
            }
        }


        template <class Function>
        void worker_process_task(Function& Iteration, size_t vectorSize)
        {
            MaskedPrimedGenerator gen(_pCommunicator->rank() - 1, _pCommunicator->size() - 1);
            ++gen;

            //
            // Resource allocation and tasks evaluations
            //

            std::vector<BlasVector<Domain>> residues;
            std::vector<Domain> domains;

            // Pre-reserve with some estimation
            size_t maxIterations = std::ceil(1.442695040889 * _workerHadamardBound / (gen.getBits() - 1));
            residues.reserve(maxIterations);
            domains.reserve(maxIterations);

            double primeLogSize = 0.0;
            while (primeLogSize < _workerHadamardBound) {
                do {
                    ++gen;
                } while (_builder.noncoprime(*gen));
                int prime = *gen;

                primeLogSize += Givaro::logtwo(prime);
                domains.emplace_back(prime);
                residues.emplace_back(domains.back(), vectorSize + 1);
            }

            uint64_t taskCount = domains.size();

            _pCommunicator->send(taskCount, 0);

            //
            // Main parallel loop
            //

            #pragma omp parallel num_threads(_threadsCount)
            #pragma omp single
            {
                // #pragma omp parallel for num_threads(_threadsCount) schedule(dynamic, 1)
                for (uint64_t j = 0; j < taskCount; j++) {
                    #pragma omp task
                    {
                        Iteration(residues[j], domains[j], _pCommunicator);
                        residues[j].push_back(domains[j].characteristic());
                    }
                }
            }

            //
            // Send back result to master
            //

            for (uint64_t j = 0; j < taskCount; j++) {
                _pCommunicator->send(residues[j].begin(), residues[j].end(), 0, 0);
            }

        }
        
        template <class Function>
        void master_process_task(Function& Iteration, size_t vectorSize)
        {
            uint64_t taskCount = 0;
            for (auto workerId = 1; workerId <= _pCommunicator->size() - 1; ++workerId) {
                uint64_t workerTaskCount;
                _pCommunicator->recv(workerTaskCount, workerId);
                taskCount += workerTaskCount;
            }

            std::cout << "Total task count: " << taskCount << std::endl;

            //
            // Main master loop waiting for results
            //

            int pp;
            Domain D(101);
            BlasVector<Domain> residue(D, vectorSize);

            while (taskCount > 0) {

                master_recv_residue(residue, pp, taskCount);

#ifdef __Detailed_Time_Measurement
                Timer chrono;
                chrono.start();
#endif

                Domain D(pp);
                if (!_builderInitialized) {
                    _builder.initialize(D, residue);
                    _builderInitialized = true;
                }
                else {
                    _builder.progress(D, residue);
                }

#ifdef __Detailed_Time_Measurement
                chrono.stop();
                std::cout << "CRT time (seconds): " << chrono.usertime() << std::endl;
#endif
            }
        }
        

        template <class Vect, class Function, class PrimeIterator>
        void para_compute(Vect& num, Function& Iteration, PrimeIterator& primeGenerator)
        {
            if (_pCommunicator->rank() == 0) {
                master_process_task(Iteration, num.size());
            }
            else {
                double starttime = omp_get_wtime();
                worker_process_task(Iteration, num.size());
                double endtime = omp_get_wtime();
                std::cout << " process(" << _pCommunicator->rank() << ") used total CPU time (seconds): " << endtime - starttime << std::endl;
            }
        }

        void master_recv_residue(BlasVector<Domain>& residue, int& pp, size_t& taskCount)
        {
            residue.resize(residue.size() + 1);

            _pCommunicator->recv(residue.begin(), residue.end(), MPI_ANY_SOURCE, 0);

            pp = residue.back();
            residue.resize(residue.size() - 1);
            taskCount -= 1;
        }


    };
}
