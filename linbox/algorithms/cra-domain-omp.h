/* linbox/algorithms/cra-domain-omp.h
 * Copyright (C) 1999-2010 The LinBox group
 *
 * Updated by Hongguang Zhu <zhuhongguang2014@gmail.com>
 *
 *
 * Naive parallel chinese remaindering
 * Launch NN iterations in parallel, where NN=omp_get_num_threads()
 * Then synchronization and termintation test.
 * Time-stamp: <13 Mar 12 13:49:58 Jean-Guillaume.Dumas@imag.fr>
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

/*! @file algorithms/cra-domain-omp.h
 * @brief Parallel (OMP) version of \ref CRA
 * @ingroup CRA
 */

#ifndef __LINBOX_omp_cra_H
#define __LINBOX_omp_cra_H
// commentator is not thread safe
#define DISABLE_COMMENTATOR
#include <omp.h>

#include "linbox/algorithms/cra-domain-seq.h"
#include "linbox/algorithms/rational-cra.h"
#include "linbox/randiter/random-prime.h"
#include <set>

namespace LinBox {

    template <class CRABase>
    struct ChineseRemainderOMP : public ChineseRemainderSeq<CRABase> {
        typedef typename CRABase::Domain Domain;
        typedef typename CRABase::DomainElement DomainElement;
        typedef ChineseRemainderSeq<CRABase> Father_t;

        double B;

        template <class Param>
        ChineseRemainderOMP(const Param& b)
            : Father_t(b)
            , B(b)
        {
        }

        ChineseRemainderOMP(const CRABase& b)
            : Father_t(b)
        {
        }

        template <class Function, class PrimeIterator>
        Integer& operator()(Integer& res, Function& Iteration, PrimeIterator& primeiter)
        {

            int NN;
            PAR_BLOCK { NN = NUM_THREADS; }

            // commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
            if (NN == 1) return Father_t::operator()(res, Iteration, primeiter);

            this->para_compute(Iteration);

            // commentator().stop ("done", NULL, "mmcrait");

            return this->Builder_.result(res);
        }

        template <class Container, class Function, class PrimeIterator>
        Container& operator()(Container& res, Function& Iteration, PrimeIterator& primeiter)
        {

            int NN;
            PAR_BLOCK { NN = NUM_THREADS; }

            // commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
            if (NN == 1) return Father_t::operator()(res, Iteration, primeiter);

            this->para_compute(Iteration);

            // commentator().stop ("done", NULL, "mmcrait");

            return this->Builder_.result(res);
        }

        template <class Function, class PrimeIterator>
        Integer& operator()(Integer& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
        {

            int NN;
            PAR_BLOCK { NN = NUM_THREADS; }

            // commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
            if (NN == 1) return Father_t::operator()(res, den, Iteration, primeiter);

            this->para_compute(Iteration);

            // commentator().stop ("done", NULL, "mmcrait");

            return this->Builder_.result(res, den);

        }

        template <class Function>
        void para_compute(Function Iteration)
        {
            int NN;
            PAR_BLOCK { NN = NUM_THREADS; }


            typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;
            // std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>> m_primeiters;
            std::vector<int> m_primeiters;
//            LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag> gen( 1, NN);;
//            LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::DeterministicTag> gen;
            LinBox::PrimeIterator<LinBox::IteratorCategories::HeuristicTag> gen;
            long Niter = std::ceil(1.442695040889 * B / (double)(gen.getBits() - 1));
            m_primeiters.reserve(Niter);
            std::vector<CRABase> vBuilders;
            vBuilders.reserve(NN);


Givaro::IntPrimeDom _IPD; //!< empty struct dealing with primality.
            for (auto j = 0; j < Niter; j++) {
                ++gen;
                while (this->Builder_.noncoprime(*gen)) ++gen;

if(!_IPD.isprime(*gen)) std::cerr<<*gen<<" is not a prime number ! ********************* "<<std::endl;
std::cout<<*gen<<std::endl;
                m_primeiters.push_back(*gen);
            }


            for (auto j = 0; j < NN; j++) {

                CRABase Builder_(B);
                vBuilders.push_back(Builder_);
            }

            std::vector<ElementContainer> ROUNDresidues;
            ROUNDresidues.resize(Niter);
            std::vector<Domain> ROUNDdomains;
            ROUNDdomains.resize(Niter);

            compute_task((this->Builder_), m_primeiters, Iteration, ROUNDdomains, ROUNDresidues, gen);

        }


        template <class Function, class Domain, class ElementContainer>
        void solve_with_prime(int m_primeiters, Function& Iteration, Domain& ROUNDdomains, ElementContainer& ROUNDresidues, bool& foundInconsistentSystem)
        {

            ROUNDdomains = Domain(m_primeiters);

            Iteration(ROUNDresidues, ROUNDdomains, foundInconsistentSystem);
        }



        template <class pFunc, class Function, class Domain, class ElementContainer, class PrimeGen>
        void compute_task(pFunc& pF, std::vector<int>& m_primeiters, Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, PrimeGen& gen)
        {

            // LinBox::PrimeIterator<LinBox::IteratorCategories::DeterministicTag> gen;
            long Niter = std::ceil(1.442695040889 * B / (double)(gen.getBits() - 1));


            std::vector<int> failedPrimesIndices; // Failed due to inconsistent systems mod p

            //            FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Grain> H;
            splitting_3(H, FFLAS::CuttingStrategy::Block, FFLAS::StrategyParameter::Grain);

            PAR_BLOCK
            {

                FOR1D(i, Niter, H, {
                    // std::cout<<"i = "<<i<<" == "<<"  _internal_iterator.blockindex()  : "<<"_internal_iterator.blockindex()%NN
                    // = "<<_internal_iterator.blockindex()%NN<<std::endl;
                    //    std::cout << "Threads: " << NUM_THREADS << ", max: " << MAX_THREADS << std::endl;
                    //    std::cerr << "OMP: " << __FFLASFFPACK_USE_OPENMP << ", max " << omp_get_max_threads() << std::endl;
                    bool foundInconsistentSystem = false;
                    solve_with_prime(m_primeiters[i], Iteration, ROUNDdomains[i], ROUNDresidues[i], foundInconsistentSystem);

                    if (foundInconsistentSystem) {
                        std::cout << "foundInconsistentSystem" << std::endl;
                        failedPrimesIndices.push_back(i);
                    }
                })
            }

            for (auto i : failedPrimesIndices) {
                ++gen;
                while(this->Builder_.noncoprime(*gen) )
                    ++gen;
                m_primeiters[i] = *gen;

                bool foundInconsistentSystem = false;
                solve_with_prime(m_primeiters[i], Iteration, ROUNDdomains[i], ROUNDresidues[i], foundInconsistentSystem);
                if (foundInconsistentSystem) {
                    // Couldn't solve it the second time, we crash
                    throw LinboxMathInconsistentSystem("Solve Paladin.");
                }
            }

            this->Builder_.initialize(ROUNDdomains[0], ROUNDresidues[0]);

            for (auto j = 1; j < Niter; j++) {

                this->Builder_.progress(ROUNDdomains[j], ROUNDresidues[j]);
            }
        }


        template <class Container, class Function, class PrimeIterator>
        Container& operator()(Container& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
        {

            int NN;
            PAR_BLOCK { NN = NUM_THREADS; }
            // commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
            if (NN == 1) return Father_t::operator()(res, den, Iteration, primeiter);

            this->para_compute(Iteration);

            // commentator().stop ("done", NULL, "mmcrait");

             return this->Builder_.result(res,den);

        }
    };
}

#endif //__LINBOX_omp_cra_H

// Local Variables:
// mode: C++
// tab-width: 4
// indent-tabs-mode: nil
// c-basic-offset: 4
// End:
// vim:sts=4:sw=4:ts=4:et:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
