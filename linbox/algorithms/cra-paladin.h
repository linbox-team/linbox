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
//#define DISABLE_COMMENTATOR
#include <omp.h>

#include <set>
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/rational-cra.h"

namespace LinBox
{

	template<class CRABase>
	struct RationalChineseRemainderPaladin {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
		typedef RationalChineseRemainder<CRABase>    Father_t;

		CRABase Builder_;
        double HB;//hadamard bound
        int Niter;

		template<class Param>
		RationalChineseRemainderPaladin(const Param& b) :
			Builder_(b), HB(b), Niter(0)
		{}

		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter)
		{

            this->para_compute<Integer>(Iteration,primeiter);
			// commentator().stop ("done", NULL, "mmcrait");

			return this->Builder_.result(res);
		}

        // @fixme REMOVE ?
		template<class Container, class Function, class PrimeIterator>
		Container& operator() (Container& res, Function& Iteration, PrimeIterator& primeiter)
		{

            this->para_compute<Container>(Iteration,primeiter);

			// commentator().stop ("done", NULL, "mmcrait");

			return this->Builder_.result(res);
		}


		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{

            this->para_compute<Integer>(Iteration,primeiter);

			// commentator().stop ("done", NULL, "mmcrait");

			return this->Builder_.result(res,den);
		}


		template<class Vect, class Function, class PrimeIterator>
		Vect& operator() (Vect& num, Vect& den, Function& Iteration, PrimeIterator& primeiter)
		{

            this->para_compute<BlasVector<Domain>>(Iteration,primeiter); // @fixme Should be rebind of Vect into Domain

			// commentator().stop ("done", NULL, "mmcrait");

			return this->Builder_.result(num,den);
		}


        template<class ResidueType, class Function, class PrimeIterator>
        void para_compute(Function &Iteration, PrimeIterator& primeiter){
            this->Niter = std::ceil(1.442695040889*HB/(double)(primeiter.getBits()-1));
            std::vector<ResidueType> ROUNDresidues;
            ROUNDresidues.resize(this->Niter);
            LinBox::PrimeIterator<LinBox::IteratorCategories::DeterministicTag> gen(primeiter.getBits());
            std::vector<int> m_primeiters;m_primeiters.reserve(this->Niter);

            for(auto j=0;j<Niter;j++){
                    ++gen;
                    while(this->Builder_.noncoprime(*gen) )
                        ++primeiter;
                    m_primeiters.push_back(*gen);
            }

		    this->compute_task(m_primeiters, Iteration, ROUNDresidues);


        }


        template< class Function, class ResidueType>
        void solve_with_prime(int m_primeiter, Function& Iteration, ResidueType& ROUNDresidue)
        {

            Domain D(m_primeiter);
            Iteration(ROUNDresidue, D);
        }


        template<class Function, class PrimeIterator, class ResidueType>
        void compute_task(std::vector<PrimeIterator>& m_primeiters,
                          Function& Iteration, std::vector<ResidueType>& ROUNDresidues)
        {

            long NN = this->Niter;
#if 0
            PAR_BLOCK{
                        auto sp=SPLITTER(NUM_THREADS,FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Threads);
                        SYNCH_GROUP({
                            FORBLOCK1D(iter, NN, sp,{
                                TASK(MODE(CONSTREFERENCE(m_primeiters,Iteration,ROUNDresidues)),{
                                    for(auto j=iter.begin(); j!=iter.end(); ++j)
                                    {

                                        solve_with_prime(m_primeiters[j], Iteration, ROUNDresidues[j]);
                                    }
                                 })
                            })

                         })
            }
#else
            PAR_BLOCK{
                        auto sp=SPLITTER(NUM_THREADS,FFLAS::CuttingStrategy::Row,FFLAS::StrategyParameter::Threads);
                            FOR1D(iter, NN, sp, MODE(CONSTREFERENCE(m_primeiters,Iteration,ROUNDresidues)),
                            {
                                        solve_with_prime(m_primeiters[iter], Iteration, ROUNDresidues[iter]);
                            })

            }
#endif
            Domain D(m_primeiters[0]);
            this->Builder_.initialize( D, ROUNDresidues[0]);
            for(auto j=1;j<NN;j++){
                Domain D(m_primeiters[j]);
                this->Builder_.progress( D, ROUNDresidues[j]);
            }
        }


		template<class Container, class Function, class PrimeIterator>
		Container& operator()  (Container& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{

		        this->para_compute<BlasVector<Domain>>(Iteration,primeiter); // @fixme SHould be container -> rebind

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
