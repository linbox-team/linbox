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


#include <set>
#include "linbox/algorithms/cra-domain-seq.h"
#include "linbox/randiter/random-prime.h"
#include "linbox/algorithms/rational-cra.h"

namespace LinBox
{
    
	template<class CRABase>
	struct ChineseRemainderOMP : public ChineseRemainderSeq<CRABase> {
		typedef typename CRABase::Domain	Domain;
		typedef typename CRABase::DomainElement	DomainElement;
		typedef ChineseRemainderSeq<CRABase>    Father_t;
        
        double B;
        
		template<class Param>
		ChineseRemainderOMP(const Param& b) :
			Father_t(b), B(b)
		{}
        
		ChineseRemainderOMP(const CRABase& b) :
			Father_t(b)
		{}
        
		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Function& Iteration, PrimeIterator& primeiter)
		{
            
			size_t NN;
#pragma omp parallel
#pragma omp single
			NN = omp_get_num_threads();

			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);
            
            this->para_compute(Iteration);
            
			// commentator().stop ("done", NULL, "mmcrait");
            
			return this->Builder_.result(res);
		}
        
		template<class Container, class Function, class PrimeIterator>
		Container& operator() (Container& res, Function& Iteration, PrimeIterator& primeiter)
		{
            
			size_t NN;
#pragma omp parallel
#pragma omp single
			NN = omp_get_num_threads();

			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);
            
            this->para_compute(Iteration);
            
			// commentator().stop ("done", NULL, "mmcrait");
            
			return this->Builder_.result(res);
		}
        
        
		template<class Function, class PrimeIterator>
		Integer& operator() (Integer& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{
			//! @bug why why why ???
			/** erreur: ‘omp_get_num_threads’ has not been declared
			 * ../linbox/algorithms/cra-domain-omp.h:152:16: note: suggested alternative:
			 * /usr/lib/gcc/x86_64-linux-gnu/4.6/include/omp.h:64:12: note:   ‘Givaro::omp_get_num_threads’
			 */
			size_t NN;
#pragma omp parallel
			NN = omp_get_num_threads();

			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res, den, Iteration,primeiter);
            
            this->para_compute(Iteration);			
			
			// commentator().stop ("done", NULL, "mmcrait");
            
			return this->Builder_.result(res,den);
		}
		

        template<class Function>
        void para_compute(Function Iteration){
			size_t NN;
#pragma omp parallel
#pragma omp single
			NN = omp_get_num_threads();


            typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;

            std::vector<ElementContainer> ROUNDresidues;ROUNDresidues.resize(NN);
            std::vector<Domain> ROUNDdomains;ROUNDdomains.resize(NN);
            std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>> m_primeiters;m_primeiters.reserve(NN);

            std::vector<CRABase> vBuilders;vBuilders.reserve(NN);


            for(auto j=0u;j<NN;j++){
                LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag> m_primeiter( j, NN);
                ++m_primeiter; m_primeiters.push_back(m_primeiter);

                CRABase Builder_(B);vBuilders.push_back(Builder_);

            }
           
		compute_task( (this->Builder_), m_primeiters, Iteration,  ROUNDdomains,
                          ROUNDresidues, vBuilders);


        }


        template< class Function, class PrimeIterator, class Domain, class ElementContainer>
        void solve_with_prime(std::vector<PrimeIterator>& m_primeiters,
                              Function& Iteration, std::vector<Domain>& ROUNDdomains,
                              std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders)
        {

            ++m_primeiters[ omp_get_thread_num()];

            while(vBuilders[ omp_get_thread_num()].noncoprime(*m_primeiters[ omp_get_thread_num()]) )
                ++m_primeiters[ omp_get_thread_num()];


            ROUNDdomains[ omp_get_thread_num()] = Domain(*m_primeiters[ omp_get_thread_num()]);

            Iteration(ROUNDresidues[ omp_get_thread_num()], ROUNDdomains[ omp_get_thread_num()]);

        }


        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<PrimeIterator>& m_primeiters,
                          Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders)
        {

			size_t NN;
#pragma omp parallel
#pragma omp single
			NN = omp_get_num_threads();
            long Niter=std::ceil(1.442695040889*B/(double)(m_primeiters[0].getBits()-1));
            
            bool need2Init=true;

#pragma omp parallel for num_threads(NN) schedule(dynamic,1)
            for(auto j=0;j<Niter;j++)
                {

                    solve_with_prime(m_primeiters, Iteration, ROUNDdomains, ROUNDresidues, vBuilders);

#pragma omp critical
                    if(need2Init){
                        this->Builder_.initialize( ROUNDdomains[ omp_get_thread_num()], ROUNDresidues[ omp_get_thread_num()]);
                        need2Init=false;
                    }else{
                        this->Builder_.progress( ROUNDdomains[ omp_get_thread_num()], ROUNDresidues[ omp_get_thread_num()]);    

                    }

                }
        }

		template<class Container, class Function, class PrimeIterator>
		Container& operator()  (Container& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{

			size_t NN;
#pragma omp parallel
#pragma omp single
			NN = omp_get_num_threads();
			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res, den, Iteration,primeiter);

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
