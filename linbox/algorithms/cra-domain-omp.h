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
			//! @bug why why why ???
			/** erreur: ‘omp_get_num_threads’ has not been declared
			 * ../linbox/algorithms/cra-domain-omp.h:152:16: note: suggested alternative:
			 * /usr/lib/gcc/x86_64-linux-gnu/4.6/include/omp.h:64:12: note:   ‘Givaro::omp_get_num_threads’
			 */
            
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 

			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res,Iteration,primeiter);
            
            this->para_compute(Iteration);
            
			// commentator().stop ("done", NULL, "mmcrait");
            
			return this->Builder_.result(res);
		}
        
		template<class Container, class Function, class PrimeIterator>
		Container& operator() (Container& res, Function& Iteration, PrimeIterator& primeiter)
		{
            
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 

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
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 

			// commentator().start ("Parallel OMP Givaro::Modular iteration", "mmcrait");
			if (NN == 1) return Father_t::operator()(res, den, Iteration,primeiter);
            
            this->para_compute(Iteration);			
			
			// commentator().stop ("done", NULL, "mmcrait");
            
			return this->Builder_.result(res,den);
		}
		
        template<class Function>
        void para_compute(Function Iteration){

            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; }

            typedef typename CRATemporaryVectorTrait<Function, Domain>::Type_t ElementContainer;
            std::vector<LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag>> m_primeiters;
            m_primeiters.reserve(NN);
            std::vector<CRABase> vBuilders;
            vBuilders.reserve(NN);


            for(auto j=0;j<NN;j++){
                LinBox::MaskedPrimeIterator<LinBox::IteratorCategories::HeuristicTag> m_primeiter( j, NN);
                m_primeiters.push_back(m_primeiter);
                
                CRABase Builder_(B);
                vBuilders.push_back(Builder_);
                
            }


            long Niter=std::ceil(1.442695040889*B/(double)(m_primeiters[0].getBits()-1));
            
            std::vector<ElementContainer> ROUNDresidues;ROUNDresidues.resize(Niter);
            std::vector<Domain> ROUNDdomains;ROUNDdomains.resize(Niter);

            compute_task( (this->Builder_), m_primeiters, Iteration,  ROUNDdomains,
                          ROUNDresidues, vBuilders);
            
            
        }


        template< class Function, class PrimeIterator, class Domain, class ElementContainer>
        void solve_with_prime(PrimeIterator& m_primeiters,
                              Function& Iteration, Domain& ROUNDdomains,
                              ElementContainer& ROUNDresidues, CRABase& vBuilders)
        {
            
            ++m_primeiters;

            while(vBuilders.noncoprime(*m_primeiters) )
                ++m_primeiters;

            
            ROUNDdomains = Domain(*m_primeiters);
            
            Iteration(ROUNDresidues, ROUNDdomains);
            
        }

#if 0 //Paladin impl
        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<PrimeIterator>& m_primeiters,
                          Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders)
        {
            
            long Niter=std::ceil(1.442695040889*B/(double)(m_primeiters[0].getBits()-1));
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 
            
            FFLAS::ParSeqHelper::Parallel<FFLAS::CuttingStrategy::Single,FFLAS::StrategyParameter::Grain> H;
            if(NN>Niter){

	            PARFORBLOCK1D(k,Niter,H,
                    solve_with_prime(m_primeiters[k], Iteration, ROUNDdomains[k], ROUNDresidues[k], vBuilders[k]);
                );

            }else{

	            PARFORBLOCK1D(k,NN,H,
	                for(auto j=k*(Niter/NN);j<(k+1)*(Niter/NN);j++)
                        solve_with_prime(m_primeiters[k], Iteration, ROUNDdomains[j], ROUNDresidues[j], vBuilders[k]);
                );
                if(Niter%NN>0){
	                PARFORBLOCK1D(k,Niter%NN,H,
                        solve_with_prime(m_primeiters[k], Iteration, ROUNDdomains[k+NN*(Niter/NN)], ROUNDresidues[k+NN*(Niter/NN)], vBuilders[k]);
                    );
                }

            }


            this->Builder_.initialize( ROUNDdomains[0], ROUNDresidues[0]);

            for(auto j=0;j<Niter;j++)
                {

                        this->Builder_.progress( ROUNDdomains[j], ROUNDresidues[j]);

                }
           
        }
#else //OMP impl
        template<class pFunc, class Function, class PrimeIterator, class Domain, class ElementContainer>
        void compute_task(pFunc& pF, std::vector<PrimeIterator>& m_primeiters,
                          Function& Iteration, std::vector<Domain>& ROUNDdomains,
                          std::vector<ElementContainer>& ROUNDresidues, std::vector<CRABase>& vBuilders)
        {
            
            long Niter=std::ceil(1.442695040889*B/(double)(m_primeiters[0].getBits()-1));
            int NN;
#pragma omp parallel
#pragma omp single 
            NN=NUM_THREADS;
            
#pragma omp parallel for schedule(dynamic,1) num_threads(NN)
	        for(auto j=0;j<Niter;j++){
                solve_with_prime(m_primeiters[omp_get_thread_num()], Iteration, ROUNDdomains[j], ROUNDresidues[j], vBuilders[omp_get_thread_num()]);
                }


            this->Builder_.initialize( ROUNDdomains[0], ROUNDresidues[0]);

            for(auto j=0;j<Niter;j++)
                {

                        this->Builder_.progress( ROUNDdomains[j], ROUNDresidues[j]);

                }
           
        }

#endif
   
        
		template<class Container, class Function, class PrimeIterator>
		Container& operator()  (Container& res, Integer& den, Function& Iteration, PrimeIterator& primeiter)
		{
            
            int NN;
            PAR_BLOCK{ NN=NUM_THREADS; } 
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
